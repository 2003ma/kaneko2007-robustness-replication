#!/usr/bin/env python3

"""Plot worst fitness vs generation from data,sigma*.csv files.

Example:
    python3 plot_worst_by_sigma.py \
        --input-dir ../../evo_sim/results/data_sigma_only/ver1 \
        --sigmas 0.005 0.1 0.4 0.6 1.0 \
        --output worst_overlay_bw.pdf
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd


SIGMA_FILE_PATTERN = re.compile(r"^data,sigma=([0-9]*\.?[0-9]+),.*\.csv$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay plot: x=generation, y=worst, line style by sigma"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory that contains data,sigma*.csv (searched recursively)",
    )
    parser.add_argument(
        "--sigmas",
        type=float,
        nargs="*",
        default=None,
        help=(
            "Sigmas to plot, e.g. --sigmas 0.005 0.01 0.04. "
            "If omitted, all sigmas are used."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("worst_vs_generation_overlay.png"),
        help="Output image path",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="",
        help="Plot title. If empty, no title is shown.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Optional max x limit",
    )
    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="Optional min y limit",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional max y limit",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image DPI",
    )
    parser.add_argument(
        "--markevery",
        type=int,
        default=20,
        help="Plot marker every N data points.",
    )
    parser.add_argument(
        "--legend-outside",
        action="store_true",
        default=True,
        help="Place legend outside the plot area.",
    )
    return parser.parse_args()


def extract_sigma_from_name(filename: str) -> float | None:
    m = SIGMA_FILE_PATTERN.match(filename)
    if not m:
        return None
    return float(m.group(1))


def collect_data_files(input_dir: Path) -> List[Tuple[float, Path]]:
    pairs: List[Tuple[float, Path]] = []

    for file_path in sorted(input_dir.rglob("data,sigma*.csv")):
        sigma = extract_sigma_from_name(file_path.name)
        if sigma is not None:
            pairs.append((sigma, file_path))

    return pairs


def deduplicate_by_sigma(
    pairs: Sequence[Tuple[float, Path]],
) -> Tuple[List[Tuple[float, Path]], Dict[float, List[Path]]]:
    by_sigma: Dict[float, List[Path]] = {}

    for sigma, path in pairs:
        by_sigma.setdefault(sigma, []).append(path)

    unique: List[Tuple[float, Path]] = []
    duplicates: Dict[float, List[Path]] = {}

    for sigma in sorted(by_sigma):
        paths = sorted(by_sigma[sigma])
        unique.append((sigma, paths[0]))

        if len(paths) > 1:
            duplicates[sigma] = paths

    return unique, duplicates


def sigma_in_requested(
    sigma: float,
    requested: Sequence[float],
    tol: float = 1e-12,
) -> bool:
    return any(abs(sigma - s) <= tol for s in requested)


def select_pairs(
    unique_pairs: Sequence[Tuple[float, Path]],
    requested_sigmas: Sequence[float] | None,
) -> List[Tuple[float, Path]]:
    if requested_sigmas is None or len(requested_sigmas) == 0:
        return list(unique_pairs)

    selected = [
        pair
        for pair in unique_pairs
        if sigma_in_requested(pair[0], requested_sigmas)
    ]

    found_sigmas = [sigma for sigma, _ in unique_pairs]
    missing = [
        sigma
        for sigma in requested_sigmas
        if not sigma_in_requested(sigma, found_sigmas)
    ]

    if missing:
        print(f"[WARN] Requested sigma not found: {missing}", file=sys.stderr)

    return selected


def setup_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "font.size": 18,
            "axes.labelsize": 22,
            "axes.titlesize": 18,
            "legend.fontsize": 16,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_for_sigma(sigma: float) -> dict:
    # 見本図っぽく、白黒でも区別しやすい設定
    style_by_sigma = {
        0.005: {"linestyle": "-"},
        0.1:  {"linestyle": "--"},
        0.4:   {"linestyle": "-."},
        0.6:   {"linestyle": ":"},
        1.0:   {"linestyle": (0, (5, 1))},
    }

    for key, style in style_by_sigma.items():
        if abs(sigma - key) <= 1e-12:
            return style

    return {
        "linestyle": "-",
        # "marker": "o",
        "mfc": "none",
        "mec": "black",
    }


def main() -> int:
    args = parse_args()

    if not args.input_dir.exists() or not args.input_dir.is_dir():
        print(f"[ERROR] input directory not found: {args.input_dir}", file=sys.stderr)
        return 1

    pairs = collect_data_files(args.input_dir)

    if not pairs:
        print("[ERROR] No files matched: data,sigma*.csv", file=sys.stderr)
        return 1

    unique_pairs, duplicates = deduplicate_by_sigma(pairs)
    selected = select_pairs(unique_pairs, args.sigmas)

    if not selected:
        print("[ERROR] No sigma selected to plot.", file=sys.stderr)
        return 1

    if duplicates:
        print(
            "[WARN] Duplicate files found for same sigma; "
            "using first path in sorted order:",
            file=sys.stderr,
        )
        for sigma, paths in sorted(duplicates.items()):
            print(f"  sigma={sigma:g}", file=sys.stderr)
            for path in paths:
                print(f"    {path}", file=sys.stderr)

    setup_plot_style()

    fig, ax = plt.subplots(figsize=(10, 7))

    for sigma, path in sorted(selected, key=lambda x: x[0]):
        try:
            df = pd.read_csv(path, usecols=["generation", "worst"])
        except ValueError as e:
            print(f"[ERROR] Failed to read required columns from: {path}", file=sys.stderr)
            print("        Required columns: generation, worst", file=sys.stderr)
            print(f"        Original error: {e}", file=sys.stderr)
            return 1

        style = style_for_sigma(sigma)

        ax.plot(
            df["generation"],
            df["worst"],
            color="black",
            linestyle=style["linestyle"],
            # marker=style["marker"],
            linewidth=0.8,
            markersize=7,
            # markerfacecolor=style["mfc"],
            # markeredgecolor=style["mec"],
            markeredgewidth=0.8,
            markevery=args.markevery,
            label=fr"$\sigma={sigma:g}$",
        )

    if args.title:
        ax.set_title(args.title)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Worst fitness")
    ax.grid(False)

    if args.legend_outside:
        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            handlelength=3.5,
            borderaxespad=0.0,
            labelspacing=0.3,
        )
    else:
        ax.legend(
            loc="center right",
            frameon=False,
            handlelength=3.5,
            borderpad=0.2,
            labelspacing=0.3,
        )

    if args.xmax is not None:
        ax.set_xlim(right=args.xmax)

    if args.ymin is not None or args.ymax is not None:
        ax.set_ylim(bottom=args.ymin, top=args.ymax)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    used_sigmas = ", ".join(
        f"{sigma:g}" for sigma, _ in sorted(selected, key=lambda x: x[0])
    )

    print(f"[OK] saved: {args.output}")
    print(f"[OK] sigmas plotted: {used_sigmas}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())