#!/usr/bin/env python3

"""Plot best fitness vs generation from data,sigma*.csv files.

Example:
    python3 plot_best_by_sigma.py \
        --input-dir ../../evo_sim/results/data_sigma_only/ver1 \
        --sigmas 0.005 0.1 0.4 0.6 1.0 \
        --output best_overlay_bw.pdf
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
        description="Overlay plot: x=generation, y=best, black-and-white style by sigma"
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory that contains data,sigma*.csv files. Files are searched recursively.",
    )
    parser.add_argument(
        "--sigmas",
        type=float,
        nargs="*",
        default=None,
        help=(
            "Sigmas to plot, e.g. --sigmas 0.005 0.04 0.1. "
            "If omitted, all sigmas are used."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("best_vs_generation_overlay_bw.pdf"),
        help="Output image path. PDF is recommended for papers.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="",
        help="Plot title. Empty by default because captions usually explain figures.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Optional max x limit.",
    )
    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="Optional min y limit.",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional max y limit.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image DPI. Mainly used for PNG output.",
    )
    parser.add_argument(
        "--markevery",
        type=int,
        default=20,
        help="Interval for drawing markers. For example, 20 means one marker every 20 points.",
    )
    parser.add_argument(
        "--legend-outside",
        action="store_true",
        help="Place legend outside the plot area.",
    )
    parser.add_argument(
        "--color",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use color in the plot. Use --color to enable colors and --no-color for black lines.",
    )

    return parser.parse_args()


def extract_sigma_from_name(filename: str) -> float | None:
    match = SIGMA_FILE_PATTERN.match(filename)
    if not match:
        return None
    return float(match.group(1))


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
    """Paper-like black-and-white style."""

    plt.rcParams.update(
        {
            # Font
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "font.size": 16,
            "axes.labelsize": 20,
            "axes.titlesize": 18,
            "legend.fontsize": 15,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,

            # Axis
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,

            # Editable fonts in vector files
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def get_style_by_sigma(sigma: float) -> dict:
    """Return black-and-white line/marker style for each sigma."""

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


def get_color_by_sigma(sigma: float, palette, sigmas_sorted: Sequence[float]):
    if sigma not in sigmas_sorted:
        return "black"

    index = sigmas_sorted.index(sigma)
    return palette(index % palette.N)


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
            "[WARN] Duplicate files found for the same sigma; "
            "using the first path in sorted order:",
            file=sys.stderr,
        )
        for sigma, paths in sorted(duplicates.items()):
            print(f"  sigma={sigma:g}", file=sys.stderr)
            for path in paths:
                print(f"    {path}", file=sys.stderr)

    setup_plot_style()

    fig, ax = plt.subplots(figsize=(10, 7))
    sigmas_sorted = [sigma for sigma, _ in sorted(selected, key=lambda x: x[0])]
    palette = plt.get_cmap("tab10")

    for sigma, path in sorted(selected, key=lambda x: x[0]):
        try:
            df = pd.read_csv(path, usecols=["generation", "best"])
        except ValueError as e:
            print(f"[ERROR] Failed to read required columns from: {path}", file=sys.stderr)
            print("        Required columns: generation, best", file=sys.stderr)
            print(f"        Original error: {e}", file=sys.stderr)
            return 1

        style = get_style_by_sigma(sigma)
        line_color = get_color_by_sigma(sigma, palette, sigmas_sorted) if args.color else "black"

        ax.plot(
            df["generation"],
            df["best"],
            color=line_color,
            linestyle=style["linestyle"],
            # marker=style["marker"],
            linewidth=1.1,
            markersize=6.5,
            # markerfacecolor=style["mfc"],
            # markeredgecolor=style["mec"],
            markeredgewidth=0.9,
            markevery=args.markevery,
            label=fr"$\sigma={sigma:g}$",
        )

    if args.title:
        ax.set_title(args.title)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Best fitness")

    # Paper-like style: no grid
    ax.grid(False)

    if args.legend_outside:
        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            borderaxespad=0.0,
            handlelength=3.5,
            labelspacing=0.35,
        )
    else:
        ax.legend(
            loc="center right",
            frameon=False,
            handlelength=3.5,
            labelspacing=0.35,
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