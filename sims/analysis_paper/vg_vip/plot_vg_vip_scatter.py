#!/usr/bin/env python3
"""
python3 ./plot_vg_vip_scatter.py \
  --input-dir ../../evo_sim/results/data_sigma_only/ver1 \
  --sigmas 0.01 0.04 0.1 0.2 0.3 0.5
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

SIGMA_FILE_PATTERN = re.compile(r"^data,sigma=([0-9]*\.?[0-9]+),.*\.csv$")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scatter plot of Vg vs Vip by sigma")
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory that contains data,sigma*.csv files (searched recursively)",
    )
    parser.add_argument(
        "--sigmas",
        type=float,
        nargs="*",
        default=None,
        help="Sigmas to plot. If omitted, all are used.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("vg_vip.png"),
        help="Output image path",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="V_g vs Vip",
        help="Plot title",
    )
    parser.add_argument("--xlim-min", type=float, default=None)
    parser.add_argument("--xlim-max", type=float, default=None)
    parser.add_argument("--ylim-min", type=float, default=None)
    parser.add_argument("--ylim-max", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=150, help="Output image DPI")
    return parser.parse_args(argv)


def import_plotting_modules():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        print("[ERROR] matplotlib is not installed. Install it with `pip install matplotlib pandas`.", file=sys.stderr)
        raise SystemExit(1) from exc

    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        print("[ERROR] pandas is not installed. Install it with `pip install matplotlib pandas`.", file=sys.stderr)
        raise SystemExit(1) from exc

    return plt, pd


def extract_sigma_from_name(filename: str) -> float | None:
    match = SIGMA_FILE_PATTERN.match(filename)
    if match is None:
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


def sigma_in_requested(sigma: float, requested: Sequence[float], tol: float = 1e-12) -> bool:
    return any(abs(sigma - requested_sigma) <= tol for requested_sigma in requested)


def resolve_axis_limits(
    preset: str,
    override_min: float | None,
    override_max: float | None,
    kind: str,
) -> Tuple[float, float]:
    if kind == "x":
        base_min, base_max = (1e-3, 1e1)
    else:
        base_min, base_max = (1e-5, 1e0)

    if override_min is not None:
        base_min = override_min
    if override_max is not None:
        base_max = override_max
    return base_min, base_max

def marker_for_sigma(sigma: float) -> str:
    if sigma in {0.01, 0.04}:
        return "o"   # small noise
    elif sigma in {0.1, 0.2}:
        return "^"   # medium noise
    elif sigma in {0.3, 0.5}:
        return "s"   # large noise
    else:
        return "o"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    plt, pd = import_plotting_modules()

    if not args.input_dir.exists() or not args.input_dir.is_dir():
        print(f"[ERROR] input directory not found: {args.input_dir}", file=sys.stderr)
        return 1

    pairs = collect_data_files(args.input_dir)
    if not pairs:
        print("[ERROR] No files matched: data,sigma*.csv", file=sys.stderr)
        return 1

    unique_pairs, duplicates = deduplicate_by_sigma(pairs)

    if args.sigmas is not None and len(args.sigmas) > 0:
        selected = [pair for pair in unique_pairs if sigma_in_requested(pair[0], args.sigmas)]
        missing = sorted(
            {
                float(requested_sigma)
                for requested_sigma in args.sigmas
                if not any(sigma_in_requested(sigma, [float(requested_sigma)]) for sigma, _ in unique_pairs)
            }
        )
        if missing:
            print(f"[WARN] Requested sigma not found: {missing}", file=sys.stderr)
    else:
        selected = unique_pairs

    if not selected:
        print("[ERROR] No sigma selected to plot.", file=sys.stderr)
        return 1

    if duplicates:
        print("[WARN] Duplicate files found for the same sigma; using the first path in sorted order:", file=sys.stderr)
        for sigma, paths in sorted(duplicates.items()):
            print(f"  sigma={sigma:.4f}", file=sys.stderr)
            for path in paths:
                print(f"    {path}", file=sys.stderr)

    selected = sorted(selected, key=lambda item: item[0])

    cmap = plt.get_cmap("turbo", len(selected) if len(selected) > 1 else 2)
    fig, ax = plt.subplots(figsize=(8.5, 7.0))
    ax.set_xscale("log")
    ax.set_yscale("log")
    cmap = plt.get_cmap("tab10")
    small_noise = {0.01, 0.04}
    medium_noise = {0.1, 0.2}
    large_noise = {0.3, 0.5}

    for index, (sigma, path) in enumerate(selected):
        color = cmap(index % 10)
        df = pd.read_csv(path, usecols=["Vip", "V_g"])
        df = df.iloc[::4]
        df = df[(df["Vip"] > 0) & (df["V_g"] > 0)]
        if df.empty:
            print(f"[WARN] No positive Vip/V_g rows in {path}; skipping.", file=sys.stderr)
            continue
        ax.scatter(
            df["Vip"],
            df["V_g"],
            s=25,
            alpha=0.8                                       ,
            color=color,
            # color="black",
            edgecolors="none",
            label=fr"$\sigma$={sigma:g}",
            marker=marker_for_sigma(sigma),
        )

    x_min, x_max = resolve_axis_limits("wide", args.xlim_min, args.xlim_max, "x")
    y_min, y_max = resolve_axis_limits("wide", args.ylim_min, args.ylim_max, "y")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    line_min = max(x_min, y_min)
    line_max = min(x_max, y_max)
    ax.plot([line_min, line_max], [line_min, line_max], color="black", linestyle="--", linewidth=1.0, alpha=1.0,label=r"$V_g = V_{ip}$",)


    ax.set_xlabel(r"$V_{ip}$",fontsize=20)
    ax.set_ylabel(r"$V_g$",fontsize=20)
    ax.grid(True, which="both", alpha=0.2)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.tick_params(axis="both", which="minor", labelsize=15)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0),fontsize=15,title=r"Evolution noise $\sigma_{\mathrm{evo}}$")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")

    used_sigmas = ", ".join(f"{sigma:g}" for sigma, _ in selected)
    print(f"[OK] saved: {args.output}")
    print(f"[OK] sigmas plotted: {used_sigmas}")
    print("[OK] preset: wide")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
