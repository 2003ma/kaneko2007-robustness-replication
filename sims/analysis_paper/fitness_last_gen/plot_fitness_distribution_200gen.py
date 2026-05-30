#!/usr/bin/env python3

"""Overlay fitness distributions for two 200-generation CSV files.

The script plots histogram counts of the `fitness` column for two CSV files.
It supports two fixed x-range presets:
  - wide: fitness in [overall min, overall max]
  - zoom: fitness in [-0.2, 0]

Example:
    python3 ./plot_fitness_distribution_200gen.py \
        --csv-a ../../evo_sim/results/ver1/sigma_0.005/evo_sim_data/gen_200_all_J_sigma_0.005_dt0.005.csv \
        --csv-b ../../evo_sim/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
        --labels sigma=0.005 sigma=0.200 \
        --preset wide \
        --output fitness_wide.pdf


    python3 ./plot_fitness_distribution_200gen.py \
        --csv-a ../../evo_sim/results/ver1/sigma_0.005/evo_sim_data/gen_200_all_J_sigma_0.005_dt0.005.csv \
        --csv-b ../../evo_sim/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
        --labels sigma=0.005 sigma=0.200 \
        --preset zoom \
        --output fitness_zoom.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence, Tuple


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Overlay fitness distributions from two CSV files")
    parser.add_argument("--csv-a", type=Path, required=True, help="First CSV file")
    parser.add_argument("--csv-b", type=Path, required=True, help="Second CSV file")
    parser.add_argument(
        "--labels",
        nargs=2,
        default=["csv-a", "csv-b"],
        help="Legend labels for the two CSV files",
    )
    parser.add_argument(
        "--preset",
        choices=["wide", "zoom"],
        default="wide",
        help="Fixed x-range preset",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fitness_distribution_overlay.pdf"),
        help="Output image path",
    )
    parser.add_argument("--bins", type=int, default=80, help="Number of histogram bins")
    parser.add_argument("--dpi", type=int, default=300, help="Output image DPI")
    return parser.parse_args(argv)


def import_plotting_modules():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        print(
            "[ERROR] matplotlib is not installed. Install it with `pip install matplotlib pandas`.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        print(
            "[ERROR] pandas is not installed. Install it with `pip install matplotlib pandas`.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    return plt, pd


def resolve_xlim(preset: str) -> Tuple[float, float]: #zoom時のみ
    return -0.2, 0.0


def load_fitness_values(pd, csv_path: Path):
    if not csv_path.exists() or not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    df = pd.read_csv(csv_path, usecols=["fitness"])
    values = df["fitness"].dropna()
    return values.to_numpy()


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    plt, pd = import_plotting_modules()

    try:
        fitness_a = load_fitness_values(pd, args.csv_a)
        fitness_b = load_fitness_values(pd, args.csv_b)
    except FileNotFoundError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    bins = args.bins
    
    # wide プリセット: データの全体範囲を使う
    # zoom プリセット: 固定範囲を使う
    if args.preset == "wide":
        x_min = min(fitness_a.min(), fitness_b.min())
        x_max = max(fitness_a.max(), fitness_b.max())
    else:  # zoom
        x_min, x_max = resolve_xlim(args.preset)
    
    # 各境界を計算
    bin_edges = [x_min + (x_max - x_min) * i / bins for i in range(bins + 1)]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(
        fitness_a,
        bins=bin_edges,
        histtype="step",
        linewidth=2.0,
        color="black",
        label=args.labels[0],
    )
    ax.hist(
        fitness_b,
        bins=bin_edges,
        histtype="step",
        linewidth=2.0,
        color="black",
        linestyle="--",
        label=args.labels[1],
    )

    ax.set_xlim(x_min, x_max)
    ax.set_yscale("log")
    ax.set_ylim(bottom=1)
    ax.set_xlabel("fitness")
    ax.set_ylabel("individual count")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)

    print(f"[OK] saved: {args.output}")
    print(f"[OK] preset: {args.preset}")
    print(f"[OK] csv-a: {args.csv_a}")
    print(f"[OK] csv-b: {args.csv_b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
