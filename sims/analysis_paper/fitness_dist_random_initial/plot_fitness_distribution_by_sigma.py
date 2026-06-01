#!/usr/bin/env python3

"""Overlay fitness distributions by sigma from CSV files.

Expected CSV format:
    trial_id,fitness,...

Typical filename example:
    trial_fitness_gen=200_ind=0_sigma=0.005_seed=12345_trials=10000_dt0.050_allx.csv

python3 plot_fitness_distribution_by_sigma.py

This script scans CSV files, groups them by sigma value parsed from filename,
and overlays histogram curves on one figure.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SIGMA_PATTERN = re.compile(r"(?:^|_)sigma=([0-9]*\.?[0-9]+)")
BIN_SIZE = 0.1


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Overlay fitness distributions by sigma")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "fitness",
        help="Directory containing fitness CSV files",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.csv",
        help="Glob pattern for input CSV files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "fitness_distribution_overlay_by_sigma.pdf",
        help="Output image path",
    )
    parser.add_argument("--x-min", type=float, default=-8.1, help="x-axis minimum")
    parser.add_argument("--x-max", type=float, default=0.1, help="x-axis maximum")
    parser.add_argument(
        "--density",
        action="store_true",
        help="Plot normalized density instead of raw counts",
    )
    parser.add_argument("--dpi", type=int, default=160, help="Output DPI")
    return parser.parse_args(argv)


def extract_sigma(path: Path) -> float | None:
    m = SIGMA_PATTERN.search(path.name)
    if m is None:
        return None
    return float(m.group(1))


def load_fitness(csv_path: Path) -> pd.Series:
    df = pd.read_csv(csv_path, usecols=["fitness"])
    return df["fitness"].dropna()


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    if not args.input_dir.exists() or not args.input_dir.is_dir():
        print(f"[ERROR] input directory not found: {args.input_dir}", file=sys.stderr)
        return 1
    if args.x_min >= args.x_max:
        print("[ERROR] x-min must be smaller than x-max", file=sys.stderr)
        return 1

    files = sorted(args.input_dir.glob(args.pattern))
    if not files:
        print(f"[ERROR] no CSV files matched: {args.input_dir / args.pattern}", file=sys.stderr)
        return 1

    grouped: Dict[float, List[pd.Series]] = {}
    skipped = 0

    for path in files:
        sigma = extract_sigma(path)
        if sigma is None:
            skipped += 1
            continue
        try:
            grouped.setdefault(sigma, []).append(load_fitness(path))
        except ValueError:
            print(f"[WARN] missing 'fitness' column, skipped: {path}", file=sys.stderr)
            skipped += 1

    if not grouped:
        print("[ERROR] no valid sigma-grouped CSV files found", file=sys.stderr)
        return 1

    fig, ax = plt.subplots(figsize=(9, 6))
    # Fixed-width bins so distributions are directly comparable across sigmas.
    bin_edges = np.arange(args.x_min, args.x_max + BIN_SIZE * 0.5, BIN_SIZE)
    style_map = {
        0.005: {"color": "tab:blue",   "linestyle": "-",  "linewidth": 2.4},
        0.04:  {"color": "tab:orange", "linestyle": "--", "linewidth": 2.2},
        0.1:   {"color": "tab:green",  "linestyle": ":",  "linewidth": 2.8},
    }
    style = style_map[sigma]

    for sigma in sorted(grouped.keys()):
        style = style_map.get(
            sigma,
            {"linestyle": "-", "linewidth": 2.0},
        )
        values = pd.concat(grouped[sigma], ignore_index=True)
        ax.hist(
            values,
            bins=bin_edges,
            histtype="step",
            linewidth=2.0,
            linestyle=style["linestyle"],
            color=style["color"],
            density=args.density,
            label=f"sigma={sigma:.3f}",
        )

    ax.set_xlim(args.x_min, args.x_max)
    ax.set_xlabel("fitness")
    ax.set_ylabel("density" if args.density else "count")
    # ax.set_title("Fitness distribution overlay by sigma")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)

    print(f"[OK] saved: {args.output}")
    print(f"[OK] sigma count: {len(grouped)}")
    if skipped > 0:
        print(f"[WARN] skipped files: {skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
