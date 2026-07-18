#!/usr/bin/env python3

"""Plot zero-fitness ratio vs simulation noise strength, overlaid by evo_sigma.

Expected input CSV format (wide):
    trial_id,0_fitness,0.01_fitness,0.02_fitness,...

Each CSV filename should include evo_sigma, e.g.:
    trial_fitness_gen=200_ind=0_evo_sigma=0.040_sim_sigma=0.000_0.300_step_0.010_seed=12345_trials=10000.csv

Example:
    ./plot_zero_fitness_ratio_vs_sim_sigma.py \
      --input-dir . \
      --output under_noise_fitness.pdf
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd

FILENAME_PATTERN = re.compile(r"evo_sigma=([0-9]*\.?[0-9]+)")
COLUMN_PATTERN = re.compile(r"^([0-9]*\.?[0-9]+)_fitness$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay plot: x=simulation sigma, y=ratio of zero fitness"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("."),
        help="Directory containing trial_fitness CSV files (searched recursively)",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="trial_fitness_gen=*_ind=*_evo_sigma=*_sim_sigma=*_seed=*_trials=*.csv",
        help="Glob pattern for target CSV files",
    )
    parser.add_argument(
        "--evo-sigmas",
        type=float,
        nargs="*",
        default=None,
        help="Optional evo_sigma values to plot. If omitted, all found values are used.",
    )
    parser.add_argument(
        "--zero-tol",
        type=float,
        default=1e-12,
        help="Absolute tolerance for treating fitness as zero",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Zero-fitness ratio vs simulation noise",
        help="Plot title",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("zero_fitness_ratio_vs_sim_sigma.pdf"),
        help="Output image path",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output image DPI",
    )
    return parser.parse_args()


def extract_evo_sigma(path: Path) -> float | None:
    match = FILENAME_PATTERN.search(path.name)
    if match is None:
        return None
    return float(match.group(1))


def sigma_requested(sigma: float, requested: Sequence[float], tol: float = 1e-12) -> bool:
    return any(abs(sigma - req) <= tol for req in requested)


def collect_zero_counts(csv_path: Path, zero_tol: float) -> Tuple[Dict[float, int], Dict[float, int]]:
    try:
        df = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        return {}, {}

    zero_counts: Dict[float, int] = {}
    total_counts: Dict[float, int] = {}

    for col in df.columns:
        match = COLUMN_PATTERN.match(col)
        if match is None:
            continue
        sim_sigma = float(match.group(1))
        values = df[col]
        is_zero = values.abs() <= zero_tol
        zero_counts[sim_sigma] = int(is_zero.sum())
        total_counts[sim_sigma] = int(values.shape[0])

    if not zero_counts:
        return {}, {}

    return zero_counts, total_counts


def main() -> int:
    args = parse_args()

    if not args.input_dir.exists() or not args.input_dir.is_dir():
        print(f"[ERROR] input directory not found: {args.input_dir}", file=sys.stderr)
        return 1

    files = sorted(args.input_dir.rglob(args.pattern))
    if not files:
        print(f"[ERROR] No files matched pattern: {args.pattern}", file=sys.stderr)
        return 1

    grouped: Dict[float, List[Path]] = {}
    skipped = 0
    for path in files:
        evo_sigma = extract_evo_sigma(path)
        if evo_sigma is None:
            skipped += 1
            continue
        if args.evo_sigmas is not None and len(args.evo_sigmas) > 0:
            if not sigma_requested(evo_sigma, args.evo_sigmas):
                continue
        grouped.setdefault(evo_sigma, []).append(path)

    if not grouped:
        print("[ERROR] No valid files found after filtering evo_sigma.", file=sys.stderr)
        return 1

    plt.figure(figsize=(8, 6))
    line_styles = ["-", "--", "-.", ":", (0, (5, 2)), (0, (3, 1, 1, 1))]
    markers = ["o", "s", "^", "D", "v", "x"]

    for i, evo_sigma in enumerate(sorted(grouped.keys())):
        paths = grouped[evo_sigma]
        zero_sum: Dict[float, int] = {}
        total_sum: Dict[float, int] = {}
        empty_or_invalid = 0

        for path in paths:
            zero_counts, total_counts = collect_zero_counts(path, args.zero_tol)
            if not zero_counts:
                empty_or_invalid += 1
                continue
            for sim_sigma, zc in zero_counts.items():
                zero_sum[sim_sigma] = zero_sum.get(sim_sigma, 0) + zc
                total_sum[sim_sigma] = total_sum.get(sim_sigma, 0) + total_counts[sim_sigma]

        sim_sigmas = sorted(total_sum.keys())
        ratios = [zero_sum[s] / total_sum[s] if total_sum[s] > 0 else 0.0 for s in sim_sigmas]

        plt.plot(
            sim_sigmas,
            ratios,
            linestyle=line_styles[i % len(line_styles)],
            marker=markers[i % len(markers)],
            linewidth=2.0,
            markersize=6,
            label=fr"$\sigma_{{evo}}$={evo_sigma:g}",
        )

        if empty_or_invalid > 0:
            print(
                f"[WARN] skipped {empty_or_invalid} empty or invalid CSV file(s) for evo_sigma={evo_sigma:g}",
                file=sys.stderr,
            )

    # plt.title(args.title)
    plt.xlabel(r"Evaluation noise strength $\sigma_{\mathrm{sim}}$",fontsize=16)
    plt.ylabel(r"Fraction of individuals with $F=0$",fontsize=16)
    plt.ylim(0.0, 1.0)
    plt.grid(True, alpha=0.25)
    plt.legend(
        title=r"Evolution noise $\sigma_{\mathrm{evo}}$",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        fontsize=15,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi)

    print(f"[OK] saved: {args.output}")
    print(f"[OK] evo_sigma count: {len(grouped)}")
    if skipped > 0:
        print(f"[WARN] skipped files without evo_sigma in name: {skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
