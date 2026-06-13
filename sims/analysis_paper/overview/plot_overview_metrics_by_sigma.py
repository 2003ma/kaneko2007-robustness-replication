#!/usr/bin/env python3

"""Create an overview plot from data,sigma*.csv files. Fig.2

This script creates two separated panels:
- top panel: mean and worst fitness averaged over generations 100-200
- bottom panel: speed defined as 1 / first generation where best fitness becomes 0

If best fitness does not become 0 by generation 200, the speed is 0.

Example:
    python3 ./plot_overview_metrics_by_sigma.py \
        --input-dir ../../evo_sim/results/data_sigma_only/ver1 \
        --sigmas 0.005 0.01 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.9 1.0 \
        --output overview.png
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

SIGMA_FILE_PATTERN = re.compile(r"^data,sigma=([0-9]*\.?[0-9]+),.*\.csv$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot overview metrics against sigma from data,sigma*.csv files"
    )
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
        default=Path("overview_metrics.png"),
        help="Output image path",
    )
    parser.add_argument(
        "--window-start",
        type=int,
        default=100,
        help="Start generation for mean/worst averaging window",
    )
    parser.add_argument(
        "--window-end",
        type=int,
        default=200,
        help="End generation for mean/worst averaging window",
    )
    parser.add_argument(
        "--zero-limit",
        type=int,
        default=200,
        help="Search for first best==0 generation up to this generation",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Overview metrics by noise strength",
        help="Plot title",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output image DPI",
    )
    return parser.parse_args()


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


def compute_metrics(pd, path: Path, window_start: int, window_end: int, zero_limit: int) -> Tuple[float, float, float]:
    df = pd.read_csv(path, usecols=["generation", "best", "mean", "worst"])

    window = df[(df["generation"] >= window_start) & (df["generation"] <= window_end)]
    if window.empty:
        raise ValueError(f"No rows found in generation window [{window_start}, {window_end}] for {path}")

    mean_average = float(window["mean"].mean())
    worst_average = float(window["worst"].mean())

    zero_window = df[(df["generation"] <= zero_limit) & (df["best"] == 0)]
    if zero_window.empty:
        speed = 0.0
    else:
        first_zero_generation = float(zero_window.iloc[0]["generation"])
        speed = 0.0 if first_zero_generation <= 0 else 1.0 / first_zero_generation

    return mean_average, worst_average, speed


def main() -> int:
    args = parse_args()
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

    sigma_values: List[float] = []
    mean_values: List[float] = []
    worst_values: List[float] = []
    speed_values: List[float] = []

    for sigma, path in selected:
        mean_average, worst_average, speed = compute_metrics(
            pd,
            path,
            args.window_start,
            args.window_end,
            args.zero_limit,
        )
        sigma_values.append(sigma)
        mean_values.append(mean_average)
        worst_values.append(worst_average)
        speed_values.append(speed)

    fig, (ax_mean_worst, ax_speed) = plt.subplots(
        2,
        1,
        figsize=(9, 8.6),
        sharex=True,
        gridspec_kw={"height_ratios": [1.45, 1.0]},
    )

    # sigma は値の間隔が不均一なので、見やすさのため等間隔の x 座標に配置する
    x_positions = list(range(len(sigma_values)))

    ax_mean_worst.plot(
        x_positions,
        mean_values,
        color="#d62728",
        marker="o",
        linewidth=2.0,
        label=r"Generation average of population mean $\bar{F}$",
    )
    ax_mean_worst.plot(
        x_positions,
        worst_values,
        color="#2ca02c",
        marker="s",
        linestyle="--",
        linewidth=2.0,
        label=r"Generation average of population worst $\bar{F}$",
    )

    ax_speed.plot(
        x_positions,
        speed_values,
        color="#1f77b4",
        marker="^",
        linestyle=":",
        linewidth=2.0,
        label=r"1 / $g_{first}$",
    )

    ax_speed.set_xticks(x_positions)
    ax_speed.set_xticklabels([f"{sigma:g}" for sigma in sigma_values])
    if x_positions:
        ax_speed.set_xlim(-0.5, x_positions[-1] + 0.5)

    ax_mean_worst.axhline(0.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax_speed.axhline(0.0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)

    # ax_mean_worst.set_title(args.title)
    ax_mean_worst.set_ylabel(r"$\bar{F}$ statistics,generations 100-200",fontsize=15)
    ax_speed.set_xlabel("noise strength σ",fontsize=15)
    ax_speed.set_ylabel(r"1 / $g_{first}$",fontsize=15)
    ax_speed.set_ylim(-0.005, 0.13)

    ax_mean_worst.grid(True, which="both", alpha=0.25)
    ax_speed.grid(True, which="both", alpha=0.25)
    ax_mean_worst.legend(frameon=False)
    ax_speed.legend(frameon=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)

    print(f"[OK] saved: {args.output}")
    print("[OK] metrics plotted: mean+worst (top), speed only (bottom)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
