#!/usr/bin/env python3

"""Plot best fitness vs generation from data,sigma*.csv files.Fig.1

Example:
    ./plot_best_by_sigma.py \
        --input-dir ../../evo_sim/results/data_sigma_only/ver1 \
        --sigmas 0.005 0.01 0.04 0.1 0.2 0.4 0.6 0.7 0.8 0.9 1.0\
        --output best_overlay_ver1.png
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
        description="Overlay plot: x=generation, y=best, color by sigma"
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
        help="Sigmas to plot (e.g. --sigmas 0.005 0.01 0.04). If omitted, all are used.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("best_vs_generation_overlay.png"),
        help="Output image path",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Top fitness vs generation",
        help="Plot title",
    )
    parser.add_argument("--xmax", type=float, default=None, help="Optional max x limit")
    parser.add_argument("--ymin", type=float, default=None, help="Optional min y limit")
    parser.add_argument("--ymax", type=float, default=None, help="Optional max y limit")
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output image DPI",
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


def deduplicate_by_sigma(pairs: Sequence[Tuple[float, Path]]) -> Tuple[List[Tuple[float, Path]], Dict[float, List[Path]]]:
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
    return any(abs(sigma - s) <= tol for s in requested)


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

    if args.sigmas is not None and len(args.sigmas) > 0:
        selected = [p for p in unique_pairs if sigma_in_requested(p[0], args.sigmas)]
        missing = sorted({float(s) for s in args.sigmas if not any(sigma_in_requested(sigma, [float(s)]) for sigma, _ in unique_pairs)})
        if missing:
            print(f"[WARN] Requested sigma not found: {missing}", file=sys.stderr)
    else:
        selected = unique_pairs

    if not selected:
        print("[ERROR] No sigma selected to plot.", file=sys.stderr)
        return 1

    if duplicates:
        print("[WARN] Duplicate files found for same sigma; using first path in sorted order:", file=sys.stderr)
        for sigma, paths in sorted(duplicates.items()):
            print(f"  sigma={sigma:.4f}", file=sys.stderr)
            for p in paths:
                print(f"    {p}", file=sys.stderr)

    plt.figure(figsize=(8, 6))

    for sigma, path in sorted(selected, key=lambda x: x[0]):
        df = pd.read_csv(path, usecols=["generation", "best"])
        plt.plot(df["generation"], df["best"], linewidth=2.0, label=fr"$\sigma$={sigma:g}")

    plt.title(args.title)
    plt.xlabel("generation")
    plt.ylabel("Top Fitness")
    plt.grid(True, alpha=0.25)
    plt.legend(frameon=False)

    if args.xmax is not None:
        plt.xlim(right=args.xmax)
    if args.ymin is not None or args.ymax is not None:
        plt.ylim(bottom=args.ymin, top=args.ymax)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi)

    used_sigmas = ", ".join(f"{s:g}" for s, _ in sorted(selected, key=lambda x: x[0]))
    print(f"[OK] saved: {args.output}")
    print(f"[OK] sigmas plotted: {used_sigmas}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
