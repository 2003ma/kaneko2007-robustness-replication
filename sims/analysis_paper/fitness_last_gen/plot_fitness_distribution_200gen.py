# #!/usr/bin/env python3

# """Overlay fitness distributions for two 200-generation CSV files.

# The script plots histogram counts of the `fitness` column for two CSV files.
# It supports two fixed x-range presets:
#   - wide: fitness in [overall min, overall max]
#   - zoom: fitness in [-0.2, 0]

# Example:
#     python3 ./plot_fitness_distribution_200gen.py \
#         --csv-a ../../evo_sim/results/ver1/sigma_0.005/evo_sim_data/gen_200_all_J_sigma_0.005_dt0.005.csv \
#         --csv-b ../../evo_sim/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
#         --labels sigma=0.005 sigma=0.200 \
#         --preset wide \
#         --output fitness_wide.pdf


#     python3 ./plot_fitness_distribution_200gen.py \
#         --csv-a ../../evo_sim/results/ver1/sigma_0.005/evo_sim_data/gen_200_all_J_sigma_0.005_dt0.005.csv \
#         --csv-b ../../evo_sim/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
#         --labels sigma=0.005 sigma=0.200 \
#         --preset zoom \
#         --output fitness_zoom.png
# """

# from __future__ import annotations

# import argparse
# import sys
# from pathlib import Path
# from typing import Sequence, Tuple


# def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
#     parser = argparse.ArgumentParser(description="Overlay fitness distributions from two CSV files")
#     parser.add_argument("--csv-a", type=Path, required=True, help="First CSV file")
#     parser.add_argument("--csv-b", type=Path, required=True, help="Second CSV file")
#     parser.add_argument(
#         "--labels",
#         nargs=2,
#         default=["csv-a", "csv-b"],
#         help="Legend labels for the two CSV files",
#     )
#     parser.add_argument(
#         "--preset",
#         choices=["wide", "zoom"],
#         default="wide",
#         help="Fixed x-range preset",
#     )
#     parser.add_argument(
#         "--output",
#         type=Path,
#         default=Path("fitness_distribution_overlay.pdf"),
#         help="Output image path",
#     )
#     parser.add_argument("--bins", type=int, default=80, help="Number of histogram bins")
#     parser.add_argument("--dpi", type=int, default=300, help="Output image DPI")
#     return parser.parse_args(argv)


# def import_plotting_modules():
#     try:
#         import matplotlib.pyplot as plt
#     except ModuleNotFoundError as exc:
#         print(
#             "[ERROR] matplotlib is not installed. Install it with `pip install matplotlib pandas`.",
#             file=sys.stderr,
#         )
#         raise SystemExit(1) from exc

#     try:
#         import pandas as pd
#     except ModuleNotFoundError as exc:
#         print(
#             "[ERROR] pandas is not installed. Install it with `pip install matplotlib pandas`.",
#             file=sys.stderr,
#         )
#         raise SystemExit(1) from exc

#     return plt, pd


# def resolve_xlim(preset: str) -> Tuple[float, float]: #zoom時のみ
#     return -0.2, 0.0


# def load_fitness_values(pd, csv_path: Path):
#     if not csv_path.exists() or not csv_path.is_file():
#         raise FileNotFoundError(f"CSV not found: {csv_path}")
#     df = pd.read_csv(csv_path, usecols=["fitness"])
#     values = df["fitness"].dropna()
#     return values.to_numpy()


# def main(argv: Sequence[str] | None = None) -> int:
#     args = parse_args(argv)
#     plt, pd = import_plotting_modules()

#     try:
#         fitness_a = load_fitness_values(pd, args.csv_a)
#         fitness_b = load_fitness_values(pd, args.csv_b)
#     except FileNotFoundError as exc:
#         print(f"[ERROR] {exc}", file=sys.stderr)
#         return 1

#     bins = args.bins
    
#     # wide プリセット: データの全体範囲を使う
#     # zoom プリセット: 固定範囲を使う
#     if args.preset == "wide":
#         x_min = min(fitness_a.min(), fitness_b.min())
#         x_max = max(fitness_a.max(), fitness_b.max())
#     else:  # zoom
#         x_min, x_max = resolve_xlim(args.preset)
    
#     # 各境界を計算
#     bin_edges = [x_min + (x_max - x_min) * i / bins for i in range(bins + 1)]

#     fig, ax = plt.subplots(figsize=(8, 6))
#     ax.hist(
#         fitness_a,
#         bins=bin_edges,
#         histtype="step",
#         linewidth=2.0,
#         color="black",
#         label=args.labels[0],
#     )
#     ax.hist(
#         fitness_b,
#         bins=bin_edges,
#         histtype="step",
#         linewidth=2.0,
#         color="black",
#         linestyle="--",
#         label=args.labels[1],
#     )

#     ax.set_xlim(x_min, x_max)
#     ax.set_yscale("log")
#     ax.set_ylim(bottom=1)
#     ax.set_xlabel("fitness")
#     ax.set_ylabel("individual count")
#     ax.grid(True, which="both", alpha=0.25)
#     ax.legend(frameon=False)

#     fig.tight_layout()
#     args.output.parent.mkdir(parents=True, exist_ok=True)
#     fig.savefig(args.output, dpi=args.dpi)

#     print(f"[OK] saved: {args.output}")
#     print(f"[OK] preset: {args.preset}")
#     print(f"[OK] csv-a: {args.csv_a}")
#     print(f"[OK] csv-b: {args.csv_b}")
#     return 0


# if __name__ == "__main__":
#     raise SystemExit(main())

#!/usr/bin/env python3

"""Overlay fitness distributions for two 200-generation CSV files.

The script plots histogram counts of the `fitness` column for two CSV files.

Presets:
  - wide  : main plot over the full fitness range
  - zoom  : main plot over [-0.2, 0]
  - inset : wide plot with an inset zooming into [-0.2, 0]

Example:
    python3 ./plot_fitness_distribution_200gen.py \
        --csv-a ../../evo_sim/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
        --csv-b ../../evo_sim/results/ver1/sigma_0.005/evo_sim_data/gen_200_all_J_sigma_0.005_dt0.005.csv \
        --labels "sigma=0.200" "sigma=0.005" \
        --preset wide \
        --output fitness_wide.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence, Tuple


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay fitness distributions from two CSV files"
    )
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
        choices=["wide", "zoom", "inset"],
        default="wide",
        help="Plot type: wide, zoom, or inset",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fitness_distribution_overlay.pdf"),
        help="Output image path",
    )
    parser.add_argument("--bins", type=int, default=80, help="Number of histogram bins")
    parser.add_argument("--dpi", type=int, default=300, help="Output image DPI")
    parser.add_argument(
        "--xlim",
        nargs=2,
        type=float,
        default=None,
        metavar=("XMIN", "XMAX"),
        help="Main panel x-range. If omitted, wide/inset use data range and zoom uses [-0.2, 0].",
    )
    parser.add_argument(
        "--inset-xlim",
        nargs=2,
        type=float,
        default=[-0.2, 0.0],
        metavar=("XMIN", "XMAX"),
        help="Inset x-range used only for --preset inset.",
    )
    parser.add_argument(
        "--color",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use red/green colors. Use --color to enable and --no-color for black-only style.",
    )
    return parser.parse_args(argv)


def import_plotting_modules():
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ModuleNotFoundError as exc:
        print(
            "[ERROR] matplotlib and numpy are required. Install them with `pip install matplotlib numpy pandas`.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        print(
            "[ERROR] pandas is required. Install it with `pip install pandas`.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    return plt, np, pd


def setup_plot_style(plt) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 15,
            "axes.labelsize": 15,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 15,
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def resolve_xlim(
    preset: str,
    fitness_a,
    fitness_b,
    xlim: Sequence[float] | None,
) -> Tuple[float, float]:
    if xlim is not None:
        return float(xlim[0]), float(xlim[1])

    if preset == "zoom":
        return -0.2, 0.0

    x_min = min(float(fitness_a.min()), float(fitness_b.min()))
    x_max = max(float(fitness_a.max()), float(fitness_b.max()))

    # 0 付近が枠に重なりすぎないように少し右側へ余白
    if x_max <= 0.0:
        x_max = 0.0
    x_pad = 0.01 * (x_max - x_min)

    return x_min, x_max + x_pad


def load_fitness_values(pd, csv_path: Path):
    if not csv_path.exists() or not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path, usecols=["fitness"])
    return df["fitness"].dropna().to_numpy()


def make_histogram(np, values, bin_edges):
    counts, edges = np.histogram(values, bins=bin_edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Keep zero-count bins so they remain visible on a symlog axis.
    counts = counts.astype(float)

    return centers, counts


def get_styles(use_color: bool):
    if use_color:
        return (
            {
                "color": "#e92b2b",
                "linestyle": "-",
                "marker": "+",
            },
            {
                "color": "#00a651",
                "linestyle": "--",
                "marker": "x",
            },
        )

    return (
        {
            "color": "black",
            "linestyle": "-",
            "marker": "+",
        },
        {
            "color": "black",
            "linestyle": "--",
            "marker": "x",
        },
    )


def plot_distribution(
    ax,
    centers,
    counts,
    *,
    label: str,
    color: str,
    linestyle,
    marker: str,
    linewidth: float = 2.4,
    markersize: float = 6.0,
):
    ax.plot(
        centers,
        counts,
        color=color,
        linestyle=linestyle,
        marker=marker,
        linewidth=linewidth,
        markersize=markersize,
        markeredgewidth=1.6,
        label=label,
    )


def decorate_axis(ax, xlim: Tuple[float, float], *, show_legend: bool = True) -> None:
    ax.set_xlim(*xlim)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Fitness",fontsize=15)
    ax.set_ylabel("Distribution",fontsize=15)
    ax.grid(False)

    if show_legend:
        ax.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.52, 0.90),
            ncol=1,
            handlelength=3.0,
            borderpad=0.2,
            labelspacing=0.2,
        )


def plot_panel(
    np,
    ax,
    fitness_a,
    fitness_b,
    xlim: Tuple[float, float],
    bins: int,
    labels: Sequence[str],
    style_a: dict,
    style_b: dict,
    *,
    show_legend: bool = True,
    linewidth: float = 2.4,
    markersize: float = 6.0,
):
    bin_edges = np.linspace(xlim[0], xlim[1], bins + 1)

    centers_a, counts_a = make_histogram(np, fitness_a, bin_edges)
    centers_b, counts_b = make_histogram(np, fitness_b, bin_edges)

    plot_distribution(
        ax,
        centers_a,
        counts_a,
        label=labels[0],
        color=style_a["color"],
        linestyle=style_a["linestyle"],
        marker=style_a["marker"],
        linewidth=linewidth,
        markersize=markersize,
    )

    plot_distribution(
        ax,
        centers_b,
        counts_b,
        label=labels[1],
        color=style_b["color"],
        linestyle=style_b["linestyle"],
        marker=style_b["marker"],
        linewidth=linewidth,
        markersize=markersize,
    )

    decorate_axis(ax, xlim, show_legend=show_legend)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    plt, np, pd = import_plotting_modules()
    setup_plot_style(plt)

    try:
        fitness_a = load_fitness_values(pd, args.csv_a)
        fitness_b = load_fitness_values(pd, args.csv_b)
    except FileNotFoundError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    if len(fitness_a) == 0 or len(fitness_b) == 0:
        print("[ERROR] One of the CSV files has no valid fitness values.", file=sys.stderr)
        return 1

    main_xlim = resolve_xlim(args.preset, fitness_a, fitness_b, args.xlim)
    style_a, style_b = get_styles(args.color)

    fig, ax = plt.subplots(figsize=(10, 7))

    plot_panel(
        np,
        ax,
        fitness_a,
        fitness_b,
        main_xlim,
        args.bins,
        args.labels,
        style_a,
        style_b,
        show_legend=True,
        linewidth=2.4,
        markersize=6.0,
    )

    if args.preset == "inset":
        inset_xlim = (float(args.inset_xlim[0]), float(args.inset_xlim[1]))

        # [left, bottom, width, height] in axes fraction
        axins = ax.inset_axes([0.18, 0.28, 0.52, 0.50])

        plot_panel(
            np,
            axins,
            fitness_a,
            fitness_b,
            inset_xlim,
            args.bins,
            args.labels,
            style_a,
            style_b,
            show_legend=False,
            linewidth=2.2,
            markersize=6.0,
        )

        axins.tick_params(labelsize=18)
        axins.xaxis.label.set_size(22)
        axins.yaxis.label.set_size(22)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] saved: {args.output}")
    print(f"[OK] preset: {args.preset}")
    print(f"[OK] csv-a: {args.csv_a}")
    print(f"[OK] csv-b: {args.csv_b}")
    print(f"[OK] main xlim: {main_xlim}")

    if args.preset == "inset":
        print(f"[OK] inset xlim: {tuple(args.inset_xlim)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())