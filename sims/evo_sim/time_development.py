"""
time_development.py

世代ごとの個体データCSVから任意の個体のJを読み込み、
x_i(0) = 0 から遺伝子発現ダイナミクスを計算する。

入力CSVの形式:
    id, fitness, v_ip, J_0_0, J_0_1, ..., J_63_63

出力:
    target_gene_trajectory.csv
    target_gene_trajectory.png

実行例:
    python3 time_development.py \
        --csv results/ver1/sigma_0.200/evo_sim_data/gen_0_all_J_sigma_0.200_dt0.005.csv \
        --individual-id 0 \
        --sigma 0.2
        --output-figure ver1/sigma_0.200/target_gene_trajectory_gen40_id0.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "CSVから1個体の遺伝子制御行列Jを読み込み、"
            "x_i(0)=0から発現ダイナミクスを計算する。"
        )
    )

    parser.add_argument(
        "--csv",
        type=Path,
        required=True,
        help="Jを保存した個体データCSV",
    )

    parser.add_argument(
        "--individual-id",
        type=int,
        default=0,
        help="使用する個体のid。デフォルトは0",
    )

    parser.add_argument(
        "--N",
        type=int,
        default=64,
        help="全遺伝子数。デフォルトは64",
    )

    parser.add_argument(
        "--k",
        type=int,
        default=8,
        help="ターゲット遺伝子数。デフォルトは8",
    )

    parser.add_argument(
        "--beta",
        type=float,
        default=7.0,
        help="tanhの応答係数。デフォルトは7.0",
    )

    parser.add_argument(
        "--sigma",
        type=float,
        default=0.2,
        help="発現ノイズ強度。デフォルトは0.2",
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=0.005,
        help="時間刻み。デフォルトは0.005",
    )

    parser.add_argument(
        "--t-relax",
        type=float,
        default=80.0,
        help="緩和時間。デフォルトは80",
    )

    parser.add_argument(
        "--t-meas",
        type=float,
        default=60.0,
        help="測定時間。デフォルトは60",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="乱数シード",
    )

    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("target_gene_trajectory.csv"),
        help="時系列CSVの出力先",
    )

    parser.add_argument(
        "--output-figure",
        type=Path,
        default=Path("target_gene_trajectory.png"),
        help="図の出力先",
    )

    return parser.parse_args()


def load_individual_J(
    csv_path: Path,
    individual_id: int,
    N: int,
) -> tuple[np.ndarray, pd.Series]:
    """
    CSVから指定したidの個体を読み込み、
    JをN×N行列として復元する。
    """

    if not csv_path.exists():
        raise FileNotFoundError(f"CSVが見つかりません: {csv_path}")

    df = pd.read_csv(csv_path)

    required_metadata = {"id", "fitness", "v_ip"}

    missing_metadata = required_metadata - set(df.columns)

    if missing_metadata:
        raise ValueError(
            f"CSVに必要な列がありません: {sorted(missing_metadata)}"
        )

    matched = df.loc[df["id"] == individual_id]

    if matched.empty:
        available_ids = df["id"].tolist()

        raise ValueError(
            f"id={individual_id} の個体が見つかりません。\n"
            f"利用可能なidの例: {available_ids[:20]}"
        )

    if len(matched) > 1:
        print(
            f"[WARNING] id={individual_id} が複数行あります。"
            "最初の行を使用します。"
        )

    row = matched.iloc[0]

    J = np.empty((N, N), dtype=np.int8)

    for i in range(N):
        for j in range(N):
            column = f"J_{i}_{j}"

            if column not in df.columns:
                raise ValueError(
                    f"Jの列が不足しています: {column}"
                )

            value = int(row[column])

            if value not in (-1, 0, 1):
                raise ValueError(
                    f"{column}={value} ですが、"
                    "Jの要素は-1, 0, 1である必要があります。"
                )

            J[i, j] = value

    return J, row


def simulate_gene_expression(
    J: np.ndarray,
    *,
    k: int,
    beta: float,
    sigma: float,
    dt: float,
    t_relax: float,
    t_meas: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    固定したJについて、x_i(0)=0から発現ダイナミクスを計算する。

    モデル:
        dx_i/dt
        = -x_i
          + tanh(beta * sum_{j=k}^{N-1} J_ij x_j)
          + sigma * eta_i(t)

    Euler-Maruyama法:
        x_i(t+dt)
        = x_i(t)
          + [-x_i + tanh(beta*h_i)]dt
          + sigma*sqrt(dt)*N(0,1)
    """

    N = J.shape[0]

    if J.shape != (N, N):
        raise ValueError("Jは正方行列である必要があります。")

    if not 0 <= k < N:
        raise ValueError("kは0以上N未満である必要があります。")

    if dt <= 0:
        raise ValueError("dtは正である必要があります。")

    total_time = t_relax + t_meas
    total_steps = int(round(total_time / dt))

    times = np.arange(total_steps + 1, dtype=float) * dt

    trajectory = np.empty(
        (total_steps + 1, N),
        dtype=float,
    )

    # 全遺伝子の初期値を-1にする
    x = np.full(N, -1.0, dtype=float)
    trajectory[0] = x

    rng = np.random.default_rng(seed)

    sqrt_dt = np.sqrt(dt)

    # 元コードと同様に、入力側にはj >= kのみを使う
    J_internal = J[:, k:]

    for step in range(1, total_steps + 1):
        h = J_internal @ x[k:]

        drift = -x + np.tanh(beta * h)

        noise = (
            sigma
            * sqrt_dt
            * rng.standard_normal(N)
        )

        x = x + drift * dt + noise

        trajectory[step] = x

    return times, trajectory


def save_target_trajectory(
    output_path: Path,
    times: np.ndarray,
    trajectory: np.ndarray,
    k: int,
    t_relax: float,
) -> None:
    """
    ターゲット遺伝子x_0,...,x_{k-1}をCSVに保存する。
    """

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    data: dict[str, np.ndarray] = {
        "time": times,
    }

    data["phase"] = np.where(
        times <= t_relax,
        "relaxation",
        "measurement",
    )

    for i in range(k):
        data[f"x_{i}"] = trajectory[:, i]

    target_values = trajectory[:, :k]

    data["number_of_on_targets"] = np.sum(
        target_values > 0,
        axis=1,
    )

    result_df = pd.DataFrame(data)

    result_df.to_csv(
        output_path,
        index=False,
    )


def plot_target_trajectory(
    output_path: Path,
    times: np.ndarray,
    trajectory: np.ndarray,
    k: int,
    t_relax: float,
    individual_id: int,
    fitness: float,
    sigma: float,
) -> None:
    """
    8つのターゲット遺伝子の時間発展を描画する。
    """

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig, ax = plt.subplots(
        figsize=(11, 6),
    )

    for i in range(k):
        ax.plot(
            times,
            trajectory[:, i],
            linewidth=1.0,
            label=rf"$x_{{{i}}}$",
        )

    ax.axhline(
        0.0,
        linestyle=":",
        linewidth=1.0,
    )

    ax.axvline(
        t_relax,
        linestyle="--",
        linewidth=1.2,
        label=rf"$t_{{\mathrm{{relax}}}}={t_relax:g}$",
    )

    ax.set_xlabel("Time")
    ax.set_ylabel("Gene expression level")

    ax.set_title(
        "Target-gene dynamics\n"
        f"id={individual_id}, "
        f"fitness={fitness:.6f}, "
        rf"$\sigma={sigma:g}$"
    )

    ax.legend(
        ncol=3,
        fontsize=9,
    )

    ax.grid(alpha=0.25)

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


def main() -> None:
    args = parse_args()

    J, individual = load_individual_J(
        csv_path=args.csv,
        individual_id=args.individual_id,
        N=args.N,
    )

    fitness = float(individual["fitness"])
    v_ip = float(individual["v_ip"])

    print("Selected individual")
    print(f"  id      = {args.individual_id}")
    print(f"  fitness = {fitness}")
    print(f"  v_ip    = {v_ip}")
    print(
        f"  nonzero elements of J = "
        f"{np.count_nonzero(J)}"
    )

    times, trajectory = simulate_gene_expression(
        J,
        k=args.k,
        beta=args.beta,
        sigma=args.sigma,
        dt=args.dt,
        t_relax=args.t_relax,
        t_meas=args.t_meas,
        seed=args.seed,
    )

    save_target_trajectory(
        output_path=args.output_csv,
        times=times,
        trajectory=trajectory,
        k=args.k,
        t_relax=args.t_relax,
    )

    plot_target_trajectory(
        output_path=args.output_figure,
        times=times,
        trajectory=trajectory,
        k=args.k,
        t_relax=args.t_relax,
        individual_id=args.individual_id,
        fitness=fitness,
        sigma=args.sigma,
    )

    print("\nFinal target-gene states")

    for i in range(args.k):
        state = "ON" if trajectory[-1, i] > 0 else "OFF"

        print(
            f"  x_{i} = "
            f"{trajectory[-1, i]: .6f} "
            f"({state})"
        )

    print(f"\nCSV saved to: {args.output_csv}")
    print(f"Figure saved to: {args.output_figure}")


if __name__ == "__main__":
    main()