#!/usr/bin/env python3
"""
invalid_count (!=0) のヒストグラムを作成するスクリプト

Usage:
    python3 invalid_count_histogram.py <stats_csv_path> [--output <output.png>]

Example:
    python3 invalid_count_histogram.py ../evo_sim/results/ver1/sigma_0.020/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/stats/all_individuals_fitness_stats_sigma=0.020.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(
        description='Create histogram for invalid_count (!=0) from fitness stats CSV'
    )
    parser.add_argument(
        'stats_csv',
        help='Path to all_individuals_fitness_stats CSV file'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output PNG file path (default: same directory as input with _invalid_histogram.png suffix)'
    )
    parser.add_argument(
        '--bin_width',
        type=int,
        default=50,
        help='Histogram bin width (default: 50)'
    )
    args = parser.parse_args()

    # CSVを読み込み
    stats_csv = Path(args.stats_csv)
    if not stats_csv.exists():
        print(f"[ERROR] File not found: {stats_csv}")
        return 1

    df = pd.read_csv(stats_csv)
    
    # invalid_count列があるか確認
    if 'invalid_count' not in df.columns:
        print(f"[ERROR] 'invalid_count' column not found in {stats_csv}")
        print(f"Available columns: {', '.join(df.columns)}")
        return 1
    
    # unique_fitness列があるか確認
    has_unique_fitness = 'unique_fitness' in df.columns

    # sigma値を抽出
    sigma = None
    m = re.search(r'sigma[_=]([\d.]+)', str(stats_csv))
    if m:
        sigma = float(m.group(1))

    # 出力ファイル名を決定
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = stats_csv.parent / (stats_csv.stem + '_invalid_count_histogram.png')

    # invalid_count != 0 のデータを抽出
    data = df['invalid_count'].dropna()
    data_nonzero = data[data != 0]
    
    if len(data_nonzero) == 0:
        print(f"[ERROR] No invalid_count != 0 data found")
        return 1

    # ヒストグラムと累積分布関数を作成
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # binの範囲を設定（0から10000まで、bin_width刻み）
    bins = np.arange(0, 10000 + args.bin_width, args.bin_width)
    
    # 左側：ヒストグラム
    ax1.hist(data_nonzero, bins=bins, color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('invalid_count')
    ax1.set_ylabel('Frequency')
    ax1.set_xlim([-200, 10200])
    
    title = f'Histogram of invalid_count (!=0, bin width={args.bin_width})'
    if sigma is not None:
        title += f'\n(sigma={sigma:.3f})'
    ax1.set_title(title)
    ax1.grid(axis='y', alpha=0.3)
    
    # 右側：累積分布関数
    sorted_data = np.sort(data_nonzero)
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    
    ax2.plot(sorted_data, cdf, color='steelblue', linewidth=2)
    ax2.set_xlabel('invalid_count')
    ax2.set_ylabel('Cumulative Probability')
    
    cdf_title = f'CDF of invalid_count (!=0)'
    if sigma is not None:
        cdf_title += f'\n(sigma={sigma:.3f})'
    ax2.set_title(cdf_title)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1])
    
    # 統計情報をテキストボックスで表示
    stats_text = (
        f"Total individuals: {len(data)}\n"
        f"invalid_count == 0: {(data == 0).sum()}\n"
        f"invalid_count != 0: {len(data_nonzero)}\n"
        f"\n"
        f"Mean (!=0): {data_nonzero.mean():.2f}\n"
        f"Median (!=0): {data_nonzero.median():.2f}\n"
        f"Min (!=0): {data_nonzero.min():.2f}\n"
        f"Max (!=0): {data_nonzero.max():.2f}\n"
        f"Q1 (!=0)(25%): {data_nonzero.quantile(0.25):.2f}\n"
        f"Q3 (!=0)(75%): {data_nonzero.quantile(0.75):.2f}"
    )
    
    ax2.text(0.98, 0.02, stats_text,
             transform=ax2.transAxes,
             fontsize=8,
             verticalalignment='bottom',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    plt.tight_layout()
    
    # 保存
    plt.savefig(output_path, dpi=200)
    plt.close()
    
    print(f"[INFO] Saved histogram to: {output_path}")
    
    # 統計情報を表示
    print(f"\n=== invalid_count (!=0) statistics ===")
    print(f"Total individuals: {len(data)}")
    print(f"invalid_count == 0: {(data == 0).sum()}")
    print(f"invalid_count != 0: {len(data_nonzero)}")
    print(f"Mean (!=0): {data_nonzero.mean():.2f}")
    print(f"Median (!=0): {data_nonzero.median():.2f}")
    print(f"Min (!=0): {data_nonzero.min():.2f}")
    print(f"Max (!=0): {data_nonzero.max():.2f}")
    print(f"Q1 (25%): {data_nonzero.quantile(0.25):.2f}")
    print(f"Q3 (75%): {data_nonzero.quantile(0.75):.2f}")

    return 0


if __name__ == '__main__':
    exit(main())
