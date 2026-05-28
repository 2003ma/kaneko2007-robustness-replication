#!/usr/bin/env python3
"""
unique_fitness_binned_0.1のヒストグラムを作成するスクリプト

Usage:
    python3 unique_fitness_histogram.py <stats_csv_path> [--output <output.png>] [--bin_width <width>]

Example:
    python3 unique_fitness_histogram.py ../evo_sim/results/ver1/sigma_0.040/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/stats/all_individuals_fitness_stats_sigma=0.040.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(
        description='Create histogram for unique_fitness_binned_0.1 from fitness stats CSV'
    )
    parser.add_argument(
        'stats_csv',
        help='Path to all_individuals_fitness_stats CSV file'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output PNG file path (default: same directory as input with _histogram.png suffix)'
    )

    parser.add_argument(
        '--bin_width',
        type=int,
        default=10,
        help='Bin width for histogram (default: 10)'
    )
    args = parser.parse_args()

    # CSVを読み込み（入力パスをそのまま受け取り、~も展開）
    stats_csv = Path(args.stats_csv).expanduser()
    if not stats_csv.exists():
        print(f"[ERROR] File not found: {stats_csv}")
        return 1

    if stats_csv.suffix.lower() != '.csv':
        print(f"[WARN] Input file does not look like CSV: {stats_csv}")

    try:
        df = pd.read_csv(stats_csv)
    except Exception as e:
        print(f"[ERROR] Failed to read CSV: {stats_csv}")
        print(f"Reason: {e}")
        return 1
    
    target_col = 'unique_fitness_binned_0.1'

    # 集計対象列があるか確認
    if target_col not in df.columns:
        print(f"[ERROR] '{target_col}' column not found in {stats_csv}")
        print(f"Available columns: {', '.join(df.columns)}")
        return 1

    # sigma値を抽出
    sigma = None
    m = re.search(r'sigma[_=]([\d.]+)', str(stats_csv))
    if m:
        sigma = float(m.group(1))

    # 出力ファイル名を決定
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = stats_csv.parent / (stats_csv.stem + '_unique_fitness_histogram.png')

    # NaNを除外
    data = df[target_col].dropna()

    # invalid_count列があれば、invalid_count != 0のデータのみに絞る
    if 'invalid_count' in df.columns:
        mask_nonzero = df['invalid_count'] != 0
        data = df.loc[mask_nonzero, target_col].dropna()
        filter_label = "invalid_count!=0"
    else:
        filter_label = "all"
    
    if len(data) == 0:
        print(f"[ERROR] No valid unique_fitness data found")
        return 1

    # ヒストグラムと累積分布関数を作成
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # binの範囲を設定（0からmax値まで、bin_width刻み）
    max_val = data.max()
    min_val = 0
    bins = np.arange(min_val, max_val + args.bin_width, args.bin_width)
    
    # 左側：ヒストグラム
    ax1.hist(data, bins=bins, color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel(target_col)
    ax1.set_ylabel('Frequency')
    
    title = f'Histogram of {target_col} ({filter_label}, bin width={args.bin_width})'
    if sigma is not None:
        title += f'\n(sigma={sigma:.3f})'
    ax1.set_title(title)
    ax1.grid(axis='y', alpha=0.3)
    
    # 右側：累積分布関数
    sorted_data = np.sort(data)
    cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    
    ax2.plot(sorted_data, cdf, color='steelblue', linewidth=2)
    ax2.set_xlabel(target_col)
    ax2.set_ylabel('Cumulative Probability')
    
    cdf_title = f'CDF of {target_col} ({filter_label})'
    if sigma is not None:
        cdf_title += f'\n(sigma={sigma:.3f})'
    ax2.set_title(cdf_title)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1])
    
    # 統計情報をテキストボックスで表示
    stats_text = (
        f"Count: {len(data)}\n"
        f"\n"
        f"Mean: {data.mean():.2f}\n"
        f"Median: {data.median():.2f}\n"
        f"Min: {data.min():.2f}\n"
        f"Max: {data.max():.2f}\n"
        f"Q1 (25%): {data.quantile(0.25):.2f}\n"
        f"Q3 (75%): {data.quantile(0.75):.2f}\n"
        f"IQR: {data.quantile(0.75) - data.quantile(0.25):.2f}"
    )
    
    ax2.text(0.98, 0.02, stats_text,
             transform=ax2.transAxes,
             fontsize=9,
             verticalalignment='bottom',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    plt.tight_layout()

    # 保存
    plt.savefig(output_path, dpi=200)
    plt.close()

    print(f"[INFO] Saved histogram to: {output_path}")
    
    # 統計情報を表示（invalid_count_histogram.pyと同じ形式）
    total_count = len(df)
    
    if 'invalid_count' in df.columns:
        filtered_count = len(data)
        zero_count = total_count - filtered_count
        print(f"\n=== {target_col} ({filter_label}) statistics ===")
        print(f"Total individuals: {total_count}")
        print(f"invalid_count == 0: {zero_count}")
        print(f"invalid_count != 0: {filtered_count}")
    else:
        print(f"\n=== {target_col} statistics ===")
        print(f"Total individuals: {total_count}")
    
    print(f"Mean ({filter_label}): {data.mean():.2f}")
    print(f"Median ({filter_label}): {data.median():.2f}")
    print(f"Min ({filter_label}): {data.min():.2f}")
    print(f"Max ({filter_label}): {data.max():.2f}")
    print(f"Q1 ({filter_label})(25%): {data.quantile(0.25):.2f}")
    print(f"Q3 ({filter_label})(75%): {data.quantile(0.75):.2f}")

    return 0


if __name__ == '__main__':
    exit(main())
