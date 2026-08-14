#!/usr/bin/env python3
"""
Script that aggregates trial_fitness results for all individuals and outputs summary statistics.

Usage:
    python3 analyze_all_J_fitness.py <trial_fitness_dir> [--output <output.csv>]

Example:
    python3 analyze_all_J_fitness.py ../evo_sim/results/ver1/sigma_0.400/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/allx/
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import re


def extract_ind_from_filename(filename):
    """Extract the individual number from the filename"""
    m = re.search(r'ind=(\d+)', filename)
    return int(m.group(1)) if m else None


def analyze_individual_fitness(csv_path):
    """Compute fitness statistics from one individual's CSV file"""
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[WARN] Failed to read {csv_path.name}: {e}")
        return None

    if 'fitness' not in df.columns:
        print(f"[WARN] No 'fitness' column in {csv_path.name}")
        return None

    fitness = df['fitness'].values

    if len(fitness) == 0:
        print(f"[WARN] Empty fitness data in {csv_path.name}")
        return None

    # invalid: fitness != 0 and != -8 (does not converge)
    # valid: fitness == 0 or == -8 (converges to a stable state)
    invalid_fitness = fitness[(fitness != 0) & (fitness != -8)]

    # Bin in steps of 0.1
    fitness_binned = np.round(fitness / 0.1) * 0.1
    invalid_fitness_binned = np.round(invalid_fitness / 0.1) * 0.1

    stats = {
        'total_trials': len(fitness),
        'invalid_count': len(invalid_fitness),
        'invalid_ratio': len(invalid_fitness) / len(fitness) if len(fitness) > 0 else 0,
        'unique_fitness': len(np.unique(fitness)),
        'unique_fitness_binned_0.1': len(np.unique(fitness_binned)),
        'unique_invalid_fitness_binned_0.1': len(np.unique(invalid_fitness_binned)) if len(invalid_fitness) > 0 else 0,
        'mean': np.mean(fitness) if len(fitness) > 0 else np.nan,
        'std': np.std(fitness) if len(fitness) > 1 else np.nan,
        'min': np.min(fitness) if len(fitness) > 0 else np.nan,
        'max': np.max(fitness) if len(fitness) > 0 else np.nan,
        'median': np.median(fitness) if len(fitness) > 0 else np.nan,
        'q25': np.percentile(fitness, 25) if len(fitness) > 0 else np.nan,
        'q75': np.percentile(fitness, 75) if len(fitness) > 0 else np.nan,
    }

    # Compute only when invalid_fitness exists
    if len(invalid_fitness) > 0:
        stats['invalid_mean'] = np.mean(invalid_fitness)
        stats['invalid_std'] = np.std(invalid_fitness) if len(invalid_fitness) > 1 else np.nan
        stats['invalid_min'] = np.min(invalid_fitness)
        stats['invalid_max'] = np.max(invalid_fitness)
        stats['invalid_median'] = np.median(invalid_fitness)
        stats['invalid_q25'] = np.percentile(invalid_fitness, 25)
        stats['invalid_q75'] = np.percentile(invalid_fitness, 75)
    else:
        stats['invalid_mean'] = np.nan
        stats['invalid_std'] = np.nan
        stats['invalid_min'] = np.nan
        stats['invalid_max'] = np.nan
        stats['invalid_median'] = np.nan
        stats['invalid_q25'] = np.nan
        stats['invalid_q75'] = np.nan

    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Analyze all individual fitness from trial_fitness directory'
    )
    parser.add_argument(
        'trial_fitness_dir',
        help='Path to trial_fitness directory'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output CSV file path (default: <trial_fitness_dir>/all_individuals_fitness_stats.csv)'
    )
    parser.add_argument(
        '--start_ind',
        type=int,
        default=0,
        help='Start individual ID (default: 0)'
    )
    parser.add_argument(
        '--end_ind',
        type=int,
        default=None,
        help='End individual ID (inclusive, default: None = all)'
    )

    args = parser.parse_args()

    # Extract sigma value from the path
    trial_fitness_dir = Path(args.trial_fitness_dir)
    sigma = None
    m = re.search(r'sigma[_=]([\d.]+)', str(trial_fitness_dir))
    if m:
        sigma = float(m.group(1))

    if not trial_fitness_dir.exists():
        print(f"[ERROR] Directory not found: {trial_fitness_dir}")
        return 1

    # Collect trial_fitness_gen=XXX_ind=YYY_..._allx.csv files
    csv_files = sorted(list(trial_fitness_dir.glob("trial_fitness_gen=*_ind=*_allx.csv")))

    if not csv_files:
        print(f"[ERROR] No CSV files found in {trial_fitness_dir}")
        return 1

    # Filter files so that only start_ind and later, and up to end_ind, are included
    csv_files = [
        f for f in csv_files
        if extract_ind_from_filename(f.name) is not None
        and extract_ind_from_filename(f.name) >= args.start_ind
    ]

    if args.end_ind is not None:
        csv_files = [
            f for f in csv_files
            if extract_ind_from_filename(f.name) <= args.end_ind
        ]

    if not csv_files:
        end_str = f" and <= {args.end_ind}" if args.end_ind is not None else ""
        print(f"[ERROR] No CSV files found for individual_id >= {args.start_ind}{end_str}")
        return 1

    end_str = f" to {args.end_ind}" if args.end_ind is not None else ""
    print(f"[INFO] Found {len(csv_files)} individual CSV files (ind >= {args.start_ind}{end_str})")

    # Collect statistics for each individual
    results = []
    for csv_path in csv_files:
        ind = extract_ind_from_filename(csv_path.name)
        if ind is None:
            print(f"[WARN] Cannot extract individual ID from {csv_path.name}")
            continue

        stats = analyze_individual_fitness(csv_path)
        if stats is None:
            print(f"[WARN] Failed to analyze {csv_path.name}")
            continue

        stats['individual_id'] = ind
        results.append(stats)

        if (ind + 1) % 50 == 0:
            print(f"[INFO] Processed {ind + 1} individuals...")

    if not results:
        print("[ERROR] No valid individual statistics were generated.")
        return 1

    # Convert to DataFrame
    df_results = pd.DataFrame(results)

    # Sort by individual_id
    df_results = df_results.sort_values('individual_id').reset_index(drop=True)

    # Add sigma column
    if sigma is not None:
        df_results['sigma'] = sigma

    # Reorder columns
    columns_order = [
        'individual_id',
        'sigma',
        'total_trials',
        'invalid_count',
        'invalid_ratio',
        'unique_fitness',
        'unique_fitness_binned_0.1',
        'unique_invalid_fitness_binned_0.1',
        'mean', 'std', 'min', 'max', 'median', 'q25', 'q75',
        'invalid_mean', 'invalid_std', 'invalid_min', 'invalid_max',
        'invalid_median', 'invalid_q25', 'invalid_q75'
    ]
    columns_order = [col for col in columns_order if col in df_results.columns]
    df_results = df_results[columns_order]

    # Create directory for summary statistics (one level above allx)
    stats_dir = trial_fitness_dir.parent / 'stats_all'
    stats_dir.mkdir(exist_ok=True)

    # Determine output filename
    if args.output:
        output_path = Path(args.output)
    else:
        if sigma is not None:
            output_path = stats_dir / f'all_individuals_fitness_stats_sigma={sigma:.3f}.csv'
        else:
            output_path = stats_dir / 'all_individuals_fitness_stats.csv'

    # Merge with existing CSV if present
    if output_path.exists() and args.start_ind > 0:
        print(f"[INFO] Loading existing CSV: {output_path}")
        df_existing = pd.read_csv(output_path)
        # Remove data from start_ind onward and replace it with new data
        df_existing = df_existing[df_existing['individual_id'] < args.start_ind]
        df_results = pd.concat([df_existing, df_results], ignore_index=True)
        df_results = df_results.sort_values('individual_id').reset_index(drop=True)
        print(f"[INFO] Merged with existing data (kept ind < {args.start_ind}, updated ind >= {args.start_ind})")

    # Save per-individual statistics
    df_results.to_csv(output_path, index=False)
    print(f"[INFO] Saved per-individual statistics to: {output_path}")

    # DataFrame for computing statistics across all individuals
    if output_path.exists():
        df_stats_all = pd.read_csv(output_path)
    else:
        df_stats_all = df_results.copy()

    metrics_to_summarize = [
        'invalid_ratio',
        'unique_fitness',
        'unique_fitness_binned_0.1',
        'unique_invalid_fitness_binned_0.1'
    ]

    # Overall stats
    overall_stats = {
        'metric': [],
        'mean': [],
        'std': [],
        'min': [],
        'max': [],
        'median': [],
        'q25': [],
        'q75': []
    }

    for metric in metrics_to_summarize:
        if metric in df_stats_all.columns:
            values = df_stats_all[metric].dropna()
            if len(values) > 0:
                overall_stats['metric'].append(metric)
                overall_stats['mean'].append(np.mean(values))
                overall_stats['std'].append(np.std(values))
                overall_stats['min'].append(np.min(values))
                overall_stats['max'].append(np.max(values))
                overall_stats['median'].append(np.median(values))
                overall_stats['q25'].append(np.percentile(values, 25))
                overall_stats['q75'].append(np.percentile(values, 75))

    # Additional count statistics
    num_invalid0 = (df_stats_all['invalid_ratio'] == 0).sum()
    num_invalid_1digit = ((df_stats_all['invalid_count'] > 0) & (df_stats_all['invalid_count'] < 10)).sum()
    num_invalid_2digit = ((df_stats_all['invalid_count'] >= 10) & (df_stats_all['invalid_count'] < 100)).sum()
    num_invalid_3digit = ((df_stats_all['invalid_count'] >= 100) & (df_stats_all['invalid_count'] < 1000)).sum()
    num_invalid_4digit = (df_stats_all['invalid_count'] >= 1000).sum()
    num_invalid_count_10000 = (df_stats_all['invalid_count'] == 10000).sum()
    num_unique_fitness_gt50 = (df_stats_all['unique_fitness'] > 50).sum()

    extra_stats = pd.DataFrame({
        'metric': [
            'num_invalid_ratio_0',
            'num_invalid_count_1digit',
            'num_invalid_count_2digit',
            'num_invalid_count_3digit',
            'num_invalid_count_4digit',
            'num_invalid_count_10000',
            'num_unique_fitness_gt50',
        ],
        'mean': [
            num_invalid0,
            num_invalid_1digit,
            num_invalid_2digit,
            num_invalid_3digit,
            num_invalid_4digit,
            num_invalid_count_10000,
            num_unique_fitness_gt50
        ],
        'std': [''] * 7,
        'min': [''] * 7,
        'max': [''] * 7,
        'median': [''] * 7,
        'q25': [''] * 7,
        'q75': [''] * 7
    })

    df_overall = pd.concat([pd.DataFrame(overall_stats), extra_stats], ignore_index=True)

    overall_output_path = output_path.parent / (output_path.stem + '_overall_stats.csv')
    df_overall.to_csv(overall_output_path, index=False)
    print(f"[INFO] Saved overall statistics to: {overall_output_path}")

    # Print summary
    print("\n=== Summary (per individual) ===")
    print(f"Invalid ratio range: {df_results['invalid_ratio'].min():.4f} - {df_results['invalid_ratio'].max():.4f}")
    print(f"Mean fitness range: {df_results['mean'].min():.4f} - {df_results['mean'].max():.4f}")

    invalid_mean_values = df_results['invalid_mean'].dropna()
    if len(invalid_mean_values) > 0:
        print(f"Invalid mean fitness range: {invalid_mean_values.min():.4f} - {invalid_mean_values.max():.4f}")
    else:
        print("Invalid mean fitness range: no invalid samples")

    print(f"num_invalid_ratio_0: {num_invalid0}")
    print(f"num_invalid_count_1digit: {num_invalid_1digit}")
    print(f"num_invalid_count_2digit: {num_invalid_2digit}")
    print(f"num_invalid_count_3digit: {num_invalid_3digit}")
    print(f"num_invalid_count_4digit: {num_invalid_4digit}")
    print(f"num_invalid_count_10000: {num_invalid_count_10000}")
    print(f"num_unique_fitness_gt50: {num_unique_fitness_gt50}")

    print("\n=== Overall Statistics (across all individuals) ===")
    print(df_overall.to_string(index=False))

    return 0


if __name__ == '__main__':
    raise SystemExit(main())