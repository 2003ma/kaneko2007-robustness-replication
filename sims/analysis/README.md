# Analysis Scripts

This directory contains utility scripts for visualizing simulation statistics from the evolutionary simulation results.
Created statistics is used in TABLE 2.

## Scripts

### `invalid_count_histogram.py`
Generates a histogram of invalid fitness counts (constraints violations).

**Purpose:** Analyzes how many constraint violations occur during evolution at different sigma (mutation) levels.

**Usage:**
```bash
python3 invalid_count_histogram.py <path_to_stats_csv> [--output output.png] [--bin_width 10]
```

**Arguments:**
- `stats_csv` — Path to the fitness statistics CSV file (required)
- `--output` — Output PNG file path (default: same directory as input with `_invalid_count_histogram.png` suffix)
- `--bin_width` — Histogram bin width (default: 50)

**Example:**
```bash
python3 invalid_count_histogram.py ../evo_sim/results/ver1/sigma_0.020/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/stats/all_individuals_fitness_stats_sigma=0.020.csv
```

**Output:** Histogram + cumulative distribution function of invalid_count values

---

### `unique_fitness_histogram.py`
Generates a histogram of unique fitness values (with 0.1 binning).

**Purpose:** Visualizes the distribution of fitness diversity in the evolved population, optionally filtered by constraint violations.

**Usage:**
```bash
python3 unique_fitness_histogram.py <path_to_stats_csv> [--output output.png] [--bin_width 10]
```

**Arguments:**
- `stats_csv` — Path to the fitness statistics CSV file (required)
- `--output` — Output PNG file path (default: same directory as input with `_unique_fitness_histogram.png` suffix)
- `--bin_width` — Bin width for histogram (default: 10)

**Example:**
```bash
python3 unique_fitness_histogram.py ../evo_sim/results/ver1/sigma_0.040/trials10000/.../stats/all_individuals_fitness_stats_sigma=0.040.csv
```

**Output:** Histogram + cumulative distribution function of unique_fitness_binned_0.1 values (filtered to invalid_count != 0 if available)

---

## Dependencies

- `pandas` — CSV reading and data manipulation
- `matplotlib` — Plot generation
- `numpy` — Numerical operations and binning

## Notes

- Sigma values are automatically extracted from the CSV filename where possible
- The scripts handle missing or malformed data gracefully with error messages
- Both scripts produce side-by-side visualizations: histogram on the left, cumulative distribution on the right
