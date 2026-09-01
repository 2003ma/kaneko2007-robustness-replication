# evolved_J_evaluation

Evaluation and analysis of evolved regulatory networks across individuals and noise conditions.

## Overview

This directory contains tools and scripts to evaluate the fitness of all individuals evolved by `evo_sim`, and to analyze the resulting fitness statistics. The evaluation can be performed under:

1. **No simulation noise** (`calc_all_J_x`): Deterministic evaluation where networks are tested under noise-free conditions.
2. **Varying simulation noise levels** (`calc_all_J_x_noise`): Stochastic evaluation where networks are tested under different noise strengths to assess their robustness.

The Python script `analyze_all_J_fitness.py` then aggregates the trial results to produce comprehensive fitness statistics for each individual.

## Directory structure

```
evolved_J_evaluation/
├── analyze_all_J_fitness.py          # Python script to aggregate fitness statistics
├── src/
│   ├── calc_all_J_x/                 # C program: evaluate all individuals (no noise)
│   │   ├── Makefile
│   │   ├── main.c
│   │   ├── all_J_trial_eval.c
│   │   ├── all_J_trial_eval.h
│   │   └── bin/
│   ├── calc_all_J_x_noise/           # C program: evaluate all individuals (with noise)
│   │   ├── Makefile
│   │   ├── main.c
│   │   ├── all_J_trial_eval_noise.c
│   │   ├── all_J_trial_eval_noise.h
│   │   └── bin/
│   └── core/                          # Shared C utilities
│       ├── dynamics.c
│       ├── dynamics.h
│       ├── j_loader.c
│       └── j_loader.h
└── README.md
```

## C Programs

### calc_all_J_x: Deterministic Evaluation

Evaluates all individuals in an evolved population under **no simulation noise** (deterministic dynamics). Outputs fitness values for each individual from a specified generation.

**Initial conditions**: Random initial gene expression state.

See [src/calc_all_J_x/](src/calc_all_J_x/) for details.

**Usage example:**
```bash
cd src/calc_all_J_x
make
./bin/all_J_trial_eval --input ../../../evo_sim/results/ver1/sigma_0.400
```

### calc_all_J_x_noise: Stochastic Evaluation with Varying Noise

Evaluates all individuals under **varying simulation noise levels** to assess robustness. For each noise condition and each individual, runs multiple trials of stochastic gene expression dynamics.

**Initial conditions**: Fixed initial gene expression state (same as used during evolution).

See [src/calc_all_J_x_noise/](src/calc_all_J_x_noise/) for details.

**Usage example:**
```bash
cd src/calc_all_J_x_noise
make
./bin/all_J_trial_eval_noise \
  --input ../../../evo_sim/results/ver1/sigma_0.005 \
  --evo-sigma 0.005 \
  --sim-sigma-start 0.0 \
  --sim-sigma-end 0.4 \
  --sim-sigma-step 0.01
```

### Core Dynamics

Shared C code located in `src/core/`:

- `dynamics.c`, `dynamics.h`: Implementation of gene regulatory network dynamics (Euler-Maruyama stochastic integration and Runge-Kutta deterministic integration).
- `j_loader.c`, `j_loader.h`: Functions to load regulatory matrices (J) and individual metadata from the evolved population files.

## Python Analysis: analyze_all_J_fitness.py

This script aggregates the trial fitness results from a single run (`calc_all_J_x` or `calc_all_J_x_noise`) and computes comprehensive fitness statistics for each individual.

### Purpose

After running `calc_all_J_x` or `calc_all_J_x_noise`, the output is organized as individual CSV files, each containing fitness values from multiple trials. This script:

1. Reads all trial fitness CSV files from a specified directory.
2. Computes statistical measures for each individual's fitness distribution.
3. Outputs a single CSV file summarizing statistics across all individuals.

### Fitness Classification

The script classifies fitness values as:

- **Valid**: Fitness = 0 or -8 (network converges to a stable state)
- **Invalid**: Fitness ≠ 0 and ≠ -8 (network fails to converge or exhibits chaotic behavior)

Invalid fitness values are analyzed separately to understand failure modes.

### Output Statistics

For each individual, the script computes:

**Overall statistics:**
- `total_trials`: Number of trials run
- `unique_fitness`: Count of unique exact fitness values
- `unique_fitness_binned_0.1`: Count of unique fitness values after binning to 0.1 intervals
- `mean`, `std`: Mean and standard deviation of all fitness values
- `min`, `max`, `median`: Minimum, maximum, and median fitness
- `q25`, `q75`: 25th and 75th percentiles

**Invalid-only statistics:**
- `invalid_count`: Number of invalid fitness trials
- `invalid_ratio`: Proportion of invalid trials (invalid_count / total_trials)
- `unique_invalid_fitness_binned_0.1`: Count of unique invalid fitness values (binned)
- `invalid_mean`, `invalid_std`: Mean and standard deviation of invalid fitness values only
- `invalid_min`, `invalid_max`, `invalid_median`: Extreme and median values for invalid trials
- `invalid_q25`, `invalid_q75`: Quartiles for invalid trials

**Metadata:**
- `individual_id`: Individual index
- `sigma`: Evolutionary noise level (extracted from path)

### Usage

```bash
python3 analyze_all_J_fitness.py <trial_fitness_dir> [options]
```

**Positional arguments:**
- `<trial_fitness_dir>`: Path to the directory containing trial fitness CSV files (e.g., `../evo_sim/results/ver1/sigma_0.400/.../.../allx/`)

**Optional arguments:**
- `--output <filename>`: Output CSV path (default: `<trial_fitness_dir>/all_individuals_fitness_stats.csv`)
- `--start_ind <N>`: Only analyze individuals with ID ≥ N (default: 0)
- `--end_ind <N>`: Only analyze individuals with ID ≤ N (inclusive, default: analyze all)

### Examples

**Analyze all individuals from a trial directory:**
```bash
python3 analyze_all_J_fitness.py ../evo_sim/results/ver1/sigma_0.400/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/allx/
```

**Analyze individuals 0-99 and save to a custom output file:**
```bash
python3 analyze_all_J_fitness.py \
  ../evo_sim/results/ver1/sigma_0.400/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/allx/ \
  --end_ind 99 \
  --output custom_stats.csv
```

**Analyze individuals 50-150:**
```bash
python3 analyze_all_J_fitness.py \
  ../evo_sim/results/ver1/sigma_0.400/trials10000/t1=80_t2=90_trials10000_dt0.050_allJ/gen_200_all_individuals/allx/ \
  --start_ind 50 \
  --end_ind 150
```

### Output Format

The script outputs a CSV file with one row per individual and columns as listed above. Example columns:

```
individual_id,sigma,total_trials,invalid_count,invalid_ratio,unique_fitness,...,invalid_q75
0,0.400,10000,5000,0.5,42,15,3.14,1.5,-10.2,5.0,...,4.5
1,0.400,10000,4800,0.48,43,16,3.20,1.6,-10.5,5.1,...,4.6
...
```

### Notes

- The script automatically extracts the `sigma` value from the directory path using regex pattern matching.
- File names must follow the pattern: `trial_fitness_gen=<gen>_ind=<id>_..._allx.csv`
- If a sigma value cannot be extracted from the path, the `sigma` column will be omitted.

## Notes

- Both C programs use OpenMP for parallel evaluation. Adjust `OMP_NUM_THREADS` for your system.
- The evaluation is memory-intensive due to storing all trial results. Consider evaluating subsets of individuals (using `--start_ind`/`--end_ind`) if memory is limited.
- These scripts are research implementations; optimization and features may be added as needed.
