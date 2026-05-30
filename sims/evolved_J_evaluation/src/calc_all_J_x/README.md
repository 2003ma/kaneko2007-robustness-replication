# calc_all_J_x

Evaluate the fitness of all evolved individuals under deterministic (noise-free) conditions.

## Overview

This program evaluates all individuals in an evolved population by simulating their gene expression dynamics without any noise (deterministic conditions). For each individual at a specified generation, it runs multiple trials with random initial gene expression states and records the fitness values.

**Initial conditions**: The initial gene expression state is **random and regenerated for each trial** (using a different random seed). This allows assessment of the network's capacity to reach stable states from various starting states despite the deterministic dynamics.


## Source files

- `all_J_trial_eval.c`, `all_J_trial_eval.h`: Core implementation for evaluating all individuals without simulation noise.
- `main.c`: Command-line interface and argument parsing.
- `Makefile`: Build configuration.

The program links against shared core dynamics code in `../core/dynamics.c`.

## Build

To compile:

```bash
make
```

To clean build artifacts:

```bash
make clean
```

## Usage

```bash
./bin/all_J_trial_eval \
  --input <input_dir> \
  [options]
```

### Required arguments

- `--input <dir>`: Path to the directory containing evolved individuals (e.g., `../../../evo_sim/results/ver1/sigma_0.400`).

### Optional arguments

- `--generation <N>`: Generation to evaluate (default: `200`).
- `--trials <N>`: Number of trials per individual (default: `10000`).
- `--t1 <N>`: Transient period start (default: `80`).
- `--t2 <N>`: Transient period end (default: `90`).
- `--beta <N>`: Inverse temperature for selection (default: `7.0`).
- `--dt <N>`: Time step (default: `0.05`).
- `--k-boundary <N>`: Gene ID boundary for tracking (default: `8`).
- `--seed <N>`: Random seed (default: `12345`).
- `--sigma <N>`: Noise level (for informational purposes, default: `0.0` for deterministic evaluation).
- `--start_ind <N>`: Start index of individuals to evaluate (default: `0`).
- `--jemk <N>`: JEMK flag (default: `1`).

### Example

```bash
./bin/all_J_trial_eval \
  --input ../../../evo_sim/results/ver1/sigma_0.400 
```

## Output

The program generates CSV files in the input directory structure. Each file contains fitness values from multiple trials for a given individual.

## Notes

- The program uses OpenMP for parallel evaluation. The default number of threads is 20 (configurable via `OMP_NUM_THREADS`).
- Compilation requires libomp (or equivalent OpenMP library) and is configured for clang on macOS with Apple Silicon optimization flags.
