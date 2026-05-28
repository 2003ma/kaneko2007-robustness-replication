# calc_all_J_x_noise

Evaluate the fitness of all evolved individuals under varying levels of simulation noise.

## Overview

This program evaluates the gene expression dynamics robustness of evolved regulatory networks by simulating them under different noise conditions. For each individual in an evolved population (at a specified generation), it runs multiple trials of stochastic gene expression dynamics with increasing simulation noise levels, and records the fitness values.

This is used to assess how the robustness of evolved networks varies with the noise environment during evaluation.

## Source files

- `all_J_trial_eval_noise.c`, `all_J_trial_eval_noise.h`: Core implementation for evaluating all individuals under varying noise levels.
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
./bin/all_J_trial_eval_noise \
  --input <input_dir> \
  --evo-sigma <sigma> \
  [options]
```

### Required arguments

- `--input <dir>`: Path to the directory containing evolved individuals (e.g., `../../../evo_sim/results/ver1/sigma_0.005`).
- `--evo-sigma <sigma>`: The noise level used during evolution (e.g., `0.005`, `0.040`, `0.100`).

### Optional arguments

- `--sim-sigma-start <N>`: Simulation noise start value (default: `0.0`).
- `--sim-sigma-end <N>`: Simulation noise end value (default: `0.4`).
- `--sim-sigma-step <N>`: Simulation noise step size (default: `0.01`).
- `--generation <N>`: Generation to evaluate (default: `200`).
- `--trials <N>`: Number of trials per individual per noise level (default: `10000`).
- `--t1 <N>`: Transient period start (default: `80`).
- `--t2 <N>`: Transient period end (default: `90`).
- `--beta <N>`: Inverse temperature for selection (default: `7.0`).
- `--dt <N>`: Time step (default: `0.05`).
- `--k-boundary <N>`: Gene ID boundary for tracking (default: `8`).
- `--seed <N>`: Random seed (default: `12345`).
- `--start_ind <N>`: Start index of individuals to evaluate (default: `0`).
- `--end_ind <N>`: End index of individuals to evaluate (inclusive, default: `1` = first individual).
- `--jemk <N>`: JEMK flag (default: `1`).

### Example

```bash
./bin/all_J_trial_eval_noise \
  --input ../../../evo_sim/results/ver1/sigma_0.005 \
  --evo-sigma 0.005 \
  --sim-sigma-start 0.0 \
  --sim-sigma-end 0.4 \
  --sim-sigma-step 0.01
```

## Output

The program generates CSV files in the input directory structure, organizing results by simulation noise levels. The output format contains fitness values for all individuals across the requested noise conditions.

## Notes

- The program uses OpenMP for parallel evaluation. The default number of threads is 20 (configurable via `OMP_NUM_THREADS`).
- Compilation requires libomp (or equivalent OpenMP library) and is configured for clang on macOS with Apple Silicon optimization flags.
- This is a research script designed for batch evaluation of many individuals; adjust `--trials` and `OMP_NUM_THREADS` based on available computational resources.
