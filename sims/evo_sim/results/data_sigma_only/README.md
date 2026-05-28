# data_sigma_only

## Overview

This directory contains simulation data organized for the analysis scripts in `sims/analysis_paper/`. 
It holds a simplified version of the full results, structured by sigma value for easy reference.

## Directory Structure

- `ver1/` — First experiment set
  - Contains subdirectories organized by sigma value: `sigma_0.005/`, `sigma_0.010/`, etc.
  - Each `sigma_xxx/` directory houses an `evo_sim_data/` folder
  - CSV files with the exact sigma value parameter are placed inside

## Example Layout

```
data_sigma_only/
├── README.md
└── ver1/
    ├── sigma_0.005/
    │   └── evo_sim_data/
    │       └── data,sigma=0.0050,dt=0.005,pop=300,L=300,meas=2000,relax=16000,seed=123456789.csv
    ├── sigma_0.010/
    │   └── evo_sim_data/
    │       └── data,sigma=0.0100,dt=0.005,pop=300,L=300,meas=2000,relax=16000,seed=123456789.csv
    └── ...
```

## Usage

The analysis scripts in `sims/analysis_paper/` reference data from this directory structure.

Make sure CSV files are placed in the corresponding `sigma_xxx/evo_sim_data/` folders with the correct sigma value in the filename.

This repository does not include a script for copying or reorganizing the CSV files into this directory. Please manually place the required files in the appropriate locations, following the structure shown above.