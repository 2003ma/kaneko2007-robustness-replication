# kaneko2007-robustness-replication

This repository contains code used in a replication study of the evolutionary gene regulatory network model proposed by Kaneko (2007).

The code was prepared as part of a university bulletin manuscript. It is mainly intended to record the implementation used in the study and to make the simulation code publicly available.

## Original paper

Kaneko, K. (2007).  
*Evolution of Robustness to Noise and Mutation in Gene Expression Dynamics.*  
PLoS ONE, 2(5), e434.  
https://doi.org/10.1371/journal.pone.0000434

## Contents

This repository includes code for:

- simulating stochastic gene expression dynamics,
- evolving gene regulatory networks under mutation and selection,
- calculating fitness values and phenotypic variances,
- generating plots used in the manuscript.

## Requirements

The code was developed with Python 3.

Required packages are listed in `requirements.txt` if available.

## Notes

This repository is not intended to be a fully packaged software project.  
The main purpose is to archive and share the code used in the replication study.

Some results may depend on random seeds and parameter settings.

## License

This repository is released under the MIT License.

## Tree
```
.
├── LICENSE
├── README.md
└── sims
    ├── analysis
    │   ├── README.md
    │   ├── invalid_count_histogram.py --TABLE 2 IC (use analyze_all_J_fitness.py results)
    │   └── unique_fitness_histogram.py --TABLE 2 UF
    ├── analysis_paper
    │   ├── README.md
    │   ├── fitness --Fugre 1 (use evo_sim results)
    │   ├── fitness_dist_random_initial --Figure 6 (use calc_all_J_x results)
    │   ├── fitness_last_gen --Figure 3 (use evo_sim results)
    │   ├── overview --Figure 2 (use evo_sim results)
    │   ├── under_noise_fitness --Figure 5 (use calc_all_J_x_nois results)
    │   └── vg_vip --Figure 4 (use evo_sim results)
    ├── evo_sim
    │   ├── Makefile
    │   ├── README.md
    │   ├── bin
    │   ├── main
    │   ├── results
    │   └── src
    └── evolved_J_evaluation
            ├── README.md
            ├── analyze_all_J_fitness.py --TABLE 2(use calc_all_J_x results)
            └── src
                ├── calc_all_J_x (use evo_sim results)
                ├── calc_all_J_x_noise (use evo_sim results)
                └── core
```