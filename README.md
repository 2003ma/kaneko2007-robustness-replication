# kaneko2007-robustness-replication

This repository contains code used in a replication study of the evolutionary gene regulatory network model proposed by Kaneko (2007).  
paper title:A Reexamination of Noise-Driven Robustness Evolution in Gene Regulatory Networks. 
url:https://arxiv.org/abs/2607.18645

## Original paper

Kaneko, K. (2007).  
*Evolution of Robustness to Noise and Mutation in Gene Expression Dynamics.*  
PLoS ONE, 2(5), e434.  
https://doi.org/10.1371/journal.pone.0000434


## License

This repository is released under the MIT License.

## Directory structure
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

## Reproducing figures and tables

Run the following in order for each figure/table.

### Figures

| Output | Run |
|---------|-----|
| **Figure 1** (Fitness transition) | `sims/evo_sim` → `sims/analysis_paper/fitness` |
| **Figure 2** (Overview) | `sims/evo_sim` → `sims/analysis_paper/overview` |
| **Figure 3** (Final-generation fitness) | `sims/evo_sim` → `sims/analysis_paper/fitness_last_gen` |
| **Figure 4** (Vg–Vip analysis) | `sims/evo_sim` → `sims/analysis_paper/vg_vip` |
| **Figure 5** (Fitness under noise) | `sims/evo_sim` → `sims/evolved_J_evaluation/src/calc_all_J_x_noise` → `sims/analysis_paper/under_noise_fitness` |
| **Figure 6** (Random initial fitness distribution) | `sims/evo_sim` → `sims/evolved_J_evaluation/src/calc_all_J_x` → `sims/analysis_paper/fitness_dist_random_initial` |

### Tables

| Output | Run |
|---------|-----|
| **Table 2** (Unique fitness, UF) | `sims/evo_sim` → `sims/evolved_J_evaluation/src/calc_all_J_x` → `sims/evolved_J_evaluation/analyze_all_J_fitness.py` → `sims/analysis/unique_fitness_histogram.py` |
| **Table 2** (Invalid count, IC) | `sims/evo_sim` → `sims/evolved_J_evaluation/src/calc_all_J_x` → `sims/evolved_J_evaluation/analyze_all_J_fitness.py` → `sims/analysis/invalid_count_histogram.py` |