#!/usr/bin/env python3
"""
./plot_fitness_distribution_200gen.py \
        --csv-a ../../evo_sim/results/ver1/sigma_0.04/evo_sim_data/gen_200_all_J_sigma_0.04_dt0.005.csv \
        --csv-b ../../evo_sim/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
        --labels sigma=0.04 sigma=0.200 \
        --preset wide \
        --output fitness_wide.png
"""

from plot_fitness_distribution_200gen import main


if __name__ == "__main__":
    raise SystemExit(main(["--preset", "wide", "--output", "fitness_distribution_wide.png"] + __import__("sys").argv[1:]))
