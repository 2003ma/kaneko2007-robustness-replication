#!/usr/bin/env python3
"""
./plot_fitness_distribution_200gen.py \
        --csv-a ../../evo_sim_v6/results/ver1/sigma_0.005/evo_sim_data/gen_200_all_J_sigma_0.005_dt0.005.csv \
        --csv-b ../../evo_sim_v6/results/ver1/sigma_0.200/evo_sim_data/gen_200_all_J_sigma_0.200_dt0.005.csv \
        --labels sigma=0.005 sigma=0.200 \
        --preset zoom \
        --output fitness_zoom.png
"""

from plot_fitness_distribution_200gen import main


if __name__ == "__main__":
    raise SystemExit(main(["--preset", "zoom", "--output", "fitness_distribution_zoom.png"] + __import__("sys").argv[1:]))
