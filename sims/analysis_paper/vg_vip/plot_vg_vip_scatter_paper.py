#!/usr/bin/env python3
"""
./plot_vg_vip_scatter_paper.py \
  --input-dir ../../evo_sim_v6/results/data_sigma_only/ver1 \
  --sigmas 0.001 0.006 0.01 0.02 0.04 0.06 0.08 0.1 0.15
"""

from plot_vg_vip_scatter import main


if __name__ == "__main__":
    raise SystemExit(main(["--preset", "paper", "--output", "vg_vip_scatter_paper.png"] + __import__("sys").argv[1:]))
