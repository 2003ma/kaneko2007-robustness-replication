#!/usr/bin/env python3
"""
./plot_vg_vip_scatter_wide.py \
  --input-dir ../../evo_sim/results/data_sigma_only/ver1 \
  --sigmas 0.01 0.04 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
"""

from plot_vg_vip_scatter import main


if __name__ == "__main__":
    raise SystemExit(main(["--preset", "wide", "--output", "vg_vip_scatter_wide.png"] + __import__("sys").argv[1:]))
