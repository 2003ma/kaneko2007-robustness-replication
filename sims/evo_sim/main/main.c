
#include "../src/kaneko_robustness.h"
#include "../src/save_all_J.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

// Default parameters
void set_default_params(Params *p)
{
    p->sigma = 0.2;
    p->N = 64; //M in the paper
    p->k = 8;
    p->pop = 300; //N in the paper
    p->gens = 201;
    p->beta = 7.0;
    p->dt = 0.005;
    p->relax = 80 / 0.005;
    p->meas = 10 / 0.005;
    p->L = 300;
    p->pedge = 0.5;
    p->elit = 0.25;
    p->mu = 1.00;
    p->seed = 123456789u;
    // p->load_J_path = "results/ver1/sigma_0.500/evo_sim_data/gen_120_all_J_sigma_0.500_dt0.005.csv"; // Use this when resuming from a midway point.
    p->load_J_path = NULL;
}

// Parse command-line arguments to set parameters (e.g., --sigma 0.01)
void parse_args(int argc, char **argv, Params *p)
{
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--sigma") == 0 && i + 1 < argc)
            p->sigma = atof(argv[++i]);
        else if (strcmp(argv[i], "--pop") == 0 && i + 1 < argc)
            p->pop = atoi(argv[++i]);
        else if (strcmp(argv[i], "--N") == 0 && i + 1 < argc)
            p->N = atoi(argv[++i]);
        else if (strcmp(argv[i], "--gens") == 0 && i + 1 < argc)
            p->gens = atoi(argv[++i]);
        else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc)
            p->seed = (uint32_t)atol(argv[++i]);
        else if (strcmp(argv[i], "--load_J") == 0 && i + 1 < argc)
            p->load_J_path = argv[++i]; // Specify CSV file path
        // ...add other parameters as needed...
    }
}

int main(int argc, char **argv)
{
    Params prm;
    set_default_params(&prm);
    parse_args(argc, argv, &prm);
    if (prm.seed == 0u)
    {
        prm.seed = (uint32_t)time(NULL);
        printf("[INFO] Seed not specified, using time-based seed: %u\n", prm.seed);
    }

    // --- Run the actual evolution simulation ---
    // (J save directory is auto-configured in kaneko_robustness.c)
    run_kaneko_robustness_simulation(&prm);
    printf("[INFO] Simulation, data saving, and J matrix saving every 20 generations completed\n");
    return 0;
}
