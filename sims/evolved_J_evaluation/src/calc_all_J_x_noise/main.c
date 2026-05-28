#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "all_J_trial_eval_noise.h"

/* 実行例:
 *  ./bin/all_J_trial_eval_noise \
 *   --input ../../../evo_sim/results/ver1/sigma_0.005 --evo-sigma 0.005
 */
static void parse_args(int argc, char **argv, AllJTrialNoiseParams *p)
{
    p->input_dir = NULL;
    p->generation = 200;
    p->n_trials = 10000;
    p->t1 = 80;
    p->t2 = 90;
    p->beta = 7.0;
    p->dt = 0.05;
    p->k_boundary = 8;
    p->jemk = 1;
    p->seed = 12345;
    p->evo_sigma = 0.300;
    p->sim_sigma_start = 0.0;
    p->sim_sigma_end = 0.4;
    p->sim_sigma_step = 0.01;
    p->start_ind = 0;
    p->end_ind = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--input") && i + 1 < argc)
            p->input_dir = argv[++i];
        else if (!strcmp(argv[i], "--generation") && i + 1 < argc)
            p->generation = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--trials") && i + 1 < argc)
            p->n_trials = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--t1") && i + 1 < argc)
            p->t1 = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--t2") && i + 1 < argc)
            p->t2 = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--beta") && i + 1 < argc)
            p->beta = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt") && i + 1 < argc)
            p->dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "--k-boundary") && i + 1 < argc)
            p->k_boundary = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc)
            p->seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--evo-sigma") && i + 1 < argc)
            p->evo_sigma = atof(argv[++i]);
        else if (!strcmp(argv[i], "--sim-sigma-start") && i + 1 < argc)
            p->sim_sigma_start = atof(argv[++i]);
        else if (!strcmp(argv[i], "--sim-sigma-end") && i + 1 < argc)
            p->sim_sigma_end = atof(argv[++i]);
        else if (!strcmp(argv[i], "--sim-sigma-step") && i + 1 < argc)
            p->sim_sigma_step = atof(argv[++i]);
        else if (!strcmp(argv[i], "--start_ind") && i + 1 < argc)
            p->start_ind = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--end_ind") && i + 1 < argc)
            p->end_ind = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--jemk") && i + 1 < argc)
            p->jemk = atoi(argv[++i]);
    }

    if (!p->input_dir || p->evo_sigma < 0.0) {
        fprintf(stderr, "ERROR: --input DIR and --evo-sigma are required.\n");
        fprintf(stderr, "Usage: %s --input <dir> --evo-sigma <sigma> [options]\n", argv[0]);
        fprintf(stderr, "  --sim-sigma-start <N>  dynamics noise sigma start (default: 0.0)\n");
        fprintf(stderr, "  --sim-sigma-end <N>    dynamics noise sigma end (default: 0.2)\n");
        fprintf(stderr, "  --sim-sigma-step <N>   dynamics noise sigma step (default: 0.01)\n");
        fprintf(stderr, "  --generation <N>    (default: 200)\n");
        fprintf(stderr, "  --trials <N>        (default: 10000)\n");
        fprintf(stderr, "  --t1 <N>            (default: 80)\n");
        fprintf(stderr, "  --t2 <N>            (default: 90)\n");
        fprintf(stderr, "  --beta <N>          (default: 7.0)\n");
        fprintf(stderr, "  --dt <N>            (default: 0.05)\n");
        fprintf(stderr, "  --k-boundary <N>    (default: 8)\n");
        fprintf(stderr, "  --seed <N>          (default: 12345)\n");
        fprintf(stderr, "  --start_ind <N>     (default: 0)\n");
        fprintf(stderr, "  --end_ind <N>       inclusive end index (default: -1 = last individual)\n");
        exit(1);
    }
}

int main(int argc, char **argv)
{
    AllJTrialNoiseParams prm;
    parse_args(argc, argv, &prm);

    printf("[INFO] all_J_trial_eval_noise start\n");
    printf("[INFO] input_dir=%s\n", prm.input_dir);
    printf(
        "[INFO] evo_sigma=%.3f sim_sigma=[%.3f:%.3f:%.3f]\n",
        prm.evo_sigma,
        prm.sim_sigma_start,
        prm.sim_sigma_end,
        prm.sim_sigma_step);
    printf("[INFO] generation=%d trials=%d\n", prm.generation, prm.n_trials);

    return run_all_J_trial_eval_noise(&prm);
}
