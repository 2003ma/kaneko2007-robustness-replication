#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "all_J_trial_eval.h"

/* 実行例:
 * ./bin/all_J_trial_eval --input ../evo_sim/results/ver1/sigma_0.400 \
 */
/*---------------------------------------
  CLI パース
----------------------------------------*/
static void parse_args(int argc, char **argv, AllJTrialParams *p)
{
    p->input_dir = NULL;
    p->generation = 200;
    p->n_trials = 10000;
    p->t1 = 80;
    p->t2 = 90;
    p->beta = 7.0;
    p->dt = 0.05;
    p->k_boundary = 8;
    p->jemk = 1; /* デフォルトでj>=kを考慮 */
    p->seed = 12345;
    p->sigma = 0.0;  /* デフォルト値 */
    p->start_ind = 0; /* デフォルト0 */

    for (int i = 1; i < argc; i++)
    {
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
        else if (!strcmp(argv[i], "--sigma") && i + 1 < argc)
            p->sigma = atof(argv[++i]);
        else if (!strcmp(argv[i], "--start_ind") && i + 1 < argc)
            p->start_ind = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--jemk") && i + 1 < argc)
            p->jemk = atoi(argv[++i]);
    }

    if (!p->input_dir)
    {
        fprintf(stderr, "ERROR: --input DIR is required.\n");
        fprintf(stderr, "Usage: %s --input <dir> [options]\n", argv[0]);
        fprintf(stderr, "  --generation <N>    (default: 200)\n");
        fprintf(stderr, "  --trials <N>        (default: 10000)\n");
        fprintf(stderr, "  --t1 <N>            (default: 150)\n");
        fprintf(stderr, "  --t2 <N>            (default: 250)\n");
        fprintf(stderr, "  --beta <N>          (default: 7.0)\n");
        fprintf(stderr, "  --dt <N>            (default: 0.05)\n");
        fprintf(stderr, "  --k-boundary <N>    (default: 8)\n");
        fprintf(stderr, "  --seed <N>          (default: 12345)\n");
        fprintf(stderr, "  --sigma <N>         (default: 0.0)\n");
        fprintf(stderr, "  --start_ind <N>     (default: 0)\n");
        exit(1);
    }
}

/*---------------------------------------
  メイン関数
----------------------------------------*/
int main(int argc, char **argv)
{
    AllJTrialParams prm;
    parse_args(argc, argv, &prm);

    printf("[INFO] Starting all_J trial evaluation...\n");
    printf("[INFO] Input directory: %s\n", prm.input_dir);
    printf("[INFO] Generation: %d\n", prm.generation);
    printf("[INFO] Trials per individual: %d\n", prm.n_trials);
    printf("[INFO] t1=%d, t2=%d, dt=%.3f\n", prm.t1, prm.t2, prm.dt);

    return run_all_J_trial_eval(&prm);
}
