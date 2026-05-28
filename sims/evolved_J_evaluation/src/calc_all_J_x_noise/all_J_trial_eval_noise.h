#ifndef ALL_J_TRIAL_EVAL_NOISE_H
#define ALL_J_TRIAL_EVAL_NOISE_H

typedef struct {
    const char *input_dir;
    int generation;
    int n_trials;
    int t1;
    int t2;
    double beta;
    double dt;
    int k_boundary;
    int seed;
    double evo_sigma;
    double sim_sigma_start;
    double sim_sigma_end;
    double sim_sigma_step;
    int start_ind;
    int end_ind;
    int jemk;
} AllJTrialNoiseParams;

int run_all_J_trial_eval_noise(const AllJTrialNoiseParams *p);

#endif
