#ifndef ALL_J_TRIAL_EVAL_H
#define ALL_J_TRIAL_EVAL_H

/* Trial evaluation parameters for all J */
typedef struct {
    const char *input_dir;   /* Directory containing gen_XXX_all_J.csv */
    int generation;          /* Generation number */
    int n_trials;            /* Number of trials per J */
    int t1;                  /* Measurement start time */
    int t2;                  /* Measurement end time */
    double beta;             /* Dynamics parameter */
    double dt;               /* Time step */
    int k_boundary;          /* Hierarchy boundary */
    int seed;                /* RNG seed */
    double sigma;            /* Mutation rate sigma */
    int start_ind;           /* Starting individual index (default 0) */
    int jemk;                /* Whether to consider j>=k in dynamics */
} AllJTrialParams;

/* Run function */
int run_all_J_trial_eval(const AllJTrialParams *p);

#endif /* ALL_J_TRIAL_EVAL_H */
