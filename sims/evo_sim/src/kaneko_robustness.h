#ifndef KANEKO_ROBUSTNESS_H
#define KANEKO_ROBUSTNESS_H
#include <stdint.h>

typedef struct
{
    int N;
    int k;
    int pop;
    int gens;
    double beta;
    double sigma;
    double dt;
    int relax;
    int meas;
    int L;
    double pedge;
    double elit;
    double mu;
    uint32_t seed;
    const char *load_J_path;
} Params;

void run_kaneko_robustness_simulation(const Params *prm);
void save_all_J_final(const double *J_list, int pop_size, int N, const char *out_path);
void set_J_save_dir(const char *dir);

#endif
