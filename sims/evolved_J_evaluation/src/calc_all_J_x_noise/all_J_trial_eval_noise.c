#include "all_J_trial_eval_noise.h"
#include "../core/dynamics.h"

#include <glob.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path, mode) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path, mode) mkdir(path, mode)
#endif

static void ensure_dir(const char *path)
{
    MKDIR(path, 0755);
}

static void fill_x0_minus_one(double *x0, int N)
{
    for (int i = 0; i < N; i++)
        x0[i] = -1.0;
}

// 計算時に入れるノイズを計算
static int build_sigma_grid(const AllJTrialNoiseParams *p, double **sigmas_out, int *count_out)
{
    if (p->sim_sigma_step <= 0.0 || p->sim_sigma_end < p->sim_sigma_start)
        return 0;

    int n_sigma = (int)floor((p->sim_sigma_end - p->sim_sigma_start) / p->sim_sigma_step + 1e-9) + 1;
    if (n_sigma <= 0)
        return 0;

    double *sigmas = (double *)malloc(sizeof(double) * n_sigma);
    if (!sigmas)
        return 0;

    for (int i = 0; i < n_sigma; i++)
        sigmas[i] = p->sim_sigma_start + (double)i * p->sim_sigma_step;

    *sigmas_out = sigmas;
    *count_out = n_sigma;
    return 1;
}

static double compute_fitness_single_noise(
    const double *J,
    int N,
    int k_boundary,
    double beta,
    double sigma,
    double dt,
    int t1,
    int t2,
    int seed,
    int trial_id,
    int jemk)
{
    int total_steps = (int)floor((double)t2 / dt);
    int start_step = (int)floor((double)t1 / dt);
    if (start_step > total_steps)
        return 0.0;

    double T_total = total_steps * dt;
    double *x0 = (double *)malloc(sizeof(double) * N);
    double *X_out = (double *)malloc((size_t)(total_steps + 1) * N * sizeof(double));
    if (!x0 || !X_out) {
        free(x0);
        free(X_out);
        return 0.0;
    }

    fill_x0_minus_one(x0, N);
    em_integrate(
        J,
        N,
        beta,
        sigma,
        dt,
        T_total,
        x0,
        seed + trial_id,
        k_boundary,
        jemk,
        X_out);

    double sum_on = 0.0;
    int n_meas = 0;
    for (int t = start_step; t <= total_steps; t++) {
        const double *x = &X_out[t * N];
        int on = 0;
        for (int a = 0; a < k_boundary; a++)
            on += (x[a] > 0.0);
        sum_on += on;
        n_meas++;
    }

    free(x0);
    free(X_out);

    if (n_meas <= 0)
        return 0.0;
    return sum_on / (double)n_meas - (double)k_boundary;
}

static int load_all_J_from_csv(const char *csv_path, double **J_out, int *N_out)
{
    FILE *fp = fopen(csv_path, "r");
    if (!fp) {
        char glob_pattern[512];
        strncpy(glob_pattern, csv_path, sizeof(glob_pattern) - 1);
        glob_pattern[sizeof(glob_pattern) - 1] = '\0';

        char *csv_ext = strstr(glob_pattern, ".csv");
        if (csv_ext)
            *csv_ext = '\0';
        strncat(glob_pattern, "*.csv", sizeof(glob_pattern) - strlen(glob_pattern) - 1);

        glob_t g;
        if (glob(glob_pattern, 0, NULL, &g) == 0 && g.gl_pathc > 0) {
            fp = fopen(g.gl_pathv[0], "r");
            if (fp)
                printf("[INFO] use %s\n", g.gl_pathv[0]);
        }
        globfree(&g);
    }

    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open %s\n", csv_path);
        return 0;
    }

    char line[65536];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return 0;
    }

    int skip_cols = strstr(line, "fitness") ? 3 : 1;
    int n_cols = 1;
    for (char *p = line; *p; p++)
        if (*p == ',')
            n_cols++;

    int n2 = n_cols - skip_cols;
    int N = (int)sqrt((double)n2);
    if (N * N != n2) {
        fprintf(stderr, "[ERROR] invalid J shape\n");
        fclose(fp);
        return 0;
    }

    int pop = 0;
    long pos = ftell(fp);
    while (fgets(line, sizeof(line), fp))
        pop++;
    fseek(fp, pos, SEEK_SET);

    double *J_all = (double *)malloc(sizeof(double) * pop * N * N);
    if (!J_all) {
        fclose(fp);
        return 0;
    }

    for (int ind = 0; ind < pop; ind++) {
        if (skip_cols == 3) {
            int id;
            double f, vip;
            if (fscanf(fp, "%d,%lf,%lf,", &id, &f, &vip) != 3) {
                free(J_all);
                fclose(fp);
                return 0;
            }
        } else {
            int id;
            if (fscanf(fp, "%d,", &id) != 1) {
                free(J_all);
                fclose(fp);
                return 0;
            }
        }

        double *J = &J_all[ind * N * N];
        for (int i = 0; i < N * N; i++) {
            if (fscanf(fp, "%lf,", &J[i]) != 1) {
                free(J_all);
                fclose(fp);
                return 0;
            }
        }
    }

    fclose(fp);
    *J_out = J_all;
    *N_out = N;
    return pop;
}

static int find_sims_root(const char *input_dir, char *out, size_t out_sz)
{
    const char *tag = "/sims/";
    const char *pos = strstr(input_dir, tag);
    if (!pos)
        return 0;
    size_t n = (size_t)(pos - input_dir) + strlen(tag) - 1;
    if (n >= out_sz)
        return 0;
    memcpy(out, input_dir, n);
    out[n] = '\0';
    return 1;
}

int run_all_J_trial_eval_noise(const AllJTrialNoiseParams *p)
{
    char input_dir_abs[PATH_MAX];
    const char *input_dir = p->input_dir;
    if (realpath(p->input_dir, input_dir_abs) != NULL)
        input_dir = input_dir_abs;

    char csv_path[512];
    snprintf(
        csv_path,
        sizeof(csv_path),
        "%s/evo_sim_data/gen_%d_all_J_sigma_%.3f_dt0.005.csv",
        input_dir,
        p->generation,
        p->evo_sigma);

    double *J_all = NULL;
    int N = 0;
    int pop_size = load_all_J_from_csv(csv_path, &J_all, &N);
    if (pop_size <= 0 || !J_all) {
        fprintf(stderr, "[ERROR] failed to load J from %s\n", csv_path);
        return 1;
    }
    if (p->k_boundary <= 0 || p->k_boundary > N) {
        fprintf(stderr, "[ERROR] invalid k_boundary\n");
        free(J_all);
        return 1;
    }

    double *sigmas = NULL;
    int n_sigma = 0;
    if (!build_sigma_grid(p, &sigmas, &n_sigma)) {
        fprintf(stderr, "[ERROR] invalid sim sigma range/step\n");
        free(J_all);
        return 1;
    }

    char sims_root[512];
    if (!find_sims_root(input_dir, sims_root, sizeof(sims_root))) {
        fprintf(stderr, "[ERROR] could not infer sims root from input_dir: %s\n", input_dir);
        free(J_all);
        free(sigmas);
        return 1;
    }

    char d0[512], d1[512], d2[512], d3[512];
    snprintf(d0, sizeof(d0), "%s/analysis_paper/under_noise_fitness", sims_root);
    ensure_dir(d0);
    snprintf(d1, sizeof(d1), "%s/evo_sigma_%.3f", d0, p->evo_sigma);
    ensure_dir(d1);
    snprintf(
        d2,
        sizeof(d2),
        "%s/sim_sigma_%.3f_%.3f_step_%.3f",
        d1,
        p->sim_sigma_start,
        p->sim_sigma_end,
        p->sim_sigma_step);
    ensure_dir(d2);
    snprintf(d3, sizeof(d3), "%s/gen_%d", d2, p->generation);
    ensure_dir(d3);

    int start_ind = p->start_ind;
    int end_ind = (p->end_ind < 0) ? (pop_size - 1) : p->end_ind;
    if (start_ind < 0 || start_ind >= pop_size || end_ind < 0 || end_ind >= pop_size || start_ind > end_ind) {
        fprintf(
            stderr,
            "[ERROR] invalid individual range: start_ind=%d end_ind=%d (pop_size=%d)\n",
            start_ind,
            end_ind,
            pop_size);
        free(J_all);
        free(sigmas);
        return 1;
    }

    int total_inds = end_ind - start_ind + 1;

    for (int ind = start_ind; ind <= end_ind; ind++) {
        const double *J = &J_all[ind * N * N];
        char out_csv[512];
        snprintf(
            out_csv,
            sizeof(out_csv),
            "%s/trial_fitness_gen=%d_ind=%d_evo_sigma=%.3f_sim_sigma=%.3f_%.3f_step_%.3f_seed=%d_trials=%d.csv",
            d3,
            p->generation,
            ind,
            p->evo_sigma,
            p->sim_sigma_start,
            p->sim_sigma_end,
            p->sim_sigma_step,
            p->seed,
            p->n_trials);

        FILE *fp = fopen(out_csv, "w");
        if (!fp) {
            fprintf(stderr, "[WARN] cannot create %s\n", out_csv);
            continue;
        }
        fprintf(fp, "trial_id");
        for (int si = 0; si < n_sigma; si++) {
            if (fabs(sigmas[si]) < 1e-12)
                fprintf(fp, ",0_fitness");
            else
                fprintf(fp, ",%.2f_fitness", sigmas[si]);
        }
        fprintf(fp, "\n");

        double *fitness_buf = (double *)malloc(sizeof(double) * p->n_trials * n_sigma);
        if (!fitness_buf) {
            fclose(fp);
            fprintf(stderr, "[WARN] malloc failed at ind=%d\n", ind);
            continue;
        }

        for (int si = 0; si < n_sigma; si++) {
            double sim_sigma = sigmas[si];
            #pragma omp parallel for schedule(static)
            for (int trial = 0; trial < p->n_trials; trial++) {
                fitness_buf[si * p->n_trials + trial] = compute_fitness_single_noise(
                    J,
                    N,
                    p->k_boundary,
                    p->beta,
                    sim_sigma,
                    p->dt,
                    p->t1,
                    p->t2,
                    p->seed,
                    trial,
                    p->jemk);
            }
        }

        for (int trial = 0; trial < p->n_trials; trial++) {
            fprintf(fp, "%d", trial);
            for (int si = 0; si < n_sigma; si++)
                fprintf(fp, ",%.8f", fitness_buf[si * p->n_trials + trial]);
            fprintf(fp, "\n");
        }

        free(fitness_buf);
        fclose(fp);

        if ((ind - start_ind + 1) % 20 == 0 || ind == end_ind)
            printf("[INFO] done ind %d/%d (global ind=%d)\n", ind - start_ind + 1, total_inds, ind);
    }

    free(J_all);
    free(sigmas);
    printf("[INFO] output root: %s\n", d3);
    return 0;
}
