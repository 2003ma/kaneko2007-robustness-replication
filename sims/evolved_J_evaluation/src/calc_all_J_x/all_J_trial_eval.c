#include "all_J_trial_eval.h"
#include "../core/dynamics.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <glob.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*-------------------------
  splitmix RNG
-------------------------*/
static inline uint32_t splitmix32(uint32_t x) {
    x += 0x9E3779B9u;
    x = (x ^ (x >> 16)) * 0x85EBCA6Bu;
    x = (x ^ (x >> 13)) * 0xC2B2AE35u;
    return x ^ (x >> 16);
}

static inline double splitmix_uniform(uint32_t *state) {
    uint32_t r = splitmix32(*state);
    *state = r;
    return (r & 0xFFFFFFu) / (double)(1u << 24);
}

/*-------------------------
    x0 generation
-------------------------*/
static void generate_x0(double *x0, int N, uint32_t seed, int trial_id)
{
    uint32_t st = seed ^ (0xA5A5A5A5u * (uint32_t)trial_id);
    for (int i = 0; i < N; i++) {
        double u = splitmix_uniform(&st);
        x0[i] = 2*u - 1;   /* -1..1 */
    }
}

/*-------------------------
    Fitness computation (1 trial)
    x_mean: time-averaged gene expression levels for genes 0..N-1 (between t1 and t2)
    Run the simulation by calling rk4_integrate
-------------------------*/
static double compute_fitness_single(
    const double *J, int N,
    int k_boundary,
    double beta, double dt,
    int t1, int t2,
    const double *x0,
    double *x_mean       /* Buffer of length N (skip if NULL) */
, int jemk // 0: input from all j, 1: input only from j >= k_boundary
)
{
    int total_steps = (int)floor((double)t2 / dt);
    int start_step  = (int)floor((double)t1 / dt);

    int measure_steps = total_steps - start_step;
    if (measure_steps < 1) return 0.0;

    double T_total = total_steps * dt;

    double *X_out =
        (double *)malloc((size_t)(total_steps + 1) * N * sizeof(double));
    if (!X_out) return 0.0;

    rk4_integrate(J, N, beta, dt, T_total, x0, k_boundary, jemk, X_out);

    double sum_on = 0.0;

    /* Temporary buffer for mean expression levels */
    double *x_sum = NULL;
    if (x_mean != NULL) {
        x_sum = (double *)calloc((size_t)N, sizeof(double));
        if (!x_sum) {
            free(X_out);
            return 0.0;
        }
    }

    for (int t = start_step; t <= total_steps; t++) {
        const double *x = &X_out[t * N];

        int on = 0;
        for (int a = 0; a < k_boundary; a++)
            on += (x[a] > 0.0);

        sum_on += on;

        /* Accumulate values of each gene (for mean expression levels) */
        if (x_sum != NULL) {
            for (int a = 0; a < N; a++) {
                x_sum[a] += x[a];
            }
        }
    }

    double mean_on = sum_on / (double)(measure_steps + 1);

    /* Write mean expression levels to x_mean */
    if (x_sum != NULL) {
        double inv_steps = 1.0 / (double)(measure_steps + 1);  // Divide by the number of steps to obtain the per-step average
        for (int a = 0; a < N; a++) {
            x_mean[a] = x_sum[a] * inv_steps;
        }
        free(x_sum);
    }

    free(X_out);

    return mean_on - (double)k_boundary;
}

/*-------------------------
    Load all J matrices from gen_XXX_all_J*.csv

    returns: number of individuals (population_size)
                     J_out: malloc'ed J matrix array (population_size × N × N)
                     N_out: number of genes
-------------------------*/
static int load_all_J_from_csv(const char *csv_path, double **J_out, int *N_out)
{
    /* First try the specified path as-is */
    FILE *fp = fopen(csv_path, "r");
    char actual_path[512];
    strncpy(actual_path, csv_path, sizeof(actual_path) - 1);
    
    /* If that fails, search for the file using a wildcard */
    if (!fp) {
        /* Remove .csv and build a wildcard pattern */
        char glob_pattern[512];
        strncpy(glob_pattern, csv_path, sizeof(glob_pattern) - 1);
        
        /* Remove trailing .csv */
        char *csv_ext = strstr(glob_pattern, ".csv");
        if (csv_ext) {
            *csv_ext = '\0';
        }
        
        /* Append *.csv */
        strncat(glob_pattern, "*.csv", sizeof(glob_pattern) - strlen(glob_pattern) - 1);
        
        glob_t globbuf;
        if (glob(glob_pattern, 0, NULL, &globbuf) == 0 && globbuf.gl_pathc > 0) {
            /* Use the first file found */
            fp = fopen(globbuf.gl_pathv[0], "r");
            if (fp) {
                printf("[INFO] Found and using: %s\n", globbuf.gl_pathv[0]);
                strncpy(actual_path, globbuf.gl_pathv[0], sizeof(actual_path) - 1);
            }
            globfree(&globbuf);
        }
    }
    
    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open %s or matching pattern\n", csv_path);
        return 0;
    }

    /* Read header */
    char header[65536];
    if (!fgets(header, sizeof(header), fp)) {
        fclose(fp);
        return 0;
    }

    /* Detect format: check whether "fitness" is included */
    int skip_cols;
    if (strstr(header, "fitness")) {
        /* evo_sim_v2 format: individual_id,fitness,v_ip,J_0_0,... */
        skip_cols = 3;
    } else {
        /* evo_sim format: individual_id,J0,J1,... */
        skip_cols = 1;
    }

    /* Count columns from the header and determine N */
    int n_cols = 1;
    for (char *p = header; *p; p++) {
        if (*p == ',') n_cols++;
    }
    int N_squared = n_cols - skip_cols;
    int N = (int)sqrt((double)N_squared);
    if (N * N != N_squared) {
        fprintf(stderr, "[ERROR] J columns not square: %d (total cols: %d, skip_cols: %d)\n", 
                N_squared, n_cols, skip_cols);
        fclose(fp);
        return 0;
    }
    *N_out = N;

    /* Count rows (individuals) */
    int pop_size = 0;
    long pos = ftell(fp);
    while (fgets(header, sizeof(header), fp)) {
        pop_size++;
    }
    fseek(fp, pos, SEEK_SET);

    /* Allocate memory */
    double *J_all = (double *)malloc(sizeof(double) * pop_size * N * N);
    if (!J_all) {
        fprintf(stderr, "[ERROR] malloc failed for J_all\n");
        fclose(fp);
        return 0;
    }

    /* Load data */
    for (int ind = 0; ind < pop_size; ind++) {
        /* Skip leading columns (such as individual_id) */
        if (skip_cols == 3) {
            /* evo_sim_v2 format: individual_id, fitness, v_ip */
            int id;
            double fitness, v_ip;
            if (fscanf(fp, "%d,%lf,%lf,", &id, &fitness, &v_ip) != 3) {
                fprintf(stderr, "[ERROR] fscanf failed at individual %d header (v2 format)\n", ind);
                free(J_all);
                fclose(fp);
                return 0;
            }
        } else {
            /* evo_sim format: individual_id */
            int id;
            if (fscanf(fp, "%d,", &id) != 1) {
                fprintf(stderr, "[ERROR] fscanf failed at individual %d header (v1 format)\n", ind);
                free(J_all);
                fclose(fp);
                return 0;
            }
        }
        
        /* Read J matrix elements */
        double *J = &J_all[ind * N * N];
        for (int i = 0; i < N * N; i++) {
            if (fscanf(fp, "%lf,", &J[i]) != 1) {
                fprintf(stderr, "[ERROR] fscanf failed at J element %d\n", i);
                free(J_all);
                fclose(fp);
                return 0;
            }
        }
    }

    fclose(fp);
    *J_out = J_all;
    return pop_size;
}

/*-------------------------
    Directory creation
-------------------------*/
#ifdef _WIN32
#include <direct.h>
#define MKDIR(path, mode) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path, mode) mkdir(path, mode)
#endif

static void ensure_dir(const char *path) {
    MKDIR(path, 0755);
}

/*-------------------------
    Main logic
-------------------------*/
int run_all_J_trial_eval(const AllJTrialParams *p)
{
    /* Extract the sigma value from the directory name */
    double dir_sigma = p->sigma;
    const char *sigma_str = strstr(p->input_dir, "sigma_");
    if (sigma_str) {
        sscanf(sigma_str, "sigma_%lf", &dir_sigma);
        printf("[INFO] Extracted sigma=%.3f from directory name\n", dir_sigma);
    }

    /* Build the path to gen_XXX_all_J.csv */
    char csv_path[512];
    snprintf(csv_path, sizeof(csv_path),
             "%s/evo_sim_data/gen_%d_all_J_sigma_%.3f_dt0.005.csv",
             p->input_dir, p->generation, dir_sigma);

    /* Load all J matrices */
    double *J_all = NULL;
    int N = 0;
    int pop_size = load_all_J_from_csv(csv_path, &J_all, &N);
    if (pop_size <= 0 || !J_all) {
        fprintf(stderr, "[ERROR] failed to load %s\n", csv_path);
        return 1;
    }

    printf("[INFO] Loaded %d individuals, N=%d from %s\n", pop_size, N, csv_path);

    if (p->k_boundary <= 0 || p->k_boundary > N) {
        fprintf(stderr, "[ERROR] invalid k_boundary\n");
        free(J_all);
        return 1;
    }

    /* Output directory structure */
    char dir2[512];
    snprintf(dir2, sizeof(dir2), "%s/trials%d", p->input_dir, p->n_trials);
    ensure_dir(dir2);

    char dir3[512];
    snprintf(dir3, sizeof(dir3),
             "%s/t1=%d_t2=%d_trials%d_dt%.3f_allJ",
             dir2, p->t1, p->t2, p->n_trials, p->dt);
    ensure_dir(dir3);

    /* Easy-to-read folder name: gen_XXX_all_individuals */
    char dir4[512];
    snprintf(dir4, sizeof(dir4), "%s/gen_%d_all_individuals", dir3, p->generation);
    ensure_dir(dir4);

    char dir5[512];
    snprintf(dir5, sizeof(dir5), "%s/allx", dir4);
    ensure_dir(dir5);

    /* Create a CSV file for each individual and perform trial evaluation */
    for (int ind = p->start_ind; ind < pop_size; ind++) {
        const double *J = &J_all[ind * N * N];

        /* CSV filename (includes sigma) */
        char csv[512];
        snprintf(csv, sizeof(csv),
                 "%s/trial_fitness_gen=%d_ind=%d_sigma=%.3f_seed=%d_trials=%d_dt%.3f_allx.csv",
                 dir5, p->generation, ind, dir_sigma, p->seed, p->n_trials, p->dt);

        FILE *fp = fopen(csv, "w");
        if (!fp) {
            fprintf(stderr, "[ERROR] cannot create %s\n", csv);
            continue;
        }

        /* Header */
        fprintf(fp, "trial_id,fitness");
        for (int i = 0; i < N; i++)
            fprintf(fp, ",x%d", i);
        fprintf(fp, "\n");
        fflush(fp);

        /* Store trial results in temporary buffers */
        double *fitness_buf = malloc(sizeof(double) * p->n_trials);
        double *x0_buf = malloc(sizeof(double) * p->n_trials * N);
        double *xmean_buf = malloc(sizeof(double) * p->n_trials * N);
        if (!fitness_buf || !x0_buf || !xmean_buf) {
            fprintf(stderr, "[ERROR] malloc failed for individual %d\n", ind);
            fclose(fp);
            if (fitness_buf) free(fitness_buf);
            if (x0_buf) free(x0_buf);
            if (xmean_buf) free(xmean_buf);
            continue;
        }

        #pragma omp parallel for schedule(static)
        for (int trial = 0; trial < p->n_trials; trial++) {
            double *x0 = &x0_buf[trial * N];
            double *x_mean = &xmean_buf[trial * N];
            generate_x0(x0, N, (uint32_t)p->seed, trial);
            double fit = compute_fitness_single(
                J, N, p->k_boundary,
                p->beta, p->dt,
                p->t1, p->t2,
                x0,
                x_mean, p->jemk
            );
            fitness_buf[trial] = fit;
        }

        /* Write CSV in ascending trial_id order */
        for (int trial = 0; trial < p->n_trials; trial++) {
            fprintf(fp, "%d,%.8f", trial, fitness_buf[trial]);

            /* Mean expression values x0..x(N-1) */
            for (int i = 0; i < N; i++)
                fprintf(fp, ",%.8f", xmean_buf[trial * N + i]);

            fprintf(fp, "\n");
        }
        
        /* Save initial values separately only for start_ind (to confirm all individuals share the same initial values) */
        if (ind == p->start_ind) {
            char csv_x0[512];
            snprintf(csv_x0, sizeof(csv_x0),
                     "%s/trial_fitness_gen=%d_sigma=%.3f_seed=%d_x0.csv",
                     dir4, p->generation, dir_sigma, p->seed);
            FILE *fp_x0 = fopen(csv_x0, "w");
            if (fp_x0) {
                fprintf(fp_x0, "trial_id");
                for (int i = 0; i < N; i++)
                    fprintf(fp_x0, ",x0_%d", i);
                fprintf(fp_x0, "\n");
                for (int trial = 0; trial < p->n_trials; trial++) {
                    fprintf(fp_x0, "%d", trial);
                    for (int i = 0; i < N; i++)
                        fprintf(fp_x0, ",%.8f", x0_buf[trial * N + i]);
                    fprintf(fp_x0, "\n");
                }
                fclose(fp_x0);
            }
        }

        free(fitness_buf);
        free(x0_buf);
        free(xmean_buf);
        fclose(fp);

        printf("[INFO] Completed individual %d/%d\n", ind + 1, pop_size);
    }

    free(J_all);
    printf("[INFO] All individuals processed successfully.\n");
    return 0;
}
