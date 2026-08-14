#include "dynamics.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*------------------------------------------------------------
 * Internal structure: precomputed indices for J ∈ {−1,0,1}
 *----------------------------------------------------------*/

typedef struct {
    int N;
    int *pos_count;   // Number of +1 entries in each row
    int *neg_count;   // Number of -1 entries in each row
    int **pos_idx;    // pos_idx[i][k] = j (J_ij = +1)
    int **neg_idx;    // neg_idx[i][k] = j (J_ij = -1)
} JIdx;

/* Create an index structure from J (must be freed) */
static JIdx *JIdx_create(const double *J, int N)
{
    JIdx *idx = (JIdx *)malloc(sizeof(JIdx));
    if (!idx)
        return NULL;

    idx->N = N;
    idx->pos_count = (int *)calloc(N, sizeof(int));
    idx->neg_count = (int *)calloc(N, sizeof(int));
    if (!idx->pos_count || !idx->neg_count) {
        free(idx->pos_count);
        free(idx->neg_count);
        free(idx);
        return NULL;
    }

    /*
    * Ideally J contains {-1, 0, +1}.
    * Using 0.5 / -0.5 thresholds yields three categories:
    *   Jij >  0.5 -> treat as +1
    *   Jij < -0.5 -> treat as -1
    *   otherwise   -> treat as 0
    * This makes +1/-1 classification explicit and robust to rounding errors.
     */
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double Jij = J[i * N + j];
            if (Jij > 0.5)         // treat as +1
                idx->pos_count[i]++;
            else if (Jij < -0.5)   // treat as -1
                idx->neg_count[i]++;
        }
    }

    idx->pos_idx = (int **)malloc(sizeof(int *) * N);
    idx->neg_idx = (int **)malloc(sizeof(int *) * N);
    if (!idx->pos_idx || !idx->neg_idx) {
        free(idx->pos_idx);
        free(idx->neg_idx);
        free(idx->pos_count);
        free(idx->neg_count);
        free(idx);
        return NULL;
    }

    for (int i = 0; i < N; ++i) {
        int pc = idx->pos_count[i];
        int nc = idx->neg_count[i];
        idx->pos_idx[i] = pc > 0 ? (int *)malloc(sizeof(int) * pc) : NULL;
        idx->neg_idx[i] = nc > 0 ? (int *)malloc(sizeof(int) * nc) : NULL;
    }

    /* Fill the indices */
    int *pc = (int *)calloc(N, sizeof(int));
    int *nc = (int *)calloc(N, sizeof(int));
    if (!pc || !nc) {
        free(pc);
        free(nc);
        for (int i = 0; i < N; ++i) {
            free(idx->pos_idx[i]);
            free(idx->neg_idx[i]);
        }
        free(idx->pos_idx);
        free(idx->neg_idx);
        free(idx->pos_count);
        free(idx->neg_count);
        free(idx);
        return NULL;
    }

    /* For each row i, store +1 columns in pos_idx and -1 columns in neg_idx */
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double Jij = J[i * N + j];
            if (Jij > 0.5) {
                idx->pos_idx[i][pc[i]++] = j;
            } else if (Jij < -0.5) {
                idx->neg_idx[i][nc[i]++] = j;
            }
        }
    }

    free(pc);
    free(nc);

    return idx;
}

/* Free JIdx */
static void JIdx_free(JIdx *idx)
{
    if (!idx)
        return;

    if (idx->pos_idx) {
        for (int i = 0; i < idx->N; ++i)
            free(idx->pos_idx[i]);
        free(idx->pos_idx);
    }
    if (idx->neg_idx) {
        for (int i = 0; i < idx->N; ++i)
            free(idx->neg_idx[i]);
        free(idx->neg_idx);
    }
    free(idx->pos_count);
    free(idx->neg_count);
    free(idx);
}

/* Compute h = Jx quickly using the index structure (OpenMP) */
static void compute_h_fast(const JIdx *idx,
                           const double *x,
                           int k_boundary,
                           double *h, 
                           int jemk)
{
    int N = idx->N;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        double hi = 0.0;
        int pc = idx->pos_count[i];
        int nc = idx->neg_count[i];

        if (jemk == 0) {
            /* Input from all j (in practice, only jemk=1 is used) */
            for (int k = 0; k < pc; ++k) {
                int j = idx->pos_idx[i][k];
                hi += x[j];
            }
            for (int k = 0; k < nc; ++k) {
                int j = idx->neg_idx[i][k];
                hi -= x[j];
            }
        } else {
            /* Only j >= k_boundary */
            for (int k = 0; k < pc; ++k) {
                int j = idx->pos_idx[i][k];
                if (j >= k_boundary)
                    hi += x[j];
            }
            for (int k = 0; k < nc; ++k) {
                int j = idx->neg_idx[i][k];
                if (j >= k_boundary)
                    hi -= x[j];
            }
        }

        h[i] = hi;
    }
}

/*------------------------------------------------------------
 * Random numbers: Box-Muller + Xorshift32
 *----------------------------------------------------------*/

double randn(unsigned int *state)
{
    // Xorshift32
    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    double u1 = ((*state >> 1) + 1) / (double)0x7FFFFFFF;

    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    double u2 = ((*state >> 1) + 1) / (double)0x7FFFFFFF;

    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/*------------------------------------------------------------
 * Euler-Maruyama method (fast version)
 *----------------------------------------------------------*/

void em_integrate(const double *J, int N, double beta, double sigma,
                  double dt, double T, const double *x0,
                  int seed, int k_boundary, int jemk, double *X_out)
{
    int steps = (int)(T / dt);

    JIdx *idx = JIdx_create(J, N);
    if (!idx) {
        fprintf(stderr, "em_integrate: failed to create JIdx\n");
        return;
    }

    double *x    = (double *)malloc(sizeof(double) * N);
    double *h    = (double *)malloc(sizeof(double) * N);
    double *dxdt = (double *)malloc(sizeof(double) * N);

    if (!x || !h || !dxdt) {
        fprintf(stderr, "em_integrate: malloc failed\n");
        free(x);
        free(h);
        free(dxdt);
        JIdx_free(idx);
        return;
    }

    memcpy(x, x0, sizeof(double) * N);

    unsigned int rng = (unsigned int)seed;

    for (int t = 0; t <= steps; ++t) {
        /* Output */
        memcpy(&X_out[t * N], x, sizeof(double) * N);

        /* Compute h = Jx (OpenMP) */
        compute_h_fast(idx, x, k_boundary, h, jemk);

        /* Compute dx/dt */
        for (int i = 0; i < N; ++i) {
            dxdt[i] = tanh(beta * h[i]) - x[i];
        }

        /* Advance one step + noise & clipping */
        for (int i = 0; i < N; ++i) {
            double noise = sigma * sqrt(dt) * randn(&rng);
            x[i] += dt * dxdt[i] + noise;

            // This clipping does not materially change the result
            if (x[i] < -3.0)
                x[i] = -3.0;
            if (x[i] > 3.0)
                x[i] = 3.0;
        }
    }

    free(x);
    free(h);
    free(dxdt);
    JIdx_free(idx);
}

/*------------------------------------------------------------
 *  4th-order Runge-Kutta method (fast version)
 *----------------------------------------------------------*/
// jemk=0: input from all j
// jemk=1: input only from j >= k_boundary
void rk4_integrate(const double *J, int N, double beta,
                   double dt, double T, const double *x0,
                   int k_boundary, int jemk, double *X_out)
{
    int steps = (int)(T / dt);

    JIdx *idx = JIdx_create(J, N);
    if (!idx) {
        fprintf(stderr, "rk4_integrate: failed to create JIdx\n");
        return;
    }

    double *x    = (double *)malloc(sizeof(double) * N);
    double *k1   = (double *)malloc(sizeof(double) * N);
    double *k2   = (double *)malloc(sizeof(double) * N);
    double *k3   = (double *)malloc(sizeof(double) * N);
    double *k4   = (double *)malloc(sizeof(double) * N);
    double *h    = (double *)malloc(sizeof(double) * N);
    double *xtmp = (double *)malloc(sizeof(double) * N);

    if (!x || !k1 || !k2 || !k3 || !k4 || !h || !xtmp) {
        fprintf(stderr, "rk4_integrate: malloc failed\n");
        free(x); free(k1); free(k2); free(k3); free(k4);
        free(h); free(xtmp);
        JIdx_free(idx);
        return;
    }

    memcpy(x, x0, sizeof(double) * N);

    for (int t = 0; t <= steps; ++t) {
        /* Output */
        memcpy(&X_out[t * N], x, sizeof(double) * N);

        /* k1 = dt * f(x) */
        compute_h_fast(idx, x, k_boundary, h, jemk);
        for (int i = 0; i < N; ++i)
            k1[i] = dt * (tanh(beta * h[i]) - x[i]);

        /* k2 = dt * f(x + k1/2) */
        for (int i = 0; i < N; ++i)
            xtmp[i] = x[i] + 0.5 * k1[i];
        compute_h_fast(idx, xtmp, k_boundary, h, jemk);
        for (int i = 0; i < N; ++i)
            k2[i] = dt * (tanh(beta * h[i]) - xtmp[i]);

        /* k3 = dt * f(x + k2/2) */
        for (int i = 0; i < N; ++i)
            xtmp[i] = x[i] + 0.5 * k2[i];
        compute_h_fast(idx, xtmp, k_boundary, h, jemk);
        for (int i = 0; i < N; ++i)
            k3[i] = dt * (tanh(beta * h[i]) - xtmp[i]);

        /* k4 = dt * f(x + k3) */
        for (int i = 0; i < N; ++i)
            xtmp[i] = x[i] + k3[i];
        compute_h_fast(idx, xtmp, k_boundary, h, jemk);
        for (int i = 0; i < N; ++i)
            k4[i] = dt * (tanh(beta * h[i]) - xtmp[i]);

        /* x <- x + (k1 + 2k2 + 2k3 + k4)/6 */
        for (int i = 0; i < N; ++i)
            x[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }

    free(x);
    free(k1); free(k2); free(k3); free(k4);
    free(h);
    free(xtmp);
    JIdx_free(idx);
}