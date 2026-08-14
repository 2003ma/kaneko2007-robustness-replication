#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Generate standard normal random numbers using the Box-Muller method
 * Internally uses Xorshift32
 */
double randn(unsigned int *seed);

/**
 * GRN dynamics by the Euler-Maruyama method
 *   dx/dt = tanh(beta * Jx) - x + noise
 *
 * J:      N×N coupling matrix (elements are -1,0,1)
 * N:      Number of genes
 * beta:   Sigmoid steepness
 * sigma:  Noise strength
 * dt:     Time step
 * T:      Total time
 * x0:     Initial values (length N)
 * seed:   Random seed
 * k_boundary: Target gene count
 * jemk:    Whether to consider j>=k (TRUE 1/FALSE 0)
 * X_out:  Array of length (steps+1)*N (stores t=0..T)
 */
void em_integrate(const double *J, int N, double beta, double sigma,
                  double dt, double T, const double *x0,
                  int seed, int k_boundary, int jemk, double *X_out);

/**
 * Deterministic dynamics using the 4th-order Runge-Kutta method
 *   dx/dt = tanh(beta * Jx) - x
 * (no noise)
 */
void rk4_integrate(const double *J, int N, double beta,
                   double dt, double T, const double *x0,
                   int k_boundary, int jemk, double *X_out);

/**
 * 4th-order Runge-Kutta method with fixed control region
 *   dx/dt = tanh(beta * Jx) - x
 * Control region [N-reg_size, N) is always fixed to x_i = 1.0
 */
void rk4_integrate_controlled(const double *J, int N, double beta,
                              double dt, double T, const double *x0,
                              int k_boundary, int reg_size, int jemk,
                              double *X_out);

/**
 * Euler-Maruyama with external force
 *   dx_i/dt = tanh(beta * (Jx)_i + f_i) - x_i + noise
 *
 * force_gene: Gene index i to which the force is applied
 * force:      Magnitude of the force added to that gene (added inside tanh)
 */
void em_integrate_with_force(const double *J, int N, double beta, double sigma,
                             double dt, double T, const double *x0,
                             int k_boundary, int jemk, int force_gene, double force,
                             double *X_out, int seed);

/**
 * 4th-order Runge-Kutta with external force
 *   dx_i/dt = tanh(beta * (Jx)_i + f_i) - x_i
 * (no noise, no clipping)
 */
void rk4_integrate_with_force(const double *J, int N, double beta,
                              double dt, double T, const double *x0,
                              int k_boundary, int jemk, int force_gene, double force,
                              double *X_out);

#ifdef __cplusplus
}
#endif

#endif // DYNAMICS_H
