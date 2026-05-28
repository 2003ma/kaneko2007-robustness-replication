#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Box-Muller 法による標準正規分布乱数生成
 * Xorshift32 を内部で利用
 */
double randn(unsigned int *seed);

/**
 * Euler-Maruyama 法による GRN ダイナミクス
 *   dx/dt = tanh(beta * Jx) - x + noise
 *
 * J:      N×N の結合行列（要素は -1,0,1）
 * N:      遺伝子数
 * beta:   シグモイドの鋭さ
 * sigma:  ノイズ強度
 * dt:     タイムステップ
 * T:      総時間
 * x0:     初期値 (長さ N)
 * seed:   乱数シード
 * k_boundary: ターゲット遺伝子数
 * jemk:    j>=k を考慮するかどうか (TRUE 1/FALSE 0)
 * X_out:  長さ (steps+1)*N の配列（t=0..T を格納）
 */
void em_integrate(const double *J, int N, double beta, double sigma,
                  double dt, double T, const double *x0,
                  int seed, int k_boundary, int jemk, double *X_out);

/**
 * 4次の Runge-Kutta 法による決定論ダイナミクス
 *   dx/dt = tanh(beta * Jx) - x
 * （ノイズなし）
 */
void rk4_integrate(const double *J, int N, double beta,
                   double dt, double T, const double *x0,
                   int k_boundary, int jemk, double *X_out);

/**
 * 制御領域固定版 4次 Runge-Kutta 法
 *   dx/dt = tanh(beta * Jx) - x
 * 制御領域 [N-reg_size, N) は常に x_i = 1.0 に固定
 */
void rk4_integrate_controlled(const double *J, int N, double beta,
                              double dt, double T, const double *x0,
                              int k_boundary, int reg_size, int jemk,
                              double *X_out);

/**
 * 外力付き Euler-Maruyama
 *   dx_i/dt = tanh(beta * (Jx)_i + f_i) - x_i + noise
 *
 * force_gene: 外力を加える遺伝子インデックス i
 * force:      その遺伝子に加える外力の大きさ（tanh の中に加算）
 */
void em_integrate_with_force(const double *J, int N, double beta, double sigma,
                             double dt, double T, const double *x0,
                             int k_boundary, int jemk, int force_gene, double force,
                             double *X_out, int seed);

/**
 * 外力付き 4次 Runge-Kutta
 *   dx_i/dt = tanh(beta * (Jx)_i + f_i) - x_i
 * （ノイズなし、クリップなし）
 */
void rk4_integrate_with_force(const double *J, int N, double beta,
                              double dt, double T, const double *x0,
                              int k_boundary, int jemk, int force_gene, double force,
                              double *X_out);

#ifdef __cplusplus
}
#endif

#endif // DYNAMICS_H
