#ifndef J_LOADER_H
#define J_LOADER_H

/* 
 * dir: best_J / worst_J があるディレクトリ（sigma_0.08 みたいな）
 * file_type: "best" or "worst"
 * generation: -1 → 最終行、それ以外 → その gen の行
 * N_out: J のサイズ N（N×N）
 *
 * 戻り値: malloc された double *J (N*N), 失敗時 NULL
 */
double *load_J_from_dir(const char *dir,
                        const char *file_type,
                        int generation,
                        int *N_out);

#endif /* J_LOADER_H */
