#ifndef J_LOADER_H
#define J_LOADER_H

/* 
 * dir: directory containing best_J / worst_J files (e.g. sigma_0.08)
 * file_type: "best" or "worst"
 * generation: -1 -> final row, otherwise -> the specified generation row
 * N_out: size N of J (N×N)
 *
 * Return value: malloc'ed double *J (N*N), or NULL on failure
 */
double *load_J_from_dir(const char *dir,
                        const char *file_type,
                        int generation,
                        int *N_out);

#endif /* J_LOADER_H */
