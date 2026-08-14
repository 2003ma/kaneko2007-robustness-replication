#include "j_loader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

/*-------------------------------------
    Search for best_J*.csv / worst_J*.csv in a directory
    If multiple files exist, return them as an array (sorted by descending number)
--------------------------------------*/
static char **find_all_csv(const char *dir, const char *file_type, int *count_out)
{
    DIR *dp = opendir(dir);
    if (!dp) return NULL;

    static char *paths[10];  /* Up to 10 files */
    static int nums[10];     /* Number extracted from each file */
    int count = 0;

    struct dirent *e;

    /* Example: search for best_J,sigma=...csv / worst_J,...csv */
    while ((e = readdir(dp)) && count < 10) {
        if (strstr(e->d_name, file_type) && strstr(e->d_name, "_J")) {
            // Extract the number in (N) format
            int num = -1;
            char *paren = strstr(e->d_name, "(");
            if (paren) {
                num = atoi(paren + 1);
            }
            
            paths[count] = (char *)malloc(1024);
            snprintf(paths[count], 1024, "%s/%s", dir, e->d_name);
            nums[count] = num;
            count++;
        }
    }

    closedir(dp);
    
    /* Sort by descending number (try larger numbers first) */
    for (int i = 0; i < count - 1; i++) {
        for (int j = i + 1; j < count; j++) {
            if (nums[i] < nums[j]) {
                int tmp_num = nums[i];
                nums[i] = nums[j];
                nums[j] = tmp_num;
                
                char *tmp_path = paths[i];
                paths[i] = paths[j];
                paths[j] = tmp_path;
            }
        }
    }
    
    *count_out = count;
    return count > 0 ? paths : NULL;
}

/*-------------------------------------
    CSV loading
--------------------------------------*/
double *load_J_from_dir(const char *dir,
                        const char *file_type,
                        int generation,
                        int *N_out)
{
    int file_count = 0;
    char **paths = find_all_csv(dir, file_type, &file_count);
    if (!paths) {
        fprintf(stderr, "[ERROR] No %s_J*.csv found in %s\n", file_type, dir);
        return NULL;
    }

    /* Try multiple files in order (starting from larger numbers) */
    for (int file_idx = 0; file_idx < file_count; file_idx++) {
        char *path = paths[file_idx];
        printf("[INFO] Trying to load J from: %s\n", path);

        FILE *fp = fopen(path, "r");
        if (!fp) {
            perror("[ERROR] fopen");
            continue;
        }

        char line[65536];

        /* Skip header row */
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "[ERROR] Empty J file: %s\n", path);
            fclose(fp);
            continue;
        }

        char *target_line = NULL;
        int last_gen = -1;

        /* Search for the target generation (-1 means the last row) */
        while (fgets(line, sizeof(line), fp)) {
            char *tmp = strdup(line);
            char *tok = strtok(tmp, ",");
            if (!tok) {
                free(tmp);
                continue;
            }

            int gen = atoi(tok);
            free(tmp);

            if (generation == -1) {
                last_gen = gen;
                free(target_line);
                target_line = strdup(line);
            } else if (gen == generation) {
                free(target_line);
                target_line = strdup(line);
                break;
            }
        }
        fclose(fp);

        if (!target_line) {
            fprintf(stderr, "[WARN] generation=%d not found in %s, trying next file...\n",
                    generation, path);
            continue;  /* Try the next file */
        }

        /* Count elements excluding the generation column */
        int count = 0;
        char *tmp = strdup(target_line);
        strtok(tmp, ","); /* generation skip */
        while (strtok(NULL, ",")) {
            count++;
        }
        free(tmp);

        int N = (int)sqrt((double)count);
        if (N * N != count) {
            fprintf(stderr, "[ERROR] J size is not a perfect square: count=%d\n", count);
            free(target_line);
            continue;
        }

        *N_out = N;
        double *J = (double *)malloc(sizeof(double) * (size_t)N * (size_t)N);
        if (!J) {
            fprintf(stderr, "[ERROR] malloc J failed\n");
            free(target_line);
            continue;
        }

        /* Read the actual values */
        tmp = strdup(target_line);
        char *tok = strtok(tmp, ","); /* generation skip */
        (void)tok;

        for (int i = 0; i < count; i++) {
            char *v = strtok(NULL, ",");
            if (!v) {
                fprintf(stderr, "[ERROR] CSV parse error (short line)\n");
                free(tmp);
                free(target_line);
                free(J);
                J = NULL;
                break;
            }
            J[i] = atof(v);  /* Should be -1, 0, or 1, but keep it as double */
        }

        free(tmp);
        free(target_line);

        if (J) {
            printf("[INFO] Loaded J from %s: N=%d\n", path, N);
            /* Free memory */
            for (int i = 0; i < file_count; i++) {
                free(paths[i]);
            }
            return J;
        }
    }

    /* Not found in any file */
    fprintf(stderr, "[ERROR] generation=%d not found in any %s_J*.csv files\n",
            generation, file_type);
    for (int i = 0; i < file_count; i++) {
        free(paths[i]);
    }
    return NULL;
}
