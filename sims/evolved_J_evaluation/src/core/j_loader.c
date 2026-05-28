#include "j_loader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

/*-------------------------------------
  best_J*.csv / worst_J*.csv を dir から探す
  複数ある場合は配列で返す（番号降順でソート）
--------------------------------------*/
static char **find_all_csv(const char *dir, const char *file_type, int *count_out)
{
    DIR *dp = opendir(dir);
    if (!dp) return NULL;

    static char *paths[10];  /* 最大10ファイル */
    static int nums[10];     /* 各ファイルの番号 */
    int count = 0;

    struct dirent *e;

    /* 例: best_J,sigma=...csv / worst_J,...csv を探す */
    while ((e = readdir(dp)) && count < 10) {
        if (strstr(e->d_name, file_type) && strstr(e->d_name, "_J")) {
            // (N)形式の番号を抽出
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
    
    /* 番号降順でソート（大きい番号から試す） */
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
  CSV 読み込み
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

    /* 複数のファイルを順番に試す（番号が大きいものから） */
    for (int file_idx = 0; file_idx < file_count; file_idx++) {
        char *path = paths[file_idx];
        printf("[INFO] Trying to load J from: %s\n", path);

        FILE *fp = fopen(path, "r");
        if (!fp) {
            perror("[ERROR] fopen");
            continue;
        }

        char line[65536];

        /* ヘッダ行スキップ */
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "[ERROR] Empty J file: %s\n", path);
            fclose(fp);
            continue;
        }

        char *target_line = NULL;
        int last_gen = -1;

        /* 目的の generation を探す（-1 なら最終行） */
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
            continue;  /* 次のファイルを試す */
        }

        /* generation カラムを除いた要素数を数える */
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

        /* 実際の値を読み込む */
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
            J[i] = atof(v);  /* -1,0,1 のはずだが double にしておく */
        }

        free(tmp);
        free(target_line);

        if (J) {
            printf("[INFO] Loaded J from %s: N=%d\n", path, N);
            /* メモリ解放 */
            for (int i = 0; i < file_count; i++) {
                free(paths[i]);
            }
            return J;
        }
    }

    /* すべてのファイルで見つからなかった */
    fprintf(stderr, "[ERROR] generation=%d not found in any %s_J*.csv files\n",
            generation, file_type);
    for (int i = 0; i < file_count; i++) {
        free(paths[i]);
    }
    return NULL;
}
