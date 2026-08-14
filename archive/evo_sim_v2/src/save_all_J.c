#include "save_all_J.h"
#include <stdio.h>
#include <stdlib.h>

// J_list: 個体数×N×Nの配列（double型）
// pop_size: 個体数
// N: 遺伝子数
// out_path: 保存先ファイルパス
void save_all_J(const double *J_list, int pop_size, int N, const char *out_path)
{
    FILE *fp = fopen(out_path, "w");
    if (!fp)
    {
        fprintf(stderr, "[ERROR] Cannot open %s\n", out_path);
        return;
    }
    // ヘッダー
    fprintf(fp, "individual_id");
    for (int i = 0; i < N * N; ++i)
    {
        fprintf(fp, ",J%d", i);
    }
    fprintf(fp, "\n");
    // 各個体
    for (int ind = 0; ind < pop_size; ++ind)
    {
        fprintf(fp, "%d", ind);
        for (int i = 0; i < N * N; ++i)
        {
            fprintf(fp, ",%.8f", J_list[ind * N * N + i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
