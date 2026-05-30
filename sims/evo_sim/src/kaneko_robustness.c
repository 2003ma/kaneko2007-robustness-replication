#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#include "save_all_J.h"
#ifdef _OPENMP //openmpが利用可能な環境であれば並列化を有効にする
#include <omp.h>
#endif
#ifndef OMP_CHUNK //並列化のチャンクサイズ（デフォルトは1、必要に応じて変更）
#define OMP_CHUNK 1
#endif

/* create a unique filename by appending (n) before extension if file exists */
static void make_unique_path(const char *dir, const char *base, const char *ext, char *out, size_t out_size)
{
    snprintf(out, out_size, "%s/%s%s", dir, base, ext);
    if (access(out, F_OK) != 0)
        return;
    for (int n = 1; n < 10000; ++n)
    {
        snprintf(out, out_size, "%s/%s(%d)%s", dir, base, n, ext);
        if (access(out, F_OK) != 0)
            return;
    }
}

// J_list: 個体数×N×Nの配列（double型）
// pop_size: 個体数
// N: 遺伝子数
// out_path: 保存先ファイルパス

/*=============================
  乱数（XorShift32 + Box-Muller）
===============================*/
typedef struct
{
    uint32_t state; // XorShift32の内部状態（次の乱数を生成するための種）
    int has_spare;  // Box-Muller法で生成した2つ目の乱数が保存されているか（1=あり, 0=なし）
    double spare;   // Box-Muller法で生成した2つ目の正規乱数（次回使用のため保存）
} RNG;

/**
 * XorShift32: 高速な疑似乱数生成アルゴリズム
 * ビット演算（XOR, シフト）だけで32ビット整数の乱数を生成
 */
static inline uint32_t xorshift32(RNG *r)
{
    uint32_t x = r->state;
    x ^= x << 13; // 左に13ビットシフトしてXOR
    x ^= x >> 17; // 右に17ビットシフトしてXOR
    x ^= x << 5;  // 左に5ビットシフトしてXOR
    r->state = x; // 次回のために状態を更新
    return x;     // 0〜2^32-1の整数を返す
}

/**
 * urand: 一様分布の乱数 [0, 1) を生成
 * XorShift32の出力（0〜2^32-1）を0〜1の範囲に正規化
 */
static inline double urand(RNG *r) { return (xorshift32(r) + 1.0) / 4294967297.0; }

/**
 * nrand: 標準正規分布 N(0,1) に従う乱数を生成
 * Box-Muller法: 一様乱数2つから正規乱数2つを生成
 * - 1回の計算で2つ生成できるので、1つは次回用に保存して効率化
 */
static inline double nrand(RNG *r)
{
    // 前回生成した2つ目の乱数が残っていればそれを返す（高速）
    if (r->has_spare)
    {
        r->has_spare = 0;
        return r->spare;
    }

    // Box-Muller法で2つの正規乱数を生成
    double u1 = urand(r), u2 = urand(r);       // 一様乱数を2つ取得
    double rad = sqrt(-2.0 * log(u1 + 1e-16)); // 極座標の半径（1e-16でlog(0)を回避）
    double z0 = rad * cos(2.0 * M_PI * u2);    // 1つ目の正規乱数（今回返す）
    r->spare = rad * sin(2.0 * M_PI * u2);     // 2つ目の正規乱数（次回用に保存）
    r->has_spare = 1;                          // 予備の乱数があることを記録
    return z0;
}

/**
 * rng_seed: 乱数生成器を初期化
 * @param r 初期化するRNG構造体
 * @param seed シード値（0の場合はデフォルト値を使用）
 *
 * XorShift32の内部状態とBox-Muller用の予備乱数をリセット
 */
static inline void rng_seed(RNG *r, uint32_t seed)
{
    if (seed == 0)
        seed = 2463534242u; // seed=0を避けるためのデフォルト値
    r->state = seed;        // XorShift32の初期状態を設定
    r->has_spare = 0;       // Box-Muller用の予備乱数をクリア
    r->spare = 0.0;
}

/**
 * splitmix32: シード値を変換・分散させる
 * @param x 元のシード値
 * @return 変換後のシード値
 *
 * 並列処理で各スレッドに異なるシードを割り当てるために使用
 * 同じシードから始めると全スレッドが同じ乱数列を生成してしまうので、
 * この関数で元のシードを変換し、スレッドごとに異なる乱数列を得る
 */
static inline uint32_t splitmix32(uint32_t x)
{
    x += 0x9E3779B9u;                  // 黄金比に基づく定数を加算
    x = (x ^ (x >> 16)) * 0x85EBCA6Bu; // ビット混合（上位ビットを下位に影響）
    x = (x ^ (x >> 13)) * 0xC2B2AE35u; // さらにビット混合
    x ^= x >> 16;                      // 最終的な分散
    return x;
}

/*=============================
  構造体：ネットワーク個体
===============================*/

/**
 * RowIdx: 行列の各行における非ゼロ要素の高速アクセス用インデックス
 *
 * J行列は{-1, 0, +1}しか持たないため、各行について+1と-1の列番号を
 * 別々に保存しておくことで、内積計算を高速化できる
 * （乗算不要で加算と減算だけで済む）
 */
typedef struct
{
    int *pos;    // J[i][j]=+1 となる列jのインデックス配列
    int *neg;    // J[i][j]=-1 となる列jのインデックス配列
    int npos;    // +1の要素数（posの実際の使用数）
    int nneg;    // -1の要素数（negの実際の使用数）
    int cap_pos; // pos配列の容量（動的拡張用）
    int cap_neg; // neg配列の容量（動的拡張用）
} RowIdx;

/**
 * Individual: 進化する個体（GRNネットワーク）
 *
 * 各個体は遺伝子制御ネットワーク（J行列）と適応度を持つ
 */
typedef struct
{
    int N;          // 遺伝子数（ネットワークのノード数）
    int8_t *J;      // N×N の相互作用行列（-1, 0, +1の3値のみ）
    RowIdx *rows;   // N行分の高速インデックス（各行の+1, -1要素の位置）
    double fitness; // 適応度（0が最大、負の値が悪い）
    double V_ip;    // 同一遺伝子型の分散（isogenic variance proxy）
} Individual;

// J(i,j) に簡単にアクセスするマクロ
#define JAT(ind, i, j) ((ind)->J[(size_t)(i) * (ind)->N + (size_t)(j)])

/*=============================
  ユーティリティ
*/
// 現在と違う-1,0,1を返す
static inline int8_t mutate_symbol(int8_t cur, RNG *rng)
{
    static const int8_t c[3] = {-1, 0, 1};
    int idx = (int)(urand(rng) * 3.0);
    int8_t v = c[idx];
    if (v == cur)
    {
        idx = (idx + 1 + (int)(urand(rng) * 2.0)) % 3;
        v = c[idx];
    }
    return v;
}

// メモリ解放
static void rowidx_free(RowIdx *ri)
{
    if (!ri)
        return;
    free(ri->pos);
    free(ri->neg);
    ri->pos = NULL;
    ri->neg = NULL;
    ri->npos = ri->nneg = ri->cap_pos = ri->cap_neg = 0;
}

/**
 * build_row_index: J行列から高速アクセス用のインデックスを構築
 * @param ind 個体（J行列とrowsを含む）
 *
 * 【目的】内積計算の高速化
 * J行列の各行について、+1と-1の要素がどの列にあるかを記録する。
 * これにより、h = Σ J[i][j] * x[j] の計算を
 * 「+1の列のxを足す、-1の列のxを引く」だけで済ませられる（乗算不要）
 *
 * 【処理の流れ】
 * 1. 各行iについて、+1と-1の個数をカウント
 * 2. その個数分の配列を確保
 * 3. もう一度走査して、+1/-1の列番号を配列に詰める
 */
static void build_row_index(Individual *ind)
{
    int N = ind->N;
    // N行分のRowIdx構造体を確保
    ind->rows = (RowIdx *)calloc((size_t)N, sizeof(RowIdx));
    if (!ind->rows)
    {
        fprintf(stderr, "[FATAL] rows alloc failed\n");
        exit(1);
    }

    for (int i = 0; i < N; ++i)
    {
        // ===== ステップ1: 各行のJ[i][j]について+1と-1の個数をカウント =====
        int npos = 0, nneg = 0;
        for (int j = 0; j < N; ++j)
        {
            int8_t v = JAT(ind, i, j);
            npos += (v == 1);  // +1なら1を加算
            nneg += (v == -1); // -1なら1を加算
        }

        // ===== ステップ2: カウント結果に基づいて配列を確保 =====
        RowIdx *ri = &ind->rows[i];
        ri->cap_pos = ri->npos = npos; // +1の個数と容量を設定
        ri->cap_neg = ri->nneg = nneg; // -1の個数と容量を設定
        if (npos > 0)
            ri->pos = (int *)malloc(sizeof(int) * npos); // +1用の配列確保
        if (nneg > 0)
            ri->neg = (int *)malloc(sizeof(int) * nneg); // -1用の配列確保
        if ((npos > 0 && !ri->pos) || (nneg > 0 && !ri->neg))
        {
            fprintf(stderr, "[FATAL] rowidx alloc failed\n");
            exit(1);
        }

        // ===== ステップ3: 再度走査して+1/-1の列番号を配列に格納 =====
        int ap = 0, an = 0; // posとneg配列のインデックス
        for (int j = 0; j < N; ++j)
        {
            int8_t v = JAT(ind, i, j);
            if (v == 1)
                ri->pos[ap++] = j; // J[i][j]=+1 なら列番号jを記録
            else if (v == -1)
                ri->neg[an++] = j; // J[i][j]=-1 なら列番号jを記録
        }
    }
}

/**
 * init_random_individual : 個体をランダムに初期化
 * @param ind 個体
 * @param N 遺伝子数
 * @param p_edge 相互作用の存在確率
 * @param rng 乱数生成器
 */
static void init_random_individual(Individual *ind, int N, double p_edge, RNG *rng)
{
    ind->N = N;
    size_t sz = (size_t)N * (size_t)N;
    ind->J = (int8_t *)malloc(sz);
    if (!ind->J)
    {
        fprintf(stderr, "malloc J failed\n");
        exit(1);
    }
    // Jを初期化
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // 確率p_edgeで相互作用を埋める
            double u = urand(rng);
            if (u < p_edge)
                // 1か-1をランダムに生成
                JAT(ind, i, j) = (urand(rng) < 0.5) ? 1 : -1;
            else
                JAT(ind, i, j) = 0;
        }
    }
    build_row_index(ind);
    ind->fitness = -1e9;
    ind->V_ip = 0.0;
}

/**
 * clone_individual : 個体をクローン（深いコピー）
 * @param src コピー元の個体
 * @return コピーされた新しい個体
 *
 * J行列とRowIdxインデックスを新たに確保してコピーする
 */
static Individual clone_individual(const Individual *src)
{
    Individual dst;
    dst.N = src->N;
    size_t sz = (size_t)dst.N * (size_t)dst.N;
    dst.J = (int8_t *)malloc(sz);
    if (!dst.J)
    {
        fprintf(stderr, "malloc clone J failed\n");
        exit(1);
    }
    memcpy(dst.J, src->J, sz);
    dst.fitness = src->fitness;
    dst.V_ip = src->V_ip;
    /* rows 再構築（Jから作る）*/
    dst.rows = NULL;
    build_row_index(&dst);
    return dst;
}

/**
 * mutate_one_edge : 個体のJ行列の1つのエッジを突然変異（参考用にコメントアウト）
 * @param ind 突然変異させる個体
 * @param rng 乱数生成器
 *
 * ランダムに(i,j)を選び、J[i][j]の値を-1,0,1のいずれかに変更
 * 変更後の値が元と同じ場合は何もしない
 * RowIdxインデックスも更新する
 */
static void mutate_one_edge(Individual *ind, RNG *rng)
{
    int N = ind->N;
    // ランダムに(i,j)を選ぶ
    int i = (int)(urand(rng) * N);
    int j = (int)(urand(rng) * N);
    int8_t cur = JAT(ind, i, j); // 現在のJijの値
    int8_t nxt = mutate_symbol(cur, rng);  //現在と違う値を返す
    if (cur == nxt)  //mutate_symbolで同じ値が返されることはないはずだが、念のためチェック
        return;
    JAT(ind, i, j) = nxt;
    RowIdx *ri = &ind->rows[i];

    /* ===== ステップ1: 古い値(cur)をRowIdxから削除 ===== */
    if (cur == 1)
    {
        // curが+1だった場合、pos配列から列番号jを探して削除
        for (int t = 0; t < ri->npos; ++t)
            if (ri->pos[t] == j)
            {
                // 削除方法: 配列の最後の要素で上書き（順序は気にしない）
                // --ri->npos でサイズを1減らしながら、最後の要素を取得
                ri->pos[t] = ri->pos[--ri->npos];
                break;
            }
    }
    else if (cur == -1)
    {
        // curが-1だった場合、neg配列から列番号jを探して削除
        for (int t = 0; t < ri->nneg; ++t)
            if (ri->neg[t] == j)
            {
                // 同様に、最後の要素で上書きして削除
                ri->neg[t] = ri->neg[--ri->nneg];
                break;
            }
    }
    // cur==0の場合は何もしない（元々インデックスに含まれていない）

    /* ===== ステップ2: 新しい値(nxt)をRowIdxに追加 ===== */
    if (nxt == 1)
    {
        // nxtが+1の場合、pos配列に列番号jを追加
        // 容量が足りなければ拡張（2倍に増やす、初期値は4）
        if (ri->npos == ri->cap_pos)
        {
            ri->cap_pos = ri->cap_pos ? ri->cap_pos * 2 : 4;
            ri->pos = (int *)realloc(ri->pos, sizeof(int) * ri->cap_pos);
        }
        // 配列の末尾に追加してサイズをインクリメント
        ri->pos[ri->npos++] = j;
    }
    else if (nxt == -1)
    {
        // nxtが-1の場合、neg配列に列番号jを追加
        // 同様に容量チェックと拡張
        if (ri->nneg == ri->cap_neg)
        {
            ri->cap_neg = ri->cap_neg ? ri->cap_neg * 2 : 4;
            ri->neg = (int *)realloc(ri->neg, sizeof(int) * ri->cap_neg);
        }
        // 配列の末尾に追加
        ri->neg[ri->nneg++] = j;
    }
    // nxt==0の場合は何もしない（インデックスに登録不要）
}

// 個体のメモリ解放
static void free_individual(Individual *ind)
{
    if (!ind)
        return;
    if (ind->rows)
    {
        for (int i = 0; i < ind->N; ++i)
            rowidx_free(&ind->rows[i]);
        free(ind->rows);
        ind->rows = NULL;
    }
    free(ind->J);
    ind->J = NULL;
}

/*=============================
  数値ダイナミクス
===============================*/
#define TANH(x) tanh(x)

// 行インデックスに基づき、x ベクトルとの内積を計算
// j>=k（内部遺伝子）からのみ入力を受ける
// RowIdx *riはi行目のインデックス r行目インデックスのpos, negを使う
static inline double dot_row_idx(const RowIdx *ri, const double *restrict x, int k)
{
    double h = 0.0;
    // +1の要素について、j>=kのもののみ加算
    for (int t = 0; t < ri->npos; ++t)
    {
        int j = ri->pos[t];
        if (j >= k)
            h += x[j];
    }
    // -1の要素について、j>=kのもののみ減算
    for (int t = 0; t < ri->nneg; ++t)
    {
        int j = ri->neg[t];
        if (j >= k)
            h -= x[j];
    }
    return h;
}

/* 1 回の試行（Lの中の1回）を評価  */
static inline double simulate_once_and_score_buf(
    const Individual *ind, const int *outputs, int k,
    double beta, double sigma, double dt, int relax_steps, int meas_steps,
    double *restrict x, double *restrict x_next, RNG *rng)
{
    const int N = ind->N;
    const double sqrt_dt = sqrt(dt);
    // 全てのスピンを-1に初期化
    for (int i = 0; i < N; ++i)
        x[i] = -1.0;

    // 緩和ステップ
    for (int t = 0; t < relax_steps; ++t)
    {
        // 全遺伝子の状態更新
        for (int i = 0; i < N; ++i)
        {
            const RowIdx *ri = &ind->rows[i];
            double h = dot_row_idx(ri, x, k);
            double drift = TANH(beta * h) - x[i];
            double noise = sigma * sqrt_dt * nrand(rng);
            double xi = x[i] + drift * dt + noise;
            x_next[i] = xi;
        }
        double *tmp = x;
        x = x_next;
        x_next = tmp;
    }

    double sum_on = 0.0;
    for (int t = 0; t < meas_steps; ++t)
    {
        for (int i = 0; i < N; ++i)
        {
            const RowIdx *ri = &ind->rows[i];
            double h = dot_row_idx(ri, x, k);
            double drift = TANH(beta * h) - x[i];
            double noise = sigma * sqrt_dt * nrand(rng);
            double xi = x[i] + drift * dt + noise;
            x_next[i] = xi;
        }
        double *tmp = x;
        x = x_next;
        x_next = tmp;

        int on_now = 0;
        for (int a = 0; a < k; ++a)
            // 発現したかどうかをカウント
            on_now += (x[outputs[a]] > 0.0);
        sum_on += (double)on_now;
    }

    double mean_on = sum_on / (double)meas_steps; // 出力遺伝子のオン数の平均
    return mean_on - (double)k;                   // 0に近いほどいい
}

/* 同じ遺伝子でL回の試行を評価
 */
static inline void evaluate_fitness_buf(
    const Individual *ind, const int *outputs, int k,
    double beta, double sigma, double dt, int relax_steps, int meas_steps,
    int L, RNG *rng, double *v_ip, double *mean_fitness,
    double *restrict x, double *restrict x_next)
{
    double one_sim_fitness = 0.0;
    double squared_one_sim_fitness = 0.0;
    for (int r = 0; r < L; ++r)
    {
        double f = simulate_once_and_score_buf(ind, outputs, k, beta, sigma, dt,
                                               relax_steps, meas_steps, x, x_next, rng);
        one_sim_fitness += f;
        squared_one_sim_fitness += f * f;


    }
    double mean = one_sim_fitness / (double)L;
    *mean_fitness = mean;
    *v_ip = (squared_one_sim_fitness / (double)L) - mean * mean;
}

/**
  個体を適応度降順に比較する関数
  @param a 比較する個体1へのポインタ
  @param b 比較する個体2へのポインタ
*/
static int cmp_ind_desc(const void *a, const void *b)
{
    const Individual *x = (const Individual *)a, *y = (const Individual *)b;
    return (x->fitness > y->fitness) ? -1 : ((x->fitness < y->fitness) ? 1 : 0);
}

/**
 * シミュレーションパラメータ構造体
 */
typedef struct
{
    int N, k, pop, gens;
    double beta, sigma, dt;
    int relax, meas, L;
    double pedge, elit;
    double mu; // 突然変異率（0.0〜1.0、確率でエッジを変異）
    uint32_t seed;
    const char *load_J_path; // J行列をロードするCSVファイルパス（NULLならランダム初期化）
} Params;

/* J保存用ディレクトリ（外部から指定） */
static const char *g_J_save_dir = NULL;

void set_J_save_dir(const char *dir)
{
    g_J_save_dir = dir;
}

/* 出力遺伝子を選ぶ（ここでは単純に最初のk個を選ぶ） */
static void choose_outputs(int N, int k, int *outputs, RNG *rng)
{
    (void)rng;
    for (int a = 0; a < k; ++a)
        outputs[a] = a;
}

/**
 * ファイル名から世代番号を抽出する
 * @param filepath ファイルパス（例: "results/ver3/sigma_0.020/gen_20_all_J.csv"）
 * @return 世代番号（見つからない場合は0）
 *
 * "gen_XX_" のパターンを探して数値部分を抽出
 */
static int extract_generation_from_path(const char *filepath)
{
    if (!filepath)
        return 0;
    
    // ファイル名部分を探す（最後の'/'以降）
    const char *filename = strrchr(filepath, '/');
    if (filename)
        filename++; // '/'の次の文字から
    else
        filename = filepath; // パス区切りがない場合は全体
    
    // "gen_" を探す
    const char *gen_pos = strstr(filename, "gen_");
    if (!gen_pos)
        return 0;
    
    // "gen_"の後の数値を取得
    int gen_num = 0;
    if (sscanf(gen_pos, "gen_%d", &gen_num) == 1)
        return gen_num;
    
    return 0;
}

/**
 * CSVファイルから集団を読み込む
 * @param filepath CSVファイルのパス（gen_XX_all_J.csv形式）
 * @param pop 個体配列（事前に確保されている必要がある）
 * @param expected_pop 期待する個体数
 * @param expected_N 期待する遺伝子数
 * @return 成功したら0、失敗したら-1
 */
static int load_population_from_csv(const char *filepath, Individual *pop, int expected_pop, int expected_N)
{
    FILE *fp = fopen(filepath, "r");
    if (!fp)
    {
        fprintf(stderr, "[ERROR] Cannot open %s for loading J matrices\n", filepath);
        return -1;
    }

    char line[65536]; // 十分大きなバッファ（N=64なら4096要素程度）

    // ヘッダー行をスキップ（id,fitness,v_ip,J_0_0,J_0_1,...）
    if (!fgets(line, sizeof(line), fp))
    {
        fprintf(stderr, "[ERROR] Cannot read header from %s\n", filepath);
        fclose(fp);
        return -1;
    }

    // 各個体を読み込む
    int loaded_count = 0;
    while (fgets(line, sizeof(line), fp) && loaded_count < expected_pop)
    {
        Individual *ind = &pop[loaded_count];
        ind->N = expected_N;
        size_t sz = (size_t)expected_N * (size_t)expected_N;
        ind->J = (int8_t *)malloc(sz);
        if (!ind->J)
        {
            fprintf(stderr, "[ERROR] malloc failed for J matrix\n");
            fclose(fp);
            return -1;
        }

        // 行をパース: id,fitness,v_ip,J_0_0,J_0_1,...
        char *token = strtok(line, ",");
        if (!token)
        {
            fprintf(stderr, "[ERROR] Invalid CSV format (no id)\n");
            free(ind->J);
            fclose(fp);
            return -1;
        }

        // fitness
        token = strtok(NULL, ",");
        if (!token)
        {
            fprintf(stderr, "[ERROR] Invalid CSV format (no fitness)\n");
            free(ind->J);
            fclose(fp);
            return -1;
        }
        ind->fitness = atof(token);

        // v_ip
        token = strtok(NULL, ",");
        if (!token)
        {
            fprintf(stderr, "[ERROR] Invalid CSV format (no v_ip)\n");
            free(ind->J);
            fclose(fp);
            return -1;
        }
        ind->V_ip = atof(token);

        // J行列（N*N個）
        for (int i = 0; i < expected_N * expected_N; ++i)
        {
            token = strtok(NULL, ",");
            if (!token)
            {
                fprintf(stderr, "[ERROR] Invalid CSV format (insufficient J elements at index %d)\n", i);
                free(ind->J);
                fclose(fp);
                return -1;
            }
            ind->J[i] = (int8_t)atoi(token);
        }

        // RowIdxを構築
        ind->rows = NULL;
        build_row_index(ind);

        ++loaded_count;
    }

    fclose(fp);

    if (loaded_count != expected_pop)
    {
        fprintf(stderr, "[WARNING] Expected %d individuals but loaded %d from %s\n",
                expected_pop, loaded_count, filepath);
    }

    printf("[INFO] Loaded %d individuals from %s\n", loaded_count, filepath);
    return 0;
}

/* スレッド私有ワークバッファ */
typedef struct
{
    double *x, *x_next;
    int capN;
} ThreadBuf;

// スレッド私有ワークバッファを N に合わせて確保・拡張
static void tb_ensure(ThreadBuf *tb, int N)
{
    if (tb->capN >= N)
        return;
    free(tb->x);
    free(tb->x_next);
    void *p1 = NULL, *p2 = NULL;
    if (posix_memalign(&p1, 64, (size_t)N * sizeof(double)) != 0)
        p1 = NULL;
    if (posix_memalign(&p2, 64, (size_t)N * sizeof(double)) != 0)
        p2 = NULL;
    tb->x = (double *)p1;
    tb->x_next = (double *)p2;
    if (!tb->x || !tb->x_next)
    {
        fprintf(stderr, "[FATAL] aligned buffer alloc failed (N=%d)\n", N);
        exit(1);
    }
    tb->capN = N;
}

/* シミュレーションを実行する関数
   @param prm_in シミュレーションパラメータへのポインタ

*/
void run_kaneko_robustness_simulation(const Params *prm_in)
{
    Params prm = *prm_in;
    if (prm.k <= 0 || prm.k > prm.N)
    {
        fprintf(stderr, "k must be in [1,N]\n");
        return;
    }
    if (prm.pop < 2)
    {
        fprintf(stderr, "pop must be >= 2\n");
        return;
    }
    // エリート率を0〜0.9の範囲に制限
    if (prm.elit < 0.0)
        prm.elit = 0.0;
    if (prm.elit > 0.9)
        prm.elit = 0.9;

    // ===== 初期化 =====
    // マスター乱数生成器を初期化（突然変異などで使用）
    RNG master;
    rng_seed(&master, prm.seed);

    // 出力遺伝子のインデックスを格納する配列を確保
    int *outputs = (int *)malloc((size_t)prm.k * sizeof(int));
    choose_outputs(prm.N, prm.k, outputs, &master); // 最初のk個を出力遺伝子に指定

    // 個体群を確保（現世代と次世代の2つ）
    Individual *pop = (Individual *)malloc((size_t)prm.pop * sizeof(Individual));
    Individual *next = (Individual *)malloc((size_t)prm.pop * sizeof(Individual));
    if (!pop || !next)
    {
        fprintf(stderr, "malloc pop failed\n");
        free(outputs);
        return;
    }

    // 初期個体群の設定（CSVからロード or ランダム生成）
    int start_gen = 0; // 開始世代番号
    if (prm.load_J_path)
    {
        // CSVファイルから読み込み
        if (load_population_from_csv(prm.load_J_path, pop, prm.pop, prm.N) != 0)
        {
            fprintf(stderr, "[ERROR] Failed to load population from %s\n", prm.load_J_path);
            free(pop);
            free(next);
            free(outputs);
            return;
        }
        // ファイル名から世代番号を抽出
        start_gen = extract_generation_from_path(prm.load_J_path);
        printf("[INFO] Starting evolution from loaded population (file: %s, generation: %d)\n", 
               prm.load_J_path, start_gen);
    }
    else
    {
        // ランダムに生成
        for (int i = 0; i < prm.pop; ++i)
            init_random_individual(&pop[i], prm.N, prm.pedge, &master);
        printf("[INFO] Starting evolution from random initial population\n");
    }

    // エリート個体数を計算（上位何個体を親とするか）
    int elit_n = (int)floor(prm.elit * prm.pop);

    setvbuf(stdout, NULL, _IOLBF, 0); /* 行バッファリング */
    printf("# Kaneko2007-like GRN evolution (C, more_fast_idx)\n");
    printf("# N=%d k=%d pop=%d gens=%d beta=%.3f sigma=%.4f dt=%.3f relax=%d meas=%d L=%d pedge=%.3f elit=%.2f mu=%.3f seed=%u\n",
           prm.N, prm.k, prm.pop, prm.gens, prm.beta, prm.sigma, prm.dt, prm.relax, prm.meas, prm.L, prm.pedge, prm.elit, prm.mu, prm.seed);
    printf("# generation,best,mean,worst,V_g,V_ip\n");

    /* ディレクトリ構造: results/verN/sigma_X.XX/ */
    /* 既存のバージョンをチェックして、同じパラメータセットがあるか確認 */
    (void)mkdir("results", 0755);

    char ver_dir[256];
    char params_path[512];
    int ver_num = 1;
    int found_matching_ver = 0;

    /* 既存のバージョンをチェック */
    for (int v = 1; v < 100; ++v)
    {
        snprintf(ver_dir, sizeof(ver_dir), "results/ver%d", v);
        snprintf(params_path, sizeof(params_path), "%s/params.txt", ver_dir);

        FILE *check_fp = fopen(params_path, "r");
        if (check_fp)
        {
            /* params.txtを読んでパラメータが一致するかチェック */
            char line[256];
            int match_pop = 0, match_L = 0, match_meas = 0, match_relax = 0, match_seed = 0;

            while (fgets(line, sizeof(line), check_fp))
            {
                int pop_val, L_val, meas_val, relax_val;
                unsigned int seed_val;

                if (sscanf(line, "Population size (pop): %d", &pop_val) == 1 && pop_val == prm.pop)
                    match_pop = 1;
                else if (sscanf(line, "Number of trials (L): %d", &L_val) == 1 && L_val == prm.L)
                    match_L = 1;
                else if (sscanf(line, "Measurement steps (meas): %d", &meas_val) == 1 && meas_val == prm.meas)
                    match_meas = 1;
                else if (sscanf(line, "Relaxation steps (relax): %d", &relax_val) == 1 && relax_val == prm.relax)
                    match_relax = 1;
                else if (sscanf(line, "Random seed: %u", &seed_val) == 1 && seed_val == prm.seed)
                    match_seed = 1;
            }
            fclose(check_fp);

            /* 全パラメータが一致したら、このバージョンを使用 */
            if (match_pop && match_L && match_meas && match_relax && match_seed)
            {
                ver_num = v;
                found_matching_ver = 1;
                break;
            }
        }
        else
        {
            /* このバージョン番号がまだ存在しない = 新規作成 */
            ver_num = v;
            break;
        }
    }

    /* バージョンディレクトリを作成 */
    snprintf(ver_dir, sizeof(ver_dir), "results/ver%d", ver_num);
    (void)mkdir(ver_dir, 0755);

    /* シグマディレクトリを作成 */
    char sigma_dir[256];
    snprintf(sigma_dir, sizeof(sigma_dir), "%s/sigma_%.3f", ver_dir, prm.sigma);
    (void)mkdir(sigma_dir, 0755);

    /* evo_sim_dataディレクトリを作成 */
    char evo_sim_data_dir[300];
    snprintf(evo_sim_data_dir, sizeof(evo_sim_data_dir), "%s/evo_sim_data", sigma_dir);
    (void)mkdir(evo_sim_data_dir, 0755);

    /* J保存ディレクトリをevo_sim_dataに設定 */
    set_J_save_dir(evo_sim_data_dir);

    /* params.txtを作成（新規バージョンの場合のみ） */
    if (!found_matching_ver)
    {
        snprintf(params_path, sizeof(params_path), "%s/params.txt", ver_dir);
        FILE *params_fp = fopen(params_path, "w");
        if (params_fp)
        {
            fprintf(params_fp, "Simulation Parameters for ver%d\n", ver_num);
            fprintf(params_fp, "================================\n");
            fprintf(params_fp, "Population size (pop): %d\n", prm.pop);
            fprintf(params_fp, "Number of trials (L): %d\n", prm.L);
            fprintf(params_fp, "Measurement steps (meas): %d\n", prm.meas);
            fprintf(params_fp, "Relaxation steps (relax): %d\n", prm.relax);
            fprintf(params_fp, "Random seed: %u\n", prm.seed);
            fprintf(params_fp, "\nVariable parameter:\n");
            fprintf(params_fp, "- sigma: (varies)\n");
            fclose(params_fp);
        }
    }

    /* データファイル */
    char summary_base[160];
    snprintf(summary_base, sizeof(summary_base), "data,sigma=%.4f,dt=%.3f,pop=%d,L=%d,meas=%d,relax=%d,seed=%u", prm.sigma, prm.dt, prm.pop, prm.L, prm.meas, prm.relax, prm.seed);
    char summary_path[512];
    make_unique_path(evo_sim_data_dir, summary_base, ".csv", summary_path, sizeof(summary_path));
    FILE *out_fp = fopen(summary_path, "w");
    if (out_fp)
    {
        /* summary header: append per-individual v_ip columns on the right */
        fprintf(out_fp, "generation,best,mean,worst,V_g,Vip");
        for (int __i = 0; __i < prm.pop; ++__i)
            fprintf(out_fp, ",v_ip_%d", __i);
        fprintf(out_fp, "\n");
        fflush(out_fp);
    }
    else
    {
        fprintf(stderr, "Warning: cannot open %s for writing. Output will only go to stdout\n", summary_path);
    }

    /* 最良個体のJを保存するファイル */
    char best_J_base[160];
    snprintf(best_J_base, sizeof(best_J_base), "best_J,sigma=%.4f,dt=%.3f,pop=%d,L=%d,meas=%d,relax=%d,seed=%u", prm.sigma, prm.dt, prm.pop, prm.L, prm.meas, prm.relax, prm.seed);
    char best_J_path[512];
    make_unique_path(evo_sim_data_dir, best_J_base, ".csv", best_J_path, sizeof(best_J_path));
    FILE *best_J_fp = fopen(best_J_path, "w");
    if (best_J_fp)
    {
        fprintf(best_J_fp, "generation");
        for (int i = 0; i < prm.N; ++i)
        {
            for (int j = 0; j < prm.N; ++j)
            {
                fprintf(best_J_fp, ",J_%d_%d", i, j);
            }
        }
        fprintf(best_J_fp, "\n");
        fflush(best_J_fp);
    }
    else
    {
        fprintf(stderr, "Warning: cannot open %s for writing. Best J will not be saved.\n", best_J_path);
    }

    /* 最悪個体のJを保存するファイル */
    char worst_J_base[160];
    snprintf(worst_J_base, sizeof(worst_J_base), "worst_J,sigma=%.4f,dt=%.3f,pop=%d,L=%d,meas=%d,relax=%d,seed=%u", prm.sigma, prm.dt, prm.pop, prm.L, prm.meas, prm.relax, prm.seed);
    char worst_J_path[512];
    make_unique_path(evo_sim_data_dir, worst_J_base, ".csv", worst_J_path, sizeof(worst_J_path));
    FILE *worst_J_fp = fopen(worst_J_path, "w");
    if (worst_J_fp)
    {
        fprintf(worst_J_fp, "generation");
        for (int i = 0; i < prm.N; ++i)
        {
            for (int j = 0; j < prm.N; ++j)
            {
                fprintf(worst_J_fp, ",J_%d_%d", i, j);
            }
        }
        fprintf(worst_J_fp, "\n");
        fflush(worst_J_fp);
    }
    else
    {
        fprintf(stderr, "Warning: cannot open %s for writing. Worst J will not be saved.\n", worst_J_path);
    }

    int nthreads = 1;
    RNG *thread_rng;
    ThreadBuf *tb_pool;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    thread_rng = (RNG *)malloc((size_t)nthreads * sizeof(RNG));
    tb_pool = (ThreadBuf *)calloc((size_t)nthreads, sizeof(ThreadBuf));
    for (int t = 0; t < nthreads; ++t)
    {
        uint32_t s = splitmix32(prm.seed ^ (uint32_t)t ^ 0xA5A5A5A5u);
        rng_seed(&thread_rng[t], s);
        tb_pool[t].x = tb_pool[t].x_next = NULL;
        tb_pool[t].capN = 0;
    }

/* ===== 常設並列チーム ===== */
#pragma omp parallel // 最大コアのスレッドを立ち上げる
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        ThreadBuf *tb = &tb_pool[tid];

        for (int g = 0; g < prm.gens; ++g)
        {
            int current_gen = start_gen + g; // 実際の世代番号

/* 個体評価：pop 方向に並列 */
#pragma omp for schedule(static, OMP_CHUNK)
            for (int i = 0; i < prm.pop; ++i)
            {
                tb_ensure(tb, pop[i].N);
                RNG rng = thread_rng[tid];
                rng.state ^= splitmix32((uint32_t)(current_gen * 2654435761u + i));

                double vip, fit;
                evaluate_fitness_buf(&pop[i], outputs, prm.k, prm.beta, prm.sigma, prm.dt, // ここのkがj>=kの入力を受け取る
                                     prm.relax, prm.meas, prm.L, &rng, &vip, &fit,
                                     tb->x, tb->x_next);
                pop[i].V_ip = vip;
                pop[i].fitness = fit;
                thread_rng[tid] = rng;
            }

            /* スレッド集約: ここでは single 節で共有累積を初期化するため、
               空の atomic 指示は不要（空文はコンパイラエラーになる）。
               以下の single 節に処理を移す。 */

#pragma omp single
            {
                /* single 節の先頭で共有累積を初期化 */
            }

/* 共有累積用の静的変数は使わず、以下の single 節で再計算するのが簡潔 */
#pragma omp barrier

/* 全個体の fitness を既に各要素に格納済み → single で集計とソート */
#pragma omp single
            {
                // Vgを計算
                //  集団の平均適応度を計算（統計出力用）
                double mean_fit_sum = 0.0, fitness_squared_sum = 0.0;
                for (int i = 0; i < prm.pop; ++i)
                {
                    mean_fit_sum += pop[i].fitness;
                    fitness_squared_sum += pop[i].fitness * pop[i].fitness;
                }
                // その世代の全遺伝子の適応度の平均
                double mean_fit = mean_fit_sum / (double)prm.pop;

                // 適応度降順にソート
                qsort(pop, (size_t)prm.pop, sizeof(Individual), cmp_ind_desc);
                double best = pop[0].fitness;
                double worst = pop[prm.pop - 1].fitness;
                double Vg = (fitness_squared_sum / (double)prm.pop) - mean_fit * mean_fit;
                double Vip_sum = 0.0;
                for (int _i = 0; _i < prm.pop; ++_i)
                    Vip_sum += pop[_i].V_ip;
                double Vip_mean = Vip_sum / (double)prm.pop;

                printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f\n", current_gen, best, mean_fit, worst, Vg, Vip_mean);

                if (out_fp)
                {
                    fprintf(out_fp, "%d,%.6f,%.6f,%.6f,%.6f,%.6f", current_gen, best, mean_fit, worst, Vg, Vip_mean);
                    for (int _pi = 0; _pi < prm.pop; ++_pi)
                        fprintf(out_fp, ",%.6f", pop[_pi].V_ip);
                    fprintf(out_fp, "\n");
                    fflush(out_fp);
                }

                /* 最良個体のJを保存 */
                if (best_J_fp)
                {
                    fprintf(best_J_fp, "%d", current_gen);
                    for (int i = 0; i < prm.N; ++i)
                    {
                        for (int j = 0; j < prm.N; ++j)
                        {
                            fprintf(best_J_fp, ",%d", (int)JAT(&pop[0], i, j));
                        }
                    }
                    fprintf(best_J_fp, "\n");
                    fflush(best_J_fp);
                }

                /* 最悪個体のJを保存 */
                if (worst_J_fp)
                {
                    fprintf(worst_J_fp, "%d", current_gen);
                    for (int i = 0; i < prm.N; ++i)
                    {
                        for (int j = 0; j < prm.N; ++j)
                        {
                            fprintf(worst_J_fp, ",%d", (int)JAT(&pop[prm.pop - 1], i, j));
                        }
                    }
                    fprintf(worst_J_fp, "\n");
                    fflush(worst_J_fp);
                }

                if (g == prm.gens - 1)
                {
                    /* 最終世代の全個体J行列を保存 (id,fitness,v_ip,J_00,... フォーマット) */
                    if (g_J_save_dir)
                    {
                        char final_J_path[512];
                        snprintf(final_J_path, sizeof(final_J_path), "%s/gen_%d_all_J_sigma_%.3f_dt%.3f.csv", g_J_save_dir, current_gen, prm.sigma, prm.dt);

                        FILE *all_fp = fopen(final_J_path, "w");
                        if (all_fp)
                        {
                            /* ヘッダ */
                            fprintf(all_fp, "id,fitness,v_ip");
                            for (int ii = 0; ii < prm.N; ++ii)
                                for (int jj = 0; jj < prm.N; ++jj)
                                    fprintf(all_fp, ",J_%d_%d", ii, jj);
                            fprintf(all_fp, "\n");

                            /* 各個体の行 */
                            for (int pi = 0; pi < prm.pop; ++pi)
                            {
                                fprintf(all_fp, "%d,%.6f,%.6f", pi, pop[pi].fitness, pop[pi].V_ip);
                                for (int j = 0; j < prm.N * prm.N; ++j)
                                    fprintf(all_fp, ",%d", (int)pop[pi].J[j]);
                                fprintf(all_fp, "\n");
                            }
                            fclose(all_fp);
                            printf("[INFO] Generation %d (final) J matrices saved to %s\n", current_gen, final_J_path);
                        }
                        else
                        {
                            fprintf(stderr, "Warning: cannot open %s for writing\n", final_J_path);
                        }
                    }
                }
                else
                {
                    /* 20世代ごとにJ行列を保存 */
                    if (g_J_save_dir && current_gen % 20 == 0)
                    {
                        char gen_J_path[512];
                        snprintf(gen_J_path, sizeof(gen_J_path), "%s/gen_%d_all_J_sigma_%.3f_dt%.3f.csv", g_J_save_dir, current_gen, prm.sigma, prm.dt);

                    FILE *all_fp = fopen(gen_J_path, "w");
                    if (all_fp)
                    {
                        /* ヘッダ */
                        fprintf(all_fp, "id,fitness,v_ip");
                        for (int ii = 0; ii < prm.N; ++ii)
                            for (int jj = 0; jj < prm.N; ++jj)
                                fprintf(all_fp, ",J_%d_%d", ii, jj);
                        fprintf(all_fp, "\n");

                        /* 各個体の行 */
                        for (int pi = 0; pi < prm.pop; ++pi)
                        {
                            fprintf(all_fp, "%d,%.6f,%.6f", pi, pop[pi].fitness, pop[pi].V_ip);
                            for (int j = 0; j < prm.N * prm.N; ++j)
                                fprintf(all_fp, ",%d", (int)pop[pi].J[j]);
                            fprintf(all_fp, "\n");
                        }
                        fclose(all_fp);
                        printf("[INFO] Generation %d J matrices saved to %s\n", current_gen, gen_J_path);
                    }
                    else
                    {
                        fprintf(stderr, "Warning: cannot open %s for writing\n", gen_J_path);
                    }
                    }

                    /* ===== 次世代の作成（エリート選択 + 突然変異）===== */
                    int w = 0; // next配列のインデックス（作成済みの子個体数）

                    // エリート選択: 上位elit_n個体を親にして子を作る
                    // 各親から複数の子を作り、合計pop個の子を均等に割り振る
                    for (int e = 0; e < elit_n; ++e)
                    {
                        // 残りの親の数と残りの子の枠から、この親が作るべき子の数を計算
                        int remain_groups = elit_n - e;                                // まだ子を作っていない親の数
                        int remain_slots = prm.pop - w;                                // まだ埋まっていない子の枠
                        int reps = (remain_slots + remain_groups - 1) / remain_groups; /* 均等割り（切り上げ）*/

                        // この親(pop[e])からreps個の子を作る
                        for (int r = 0; r < reps && w < prm.pop; ++r)
                        {
                            next[w] = clone_individual(&pop[e]); // e番目の親（上位個体）をコピー

                            // 確率muで突然変異を起こす（基本的にはmu=1で設定）
                            if (urand(&master) < prm.mu)
                            {
                                mutate_one_edge(&next[w], &master); // 1箇所突然変異
                            }
                            ++w;
                        }
                    }

                    // 現世代(pop)を解放して、次世代(next)と入れ替え
                    for (int i = 0; i < prm.pop; ++i)
                    {
                        free_individual(&pop[i]); // 古い個体を削除
                        pop[i] = next[i];         // 新しい個体に置き換え
                    }
                }
            }
#pragma omp barrier
        } /* gens */
    } /* parallel */

    if (best_J_fp)
    {
        fclose(best_J_fp);
    }

    if (worst_J_fp)
    {
        fclose(worst_J_fp);
    }

    for (int i = 0; i < prm.pop; ++i)
        free_individual(&pop[i]);
    free(pop);
    free(next);
    free(outputs);
    for (int t = 0; t < nthreads; ++t)
    {
        free(tb_pool[t].x);
        free(tb_pool[t].x_next);
    }
    free(tb_pool);
    free(thread_rng);
}