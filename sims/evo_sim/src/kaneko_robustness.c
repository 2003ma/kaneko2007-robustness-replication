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
  Random Number Generation (XorShift32 + Box-Muller)
*============================*/
typedef struct
{
    uint32_t state; // XorShift32 internal state (seed for next random number generation)
    int has_spare;  // Whether the second random number from Box-Muller is stored (1=yes, 0=no)
    double spare;   // Second normal random number from Box-Muller (saved for next use)
} RNG;

/**
 * XorShift32: Fast pseudo-random number generation algorithm
 * Generates 32-bit integer random numbers using only bitwise operations (XOR, shift)
 */

static inline uint32_t xorshift32(RNG *r)
{
    uint32_t x = r->state;
    x ^= x << 13; // XOR with left 13-bit shift
    x ^= x >> 17; // XOR with right 17-bit shift
    x ^= x << 5;  // XOR with left 5-bit shift
    r->state = x; // Update state for next generation
    return x;     // Return integer in range [0, 2^32-1]
}

/**
 * urand: Generate uniform random number in range [0, 1)
 * Normalize XorShift32 output (0 to 2^32-1) to range [0, 1)
 */
static inline double urand(RNG *r) { return (xorshift32(r) + 1.0) / 4294967297.0; }

/**
 * nrand: Generate random number following standard normal distribution N(0,1)
 * Box-Muller method: Generate two normal random numbers from two uniform random numbers
 * - Since two are generated per calculation, one is saved for next use for efficiency
 */
static inline double nrand(RNG *r)
{
    // Return the saved second random number from previous generation (fast path)
    if (r->has_spare)
    {
        r->has_spare = 0;
        return r->spare;
    }

    // Generate two normal random numbers using Box-Muller method
    double u1 = urand(r), u2 = urand(r);       // Get two uniform random numbers
    double rad = sqrt(-2.0 * log(u1 + 1e-16)); // Polar radius (1e-16 avoids log(0))
    double z0 = rad * cos(2.0 * M_PI * u2);    // First normal random number (return this time)
    r->spare = rad * sin(2.0 * M_PI * u2);     // Second normal random number (save for next)
    r->has_spare = 1;                          // Record that spare random number is available
    return z0;
}

/**
 * rng_seed: Initialize random number generator
 * @param r RNG structure to initialize
 * @param seed Seed value (uses default if 0)
 *
 * Reset XorShift32 internal state and Box-Muller spare random number
 */
static inline void rng_seed(RNG *r, uint32_t seed)
{
    if (seed == 0)
        seed = 2463534242u; // Default value to avoid seed=0
    r->state = seed;        // Set XorShift32 initial state
    r->has_spare = 0;       // Clear Box-Muller spare random number
    r->spare = 0.0;
}

/**
 * splitmix32: Transform and distribute seed value
 * @param x Original seed value
 * @return Transformed seed value
 *
 * Used to assign different seeds to each thread in parallel processing.
 * If all threads start with the same seed, they will generate the same
 * random number sequence, so this function transforms the seed to get
 * different sequences per thread.
 */
static inline uint32_t splitmix32(uint32_t x)
{
    x += 0x9E3779B9u;                  // Add constant based on golden ratio
    x = (x ^ (x >> 16)) * 0x85EBCA6Bu; // Bit mixing (affect lower bits with upper bits)
    x = (x ^ (x >> 13)) * 0xC2B2AE35u; // Further bit mixing
    x ^= x >> 16;                      // Final distribution
    return x;
}

/*=============================
  Data Structures: Network Individual
*============================*/

/**
 * RowIdx: Index for fast access of non-zero elements in each row of a matrix
 *
 * Since J matrix only contains {-1, 0, +1}, storing column indices of +1 and -1
 * separately for each row enables fast dot product calculation
 * (only addition and subtraction needed, no multiplication)
 */
typedef struct
{
    int *pos;    // Index array of columns j where J[i][j]=+1
    int *neg;    // Index array of columns j where J[i][j]=-1
    int npos;    // Number of +1 elements (actual usage of pos)
    int nneg;    // Number of -1 elements (actual usage of neg)
    int cap_pos; // Capacity of pos array (for dynamic expansion)
    int cap_neg; // Capacity of neg array (for dynamic expansion)
} RowIdx;

/**
 * Individual: Evolving individual (GRN network)
 *
 * Each individual has a gene regulatory network (J matrix) and fitness value
 */
typedef struct
{
    int N;          // Number of genes (network nodes)
    int8_t *J;      // N×N interaction matrix (only values -1, 0, +1)
    RowIdx *rows;   // Fast index for N rows (positions of +1, -1 elements in each row)
    double fitness; // Fitness value (0 is best, negative values are worse)
    double V_ip;    // Isogenic variance proxy
} Individual;

// Macro to easily access J(i,j)
#define JAT(ind, i, j) ((ind)->J[(size_t)(i) * (ind)->N + (size_t)(j)])

/*=============================
  Utility Functions
*/
// Return a value from {-1, 0, 1} different from current value
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

// Free memory
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
 * build_row_index: Build index for fast access from J matrix
 * @param ind Individual (containing J matrix and rows)
 *
 * [Purpose] Speed up dot product calculation
 * For each row of J matrix, record which columns have +1 and -1 elements.
 * This allows computing h = Σ J[i][j] * x[j] by simply
 * adding x values at +1 columns and subtracting at -1 columns (no multiplication)
 *
 * [Processing flow]
 * 1. For each row i, count the number of +1 and -1 elements
 * 2. Allocate arrays for those counts
 * 3. Scan again and store column numbers of +1/-1 elements in arrays
 */
static void build_row_index(Individual *ind)
{
    int N = ind->N;
    // Allocate RowIdx structures for N rows
    ind->rows = (RowIdx *)calloc((size_t)N, sizeof(RowIdx));
    if (!ind->rows)
    {
        fprintf(stderr, "[FATAL] rows alloc failed\n");
        exit(1);
    }

    for (int i = 0; i < N; ++i)
    {
        // ===== Step 1: Count +1 and -1 elements in each row =====
        int npos = 0, nneg = 0;
        for (int j = 0; j < N; ++j)
        {
            int8_t v = JAT(ind, i, j);
            npos += (v == 1);  // Increment if +1
            nneg += (v == -1); // Increment if -1
        }

        // ===== Step 2: Allocate arrays based on element counts =====
        RowIdx *ri = &ind->rows[i];
        ri->cap_pos = ri->npos = npos; // Set count and capacity for +1
        ri->cap_neg = ri->nneg = nneg; // Set count and capacity for -1
        if (npos > 0)
            ri->pos = (int *)malloc(sizeof(int) * npos); // Allocate array for +1
        if (nneg > 0)
            ri->neg = (int *)malloc(sizeof(int) * nneg); // Allocate array for -1
        if ((npos > 0 && !ri->pos) || (nneg > 0 && !ri->neg))
        {
            fprintf(stderr, "[FATAL] rowidx alloc failed\n");
            exit(1);
        }

        // ===== Step 3: Scan again and store column indices in arrays =====
        int ap = 0, an = 0; // Indices for pos and neg arrays
        for (int j = 0; j < N; ++j)
        {
            int8_t v = JAT(ind, i, j);
            if (v == 1)
                ri->pos[ap++] = j; // Record column j if J[i][j]=+1
            else if (v == -1)
                ri->neg[an++] = j; // J[i][j]=-1 なら列番号jを記録
        }
    }
}

/**
 * init_random_individual: Initialize individual randomly
 * @param ind Individual structure
 * @param N Number of genes
 * @param p_edge Probability of interaction existence
 * @param rng Random number generator
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
    // Initialize J
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // Fill interactions with probability p_edge
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
 * clone_individual: Clone individual (deep copy)
 * @param src Source individual for copying
 * @return Newly copied individual
 *
 * Allocate new J matrix and RowIdx index and copy them
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
    /* Rebuild rows (construct from J) */
    dst.rows = NULL;
    build_row_index(&dst);
    return dst;
}

/**
 * mutate_one_edge: Mutate one edge of individual's J matrix
 * @param ind Individual to mutate
 * @param rng Random number generator
 *
 * Randomly select (i,j) and change J[i][j] to one of {-1, 0, +1}
 * Do nothing if new value equals old value
 * Also update RowIdx index
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
        // If cur was -1, find and remove column number j from neg array
        for (int t = 0; t < ri->nneg; ++t)
            if (ri->neg[t] == j)
            {
                // 同様に、最後の要素で上書きして削除
                ri->neg[t] = ri->neg[--ri->nneg];
                break;
            }
    }
    // Do nothing if cur==0 (not already in index)

    /* ===== Step 2: Add new value (nxt) to RowIdx ===== */
    if (nxt == 1)
    {
        // If nxt is +1, add column j to pos array
        // Expand if capacity insufficient (double size, initial 4)
        if (ri->npos == ri->cap_pos)
        {
            ri->cap_pos = ri->cap_pos ? ri->cap_pos * 2 : 4;
            ri->pos = (int *)realloc(ri->pos, sizeof(int) * ri->cap_pos);
        }
        // Append to array end and increment size
        ri->pos[ri->npos++] = j;
    }
    else if (nxt == -1)
    {
        // If nxt is -1, add column j to neg array
        // Similarly check capacity and expand
        if (ri->nneg == ri->cap_neg)
        {
            ri->cap_neg = ri->cap_neg ? ri->cap_neg * 2 : 4;
            ri->neg = (int *)realloc(ri->neg, sizeof(int) * ri->cap_neg);
        }
        // Append to array end
        ri->neg[ri->nneg++] = j;
    }
    // Do nothing if nxt==0 (no need to register in index)
}

// Free individual memory
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
  Numerical Dynamics
===============================*/
#define TANH(x) tanh(x)

// Compute dot product with vector x based on row index
// Accept input only from j>=k (internal genes)
// RowIdx *ri is the index i-th row, use pos, neg of i-th row index
static inline double dot_row_idx(const RowIdx *ri, const double *restrict x, int k)
{
    double h = 0.0;
    // For +1 elements, only add x when j>=k
    for (int t = 0; t < ri->npos; ++t)
    {
        int j = ri->pos[t];
        if (j >= k)
            h += x[j];
    }
    // For -1 elements, only subtract x when j>=k
    for (int t = 0; t < ri->nneg; ++t)
    {
        int j = ri->neg[t];
        if (j >= k)
            h -= x[j];
    }
    return h;
}

/* Evaluate one trial (one iteration within L) */
static inline double simulate_once_and_score_buf(
    const Individual *ind, const int *outputs, int k,
    double beta, double sigma, double dt, int relax_steps, int meas_steps,
    double *restrict x, double *restrict x_next, RNG *rng)
{
    const int N = ind->N;
    const double sqrt_dt = sqrt(dt);
    // Initialize all spins to -1
    for (int i = 0; i < N; ++i)
        x[i] = -1.0;

    // Relaxation steps
    for (int t = 0; t < relax_steps; ++t)
    {
        // Update state of all genes
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
            // Count whether gene is expressed
            on_now += (x[outputs[a]] > 0.0);
        sum_on += (double)on_now;
    }

    double mean_on = sum_on / (double)meas_steps; // Average number of ON output genes
    return mean_on - (double)k;                   // Closer to 0 is better
}

/* Evaluate same individual over L trials
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
  Compare individuals in descending order of fitness
  @param a Pointer to individual 1 for comparison
  @param b Pointer to individual 2 for comparison
*/
static int cmp_ind_desc(const void *a, const void *b)
{
    const Individual *x = (const Individual *)a, *y = (const Individual *)b;
    return (x->fitness > y->fitness) ? -1 : ((x->fitness < y->fitness) ? 1 : 0);
}

/**
 * Simulation parameter structure
 */
typedef struct
{
    int N, k, pop, gens;
    double beta, sigma, dt;
    int relax, meas, L;
    double pedge, elit;
    double mu; // Mutation rate (0.0-1.0, probability of mutating each edge)
    uint32_t seed;
    const char *load_J_path; // CSV file path to load J matrix (NULL for random initialization)
} Params;

/* Directory for saving J (specified externally) */
static const char *g_J_save_dir = NULL;

void set_J_save_dir(const char *dir)
{
    g_J_save_dir = dir;
}

/* Select output genes (simply select first k genes here) */
static void choose_outputs(int N, int k, int *outputs, RNG *rng)
{
    (void)rng;
    for (int a = 0; a < k; ++a)
        outputs[a] = a;
}

/**
 * Extract generation number from file path
 * @param filepath File path (e.g., "results/ver3/sigma_0.020/gen_20_all_J.csv")
 * @return Generation number (0 if not found)
 *
 * Search for "gen_XX_" pattern and extract numeric part
 */
static int extract_generation_from_path(const char *filepath)
{
    if (!filepath)
        return 0;
    
    // Find filename part (after last '/')
    const char *filename = strrchr(filepath, '/');
    if (filename)
        filename++; // From character after '/'
    else
        filename = filepath; // If no path separator, use entire string
    
    // Search for "gen_"
    const char *gen_pos = strstr(filename, "gen_");
    if (!gen_pos)
        return 0;
    
    // Get numeric value after "gen_"
    int gen_num = 0;
    if (sscanf(gen_pos, "gen_%d", &gen_num) == 1)
        return gen_num;
    
    return 0;
}

/**
 * Load population from CSV file
 * @param filepath Path to CSV file (gen_XX_all_J.csv format)
 * @param pop Individual array (must be pre-allocated)
 * @param expected_pop Expected number of individuals
 * @param expected_N Expected number of genes
 * @return 0 on success, -1 on failure
 */
static int load_population_from_csv(const char *filepath, Individual *pop, int expected_pop, int expected_N)
{
    FILE *fp = fopen(filepath, "r");
    if (!fp)
    {
        fprintf(stderr, "[ERROR] Cannot open %s for loading J matrices\n", filepath);
        return -1;
    }

    char line[65536]; // Sufficiently large buffer (for N=64, ~4096 elements)

    // Skip header row (id,fitness,v_ip,J_0_0,J_0_1,...)
    if (!fgets(line, sizeof(line), fp))
    {
        fprintf(stderr, "[ERROR] Cannot read header from %s\n", filepath);
        fclose(fp);
        return -1;
    }

    // Load each individual
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

        // Parse line: id,fitness,v_ip,J_0_0,J_0_1,...
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

        // J matrix (N*N elements)
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

        // Build RowIdx
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

/* Thread-private work buffer */
typedef struct
{
    double *x, *x_next;
    int capN;
} ThreadBuf;

// Allocate and expand thread-private work buffer to match N
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

/* Function to run simulation
   @param prm_in Pointer to simulation parameters
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

    // ===== Initialization =====
    // Initialize master random number generator (for mutations, etc.)
    RNG master;
    rng_seed(&master, prm.seed);

    // Allocate array to store output gene indices
    int *outputs = (int *)malloc((size_t)prm.k * sizeof(int));
    choose_outputs(prm.N, prm.k, outputs, &master); // Designate first k genes as output

    // Allocate population (current and next generation)
    Individual *pop = (Individual *)malloc((size_t)prm.pop * sizeof(Individual));
    Individual *next = (Individual *)malloc((size_t)prm.pop * sizeof(Individual));
    if (!pop || !next)
    {
        fprintf(stderr, "malloc pop failed\n");
        free(outputs);
        return;
    }

    // Set initial population (load from CSV or generate randomly)
    int start_gen = 0; // Starting generation number
    if (prm.load_J_path)
    {
        // Load from CSV file
        if (load_population_from_csv(prm.load_J_path, pop, prm.pop, prm.N) != 0)
        {
            fprintf(stderr, "[ERROR] Failed to load population from %s\n", prm.load_J_path);
            free(pop);
            free(next);
            free(outputs);
            return;
        }
        // Extract generation number from filename
        start_gen = extract_generation_from_path(prm.load_J_path);
        printf("[INFO] Starting evolution from loaded population (file: %s, generation: %d)\n", 
               prm.load_J_path, start_gen);
    }
    else
    {
        // Generate randomly
        for (int i = 0; i < prm.pop; ++i)
            init_random_individual(&pop[i], prm.N, prm.pedge, &master);
        printf("[INFO] Starting evolution from random initial population\n");
    }

    // Calculate elite size (how many top individuals as parents)
    int elit_n = (int)floor(prm.elit * prm.pop);

    setvbuf(stdout, NULL, _IOLBF, 0); /* Line buffering */
    printf("# Kaneko2007-like GRN evolution (C, more_fast_idx)\n");
    printf("# N=%d k=%d pop=%d gens=%d beta=%.3f sigma=%.4f dt=%.3f relax=%d meas=%d L=%d pedge=%.3f elit=%.2f mu=%.3f seed=%u\n",
           prm.N, prm.k, prm.pop, prm.gens, prm.beta, prm.sigma, prm.dt, prm.relax, prm.meas, prm.L, prm.pedge, prm.elit, prm.mu, prm.seed);
    printf("# generation,best,mean,worst,V_g,V_ip\n");

    /* Directory structure: results/verN/sigma_X.XX/ */
    /* Check existing versions to see if matching parameter set exists */
    (void)mkdir("results", 0755);

    char ver_dir[256];
    char params_path[512];
    int ver_num = 1;
    int found_matching_ver = 0;

    /* Check existing versions */
    for (int v = 1; v < 100; ++v)
    {
        snprintf(ver_dir, sizeof(ver_dir), "results/ver%d", v);
        snprintf(params_path, sizeof(params_path), "%s/params.txt", ver_dir);

        FILE *check_fp = fopen(params_path, "r");
        if (check_fp)
        {
            /* Read params.txt and check if parameters match */
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

    /* Data file */
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

    /* File to save best individual J */
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

    /* File to save worst individual J */
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
        uint32_t s = splitmix32(prm.seed ^ (uint32_t)t ^ 0xA5A5A5A5u); // Generate different seed per thread (scattered)
        rng_seed(&thread_rng[t], s);
        tb_pool[t].x = tb_pool[t].x_next = NULL;
        tb_pool[t].capN = 0;
    }

/* ===== Permanent Parallel Team ===== */
#pragma omp parallel // Spawn threads for maximum cores
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        ThreadBuf *tb = &tb_pool[tid];

        for (int g = 0; g < prm.gens; ++g)
        {
            int current_gen = start_gen + g; // Actual generation number

/* Individual evaluation: parallel over population */
#pragma omp for schedule(static, OMP_CHUNK)
            for (int i = 0; i < prm.pop; ++i)
            {
                tb_ensure(tb, pop[i].N);
                RNG rng = thread_rng[tid];
                rng.state ^= splitmix32((uint32_t)(current_gen * 2654435761u + i)); // Further vary rng by generation and individual index

                double vip, fit;
                evaluate_fitness_buf(&pop[i], outputs, prm.k, prm.beta, prm.sigma, prm.dt, // k receives input from j>=k
                                     prm.relax, prm.meas, prm.L, &rng, &vip, &fit,
                                     tb->x, tb->x_next);
                pop[i].V_ip = vip;
                pop[i].fitness = fit;
                thread_rng[tid] = rng;
            }

            /* Thread synchronization: Initialize shared accumulation in single section,
               no empty atomic directive needed (empty statement causes compiler error).
               Move processing to following single section. */

#pragma omp single
            {
                /* Initialize shared accumulation at start of single section */
            }

/* Recompute in following single section instead of using static variables for sharing */
#pragma omp barrier

/* All individual fitness already stored in each element → aggregate and sort in single */
#pragma omp single
            {
                // Calculate Vg
                // Calculate average fitness of population (for statistics output)
                double mean_fit_sum = 0.0, fitness_squared_sum = 0.0;
                for (int i = 0; i < prm.pop; ++i)
                {
                    mean_fit_sum += pop[i].fitness;
                    fitness_squared_sum += pop[i].fitness * pop[i].fitness;
                }
                // Average fitness of all individuals in this generation
                double mean_fit = mean_fit_sum / (double)prm.pop;

                // Sort in descending order of fitness
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

                /* Save best individual J */
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

                /* Save worst individual J */
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
                    /* Save all individual J matrices in final generation (id,fitness,v_ip,J_00,... format) */
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
                    /* Save J matrix every 20 generations */
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