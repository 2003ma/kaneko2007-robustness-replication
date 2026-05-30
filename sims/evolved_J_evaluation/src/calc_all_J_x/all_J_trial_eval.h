#ifndef ALL_J_TRIAL_EVAL_H
#define ALL_J_TRIAL_EVAL_H

/* 全Jに対するtrial評価パラメータ */
typedef struct {
    const char *input_dir;   /* gen_XXX_all_J.csv があるディレクトリ */
    int generation;          /* 世代番号 */
    int n_trials;            /* 各Jあたりの試行回数 */
    int t1;                  /* 測定開始時刻 */
    int t2;                  /* 測定終了時刻 */
    double beta;             /* ダイナミクスパラメータ */
    double dt;               /* 時間刻み */
    int k_boundary;          /* 階層境界 */
    int seed;                /* RNG seed */
    double sigma;            /* 変異率σ */
    int start_ind;           /* 個体開始インデックス（デフォルト0） */
    int jemk;                /* ダイナミクスでj>=kを考慮するかどうか*/
} AllJTrialParams;

/* 実行関数 */
int run_all_J_trial_eval(const AllJTrialParams *p);

#endif /* ALL_J_TRIAL_EVAL_H */
