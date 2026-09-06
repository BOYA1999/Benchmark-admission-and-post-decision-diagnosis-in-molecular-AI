# MXC archived-vs-recomputed sentinel reconciliation

## Summary

- The 40 sentinel cells include 6 exact matches and 15 matches within an absolute tolerance of 0.02.
- All three construct-gate inputs were recomputed from the saved MSE attribution rows. The recorded thresholds yield **STOP**.

## Sentinel differences

| Class | Count | Cells |
|---|---:|---|
| Exact match | 6 | global_dir_test@70, global_dir_test@80, global_dir_test@95, local_dir_test@65, local_dir_test@90, local_dir_test@95 |
| Nonzero difference ≤0.02 | 9 | global_dir_test@50, global_dir_test@55, global_dir_test@60, global_dir_test@65, global_dir_test@75, local_dir_test@50, local_dir_test@55, local_dir_test@60, local_dir_test@70 |
| Difference >0.02 | 25 | acc_test@50–95 (10), f1_test@50–95 (10), global_dir_test@85, global_dir_test@90, local_dir_test@75, local_dir_test@80, local_dir_test@85 |

## Recomputed gate inputs

| gate 输入 | signal 阈值 | archived | raw-recomputed | 绝对差 | raw signal | 分母/范围 |
|---|---:|---:|---:|---:|---|---|
| minimum_Kendall_tau_MCS50_vs_55_95 | <0.80 | 0.866666666667 | 1 | 0.133333333333 | 否 | 9 threshold comparisons; pairwise-complete targets reported in mxc_raw_rank_stability.csv |
| definition_to_method_marginal_variance_ratio | >=1.0 | 0.127583276491 | 0.114934948848 | 0.0126483276425 | 否 | 232 complete targets x 6 methods x 10 thresholds |
| pooled_MSE_normalization_verdict_flip_rate | >=0.20 | 0.148895815249 | 0.148895815249 | 0 | 否 | 96402 flipped / 647446 eligible pair-method rows |

The ranking input is target-median method rank based on `MSE / acc_test` across MCS 50 and 55–95. The variance input is the marginal definition-to-method variance ratio across the target × method × threshold cube. The flip input is the pooled pair-method proportion whose global-direction verdict changes after L1 or signed-rank normalization.

The archived values are minimum tau 0.866666666667, variance ratio 0.127583276491, and flip rate 96402/647446 = 0.148895815249. Recomputed marginal method variance is 0.00185201161371 and marginal definition variance is 0.000212860860088.

Across 19854 comparable primary-accuracy cells, 1662 match exactly, 4925 are within 0.02, and the median absolute difference is 0.0570808910362.

## Sentinel cells

| metric | threshold | archived | recomputed | absolute difference | ≤0.02 | status | gate use |
|---|---:|---:|---:|---:|---|---|---|
| acc_test | 50 | 0.598915903218 | 0.639414074298 | 0.0404981710796 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 55 | 0.598915903218 | 0.639414074298 | 0.0404981710796 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 60 | 0.598915903218 | 0.639414074298 | 0.0404981710796 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 65 | 0.585265087422 | 0.644730481025 | 0.0594653936024 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 70 | 0.585 | 0.643306122449 | 0.058306122449 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 75 | 0.482323232323 | 0.573063973064 | 0.0907407407407 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 80 | 0.434615384615 | 0.551282051282 | 0.116666666667 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 85 | 0.405405405405 | 0.630630630631 | 0.225225225225 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 90 | 0.527777777778 | 0.75 | 0.222222222222 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| acc_test | 95 | 0.5 | 0.833333333333 | 0.333333333333 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 50 | 0.584718392942 | 0.528149301119 | 0.0565690918229 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 55 | 0.584718392942 | 0.528149301119 | 0.0565690918229 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 60 | 0.584718392942 | 0.528149301119 | 0.0565690918229 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 65 | 0.569147321939 | 0.534198958535 | 0.0349483634046 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 70 | 0.569916623059 | 0.532063682349 | 0.0378529407101 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 75 | 0.41341991342 | 0.469691139388 | 0.0562712259682 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 80 | 0.359340659341 | 0.477118437118 | 0.117777777778 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 85 | 0.328828828829 | 0.574774774775 | 0.245945945946 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 90 | 0.5 | 0.666666666667 | 0.166666666667 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| f1_test | 95 | 0.5 | 0.833333333333 | 0.333333333333 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 50 | 0.825688073394 | 0.834862385321 | 0.00917431192661 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 55 | 0.825688073394 | 0.834862385321 | 0.00917431192661 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 60 | 0.825688073394 | 0.834862385321 | 0.00917431192661 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 65 | 0.82 | 0.83 | 0.01 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 70 | 0.808988764045 | 0.808988764045 | 0 | 是 | established | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 75 | 0.666666666667 | 0.686274509804 | 0.0196078431373 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 80 | 0.647058823529 | 0.647058823529 | 0 | 是 | established | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 85 | 0.55 | 0.75 | 0.2 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 90 | 0.7 | 0.8 | 0.1 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| global_dir_test | 95 | 0.75 | 0.75 | 0 | 是 | established | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 50 | 0.790960451977 | 0.779661016949 | 0.0112994350282 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 55 | 0.790960451977 | 0.779661016949 | 0.0112994350282 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 60 | 0.790960451977 | 0.779661016949 | 0.0112994350282 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 65 | 0.779874213836 | 0.779874213836 | 0 | 是 | established | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 70 | 0.768939393939 | 0.761363636364 | 0.00757575757576 | 是 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 75 | 0.688405797101 | 0.659420289855 | 0.0289855072464 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 80 | 0.666666666667 | 0.62962962963 | 0.037037037037 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 85 | 0.645833333333 | 0.6875 | 0.0416666666667 | 否 | difference_recorded | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 90 | 0.8 | 0.8 | 0 | 是 | established | 基线 gate：集体进入；construct gate：不进入 |
| local_dir_test | 95 | 1 | 1 | 0 | 是 | established | 基线 gate：集体进入；construct gate：不进入 |

## Outputs

- `mxc_sentinel_reconciliation.csv`：40 个 sentinel 单元。
- `mxc_gate_reconciliation.csv`：三个 gate 输入 archived-vs-raw-recomputed 对照。
- `mxc_raw_gate_target_metrics.csv`：原始 attribution 行重算的 target × method × threshold `acc_test` 审计轨迹。
- `mxc_raw_vs_archived_gate_cells.csv`：primary accuracy 原始行复算与 archived 单元级对照。
- `mxc_raw_rank_stability.csv`：9 个 MCS 对照的 Kendall tau。
- `mxc_raw_flip_counts.csv`：target × method 的 flip 分子/分母。
- `mxc_raw_missing_inputs.csv`：缺失 MSE color 输入。
- `mxc_reconciliation_input_hashes.csv`：关键输入及评分代码 SHA-256。
- `runtime_environment.json`：本次复算运行环境。
- `recompute_mxc_sentinel_reconciliation.py`：可重复执行脚本。
