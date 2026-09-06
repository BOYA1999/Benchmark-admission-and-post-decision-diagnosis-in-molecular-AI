import hashlib
import json
import platform
import sys
import warnings
from pathlib import Path

import torch
import torch_geometric
import dill
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.stats import kendalltau, rankdata


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASE = ROOT / "molecular_xai_construct"
RAW = CASE / "source_data" / "molucn_rawtest"
FULL = CASE / "source_data" / "molucn"
RELEASED = FULL / "results"
ARCHIVED_RESULTS = CASE / "results"
REFERENCE = CASE / "reference" / "molucn"

METHODS = ["cam", "diff", "gradcam", "gradinput", "ig", "shap"]
RAW_METHOD = {
    "cam": "cam",
    "diff": "diff",
    "gradcam": "gradcam",
    "gradinput": "gradinput",
    "ig": "ig",
    "shap": "graphsvx",
}
THRESHOLDS = list(range(50, 100, 5))
METRICS = ["acc_test", "f1_test", "global_dir_test", "local_dir_test"]
SENTINEL_TARGET = "1D3G-BRE"
SENTINEL_METHOD = "diff"
SENTINEL_LOSS = "MSE"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def load_pickle(path):
    with path.open("rb") as handle:
        return dill.load(handle)


def transform(values, mode):
    values = np.asarray(values, dtype=float)
    if mode == "raw":
        return values
    if mode == "l1":
        denominator = np.abs(values).sum()
        return values if denominator == 0 else values / denominator
    ranks = rankdata(np.abs(values), method="average") / max(len(values), 1)
    out = np.sign(values) * ranks
    out[values == 0] = 0
    return out


def target_accuracy_by_threshold(dataset, colors):
    sums = np.zeros(len(THRESHOLDS), dtype=float)
    counts = np.zeros(len(THRESHOLDS), dtype=int)
    pair_counts = np.zeros(len(THRESHOLDS), dtype=int)
    for pair, predicted in zip(dataset, colors):
        eligible = np.asarray(pair["mcs"], dtype=bool)
        pair_counts += eligible.astype(int)
        pair_scores = []
        for key, pred in zip(["data_i", "data_j"], predicted):
            true = np.asarray(pair[key].mask).reshape(-1)
            idx = np.flatnonzero(true != 0)
            if not len(idx):
                pair_scores.append(None)
                continue
            pred_noncommon = np.asarray(pred, dtype=float).reshape(-1)[idx].copy()
            zeros = pred_noncommon == 0
            if zeros.any():
                pred_noncommon[zeros] += np.random.uniform(
                    low=-1e-8, high=1e-8, size=int(zeros.sum())
                )
            score = float(np.mean(true[idx] == np.sign(pred_noncommon)))
            pair_scores.append(score)
        for key, pred in zip(["data_i", "data_j"], predicted):
            true = np.asarray(pair[key].mask).reshape(-1)
            idx = np.flatnonzero(true != 0)
            if not len(idx):
                continue
            pred_noncommon = np.asarray(pred, dtype=float).reshape(-1)[idx].copy()
            zeros = pred_noncommon == 0
            if zeros.any():
                np.random.uniform(low=-1e-8, high=1e-8, size=int(zeros.sum()))
        for score in pair_scores:
            if score is None:
                continue
            sums[eligible] += score
            counts[eligible] += 1
    values = np.divide(
        sums,
        counts,
        out=np.full(len(THRESHOLDS), np.nan, dtype=float),
        where=counts != 0,
    )
    return values, counts, pair_counts


def global_verdict(pair, predicted, mode):
    data_i, data_j = pair["data_i"], pair["data_j"]
    pred_i = transform(predicted[0], mode)
    pred_j = transform(predicted[1], mode)
    mask_i = np.asarray(data_i.mask)
    mask_j = np.asarray(data_j.mask)
    idx_i = np.flatnonzero(mask_i != 0)
    idx_j = np.flatnonzero(mask_j != 0)
    mean_i = pred_i[idx_i].mean() if len(idx_i) else 0.0
    mean_j = pred_j[idx_j].mean() if len(idx_j) else 0.0
    activity_difference = float(
        np.asarray(data_i.y).ravel()[0] - np.asarray(data_j.y).ravel()[0]
    )
    higher = mean_i > mean_j
    return float(
        (activity_difference > 0 and higher)
        or (activity_difference < 0 and not higher)
    )


def flip_counts_at_mcs50(dataset, colors):
    n_pairs = 0
    n_flipped = 0
    for pair, predicted in zip(dataset, colors):
        if not bool(pair["mcs"][0]):
            continue
        verdicts = [
            global_verdict(pair, predicted, mode)
            for mode in ["raw", "l1", "signed_rank"]
        ]
        n_pairs += 1
        n_flipped += int(any(value != verdicts[0] for value in verdicts[1:]))
    return n_pairs, n_flipped


def fixed_effect_ratio(cube):
    method_mean = cube.mean(axis=(0, 2))
    definition_mean = cube.mean(axis=(0, 1))
    method_variance = float(np.var(method_mean))
    definition_variance = float(np.var(definition_mean))
    ratio = definition_variance / method_variance if method_variance else np.inf
    return method_variance, definition_variance, float(ratio)


def sentinel_reconciliation():
    sys.path.insert(0, str(REFERENCE))
    from molucn.evaluation.explain_color import get_scores
    from molucn.evaluation.explain_direction_global import get_global_directions
    from molucn.evaluation.explain_direction_local import get_local_directions

    data_path = (
        RAW
        / "data"
        / SENTINEL_TARGET
        / f"{SENTINEL_TARGET}_seed_1337_test.pt"
    )
    color_path = (
        RAW
        / "colors"
        / SENTINEL_METHOD
        / SENTINEL_TARGET
        / f"{SENTINEL_TARGET}_seed_1337_nn_{SENTINEL_LOSS}_mean_1.0_{SENTINEL_METHOD}_test.pt"
    )
    dataset = load_pickle(data_path)
    colors = load_pickle(color_path)
    np.random.seed(1337)
    acc, f1 = get_scores(dataset, colors, set="test")
    global_direction = get_global_directions(dataset, colors, set="test")
    local_direction = get_local_directions(dataset, colors, set="test")
    current = pd.DataFrame(
        {
            "mcs": THRESHOLDS,
            "acc_test": np.nanmean(acc, axis=0),
            "f1_test": np.nanmean(f1, axis=0),
            "global_dir_test": np.nanmean(global_direction, axis=0),
            "local_dir_test": np.nanmean(local_direction, axis=0),
        }
    )

    np.random.seed(1337)
    manual_acc, _, _ = target_accuracy_by_threshold(dataset, colors)
    if not np.allclose(manual_acc, current.acc_test, atol=1e-12, rtol=0):
        raise RuntimeError("Independent accuracy recomputation disagrees with current evaluator.")

    archived = pd.read_csv(RELEASED / "mcs_attr_scores_350_stats_lin2.csv")
    archived = archived[
        (archived.target == SENTINEL_TARGET)
        & (archived.loss == SENTINEL_LOSS)
        & (archived.explainer == SENTINEL_METHOD)
    ]
    rows = []
    for metric in METRICS:
        merged = current[["mcs", metric]].merge(
            archived[["mcs", metric]], on="mcs", suffixes=("_recomputed", "_archived")
        )
        for row in merged.itertuples(index=False):
            recomputed = float(row[1])
            archived_value = float(row[2])
            difference = abs(archived_value - recomputed)
            cause_status = "exact_match" if difference <= 1e-12 else "difference_recorded"
            rows.append(
                {
                    "target": SENTINEL_TARGET,
                    "loss": SENTINEL_LOSS,
                    "explainer": SENTINEL_METHOD,
                    "metric": metric,
                    "threshold": int(row.mcs),
                    "archived_value": archived_value,
                    "recomputed_value": recomputed,
                    "absolute_difference": difference,
                    "within_0_02": bool(difference <= 0.02),
                    "cause_status": cause_status,
                    "baseline_gate_use": "yes_collectively_via_any_cell_with_abs_diff_le_0.02",
                    "construct_gate_use": "no",
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(HERE / "mxc_sentinel_reconciliation.csv", index=False)
    return frame, data_path, color_path


def recompute_raw_gate_inputs():
    metric_rows = []
    flip_rows = []
    missing_rows = []
    data_files = sorted((RAW / "data").glob("*/*_test.pt"))
    for index, data_path in enumerate(data_files, start=1):
        target = data_path.parent.name
        dataset = load_pickle(data_path)
        for method in METHODS:
            raw_method = RAW_METHOD[method]
            color_path = (
                RAW
                / "colors"
                / raw_method
                / target
                / f"{target}_seed_1337_nn_MSE_mean_1.0_{raw_method}_test.pt"
            )
            if not color_path.exists():
                missing_rows.append(
                    {"target": target, "loss": "MSE", "explainer": method}
                )
                continue
            colors = load_pickle(color_path)
            if len(colors) != len(dataset):
                raise RuntimeError(f"Length mismatch: {color_path}")
            np.random.seed(1337)
            accuracy, molecule_counts, pair_counts = target_accuracy_by_threshold(
                dataset, colors
            )
            for threshold, value, n_molecules, n_pairs in zip(
                THRESHOLDS, accuracy, molecule_counts, pair_counts
            ):
                metric_rows.append(
                    {
                        "target": target,
                        "loss": "MSE",
                        "explainer": method,
                        "threshold": threshold,
                        "acc_test": value,
                        "n_molecule_scores": int(n_molecules),
                        "n_pairs": int(n_pairs),
                    }
                )
            n_pairs, n_flipped = flip_counts_at_mcs50(dataset, colors)
            flip_rows.append(
                {
                    "target": target,
                    "loss": "MSE",
                    "explainer": method,
                    "n_pairs": n_pairs,
                    "n_flipped": n_flipped,
                    "flip_rate": n_flipped / n_pairs if n_pairs else np.nan,
                }
            )
        if index % 25 == 0 or index == len(data_files):
            print(f"RAW_GATE_PROGRESS {index}/{len(data_files)}", flush=True)

    metrics = pd.DataFrame(metric_rows)
    flips = pd.DataFrame(flip_rows)
    missing = pd.DataFrame(missing_rows)
    metrics.to_csv(HERE / "mxc_raw_gate_target_metrics.csv", index=False)
    flips.to_csv(HERE / "mxc_raw_flip_counts.csv", index=False)
    missing.to_csv(HERE / "mxc_raw_missing_inputs.csv", index=False)

    archived_flips = pd.read_csv(ARCHIVED_RESULTS / "verdict_flips.csv")
    archived_flips = archived_flips[archived_flips.loss == "MSE"][
        ["target", "explainer", "n_pairs", "n_flipped"]
    ].rename(
        columns={
            "n_pairs": "archived_n_pairs",
            "n_flipped": "archived_n_flipped",
        }
    )
    flip_comparison = flips.merge(
        archived_flips, on=["target", "explainer"], how="outer", indicator=True
    )
    flip_target_method_exact = bool(
        (flip_comparison._merge == "both").all()
        and (
            flip_comparison.n_pairs == flip_comparison.archived_n_pairs
        ).all()
        and (
            flip_comparison.n_flipped == flip_comparison.archived_n_flipped
        ).all()
    )

    archived_cells = pd.read_csv(RELEASED / "mcs_attr_scores_350_stats_lin2.csv")
    archived_cells = archived_cells[
        (archived_cells.loss == "MSE") & archived_cells.explainer.isin(METHODS)
    ][["target", "explainer", "mcs", "acc_test"]].rename(
        columns={"mcs": "threshold", "acc_test": "archived_acc_test"}
    )
    cell_comparison = metrics.merge(
        archived_cells,
        on=["target", "explainer", "threshold"],
        how="outer",
        indicator=True,
    )
    cell_comparison["absolute_difference"] = (
        cell_comparison.acc_test - cell_comparison.archived_acc_test
    ).abs()
    cell_comparison["source_presence"] = cell_comparison.pop("_merge").astype(str)
    cell_comparison.to_csv(HERE / "mxc_raw_vs_archived_gate_cells.csv", index=False)

    base = metrics[metrics.threshold == 50].pivot(
        index="target", columns="explainer", values="acc_test"
    ).reindex(columns=METHODS)
    rank_rows = []
    for threshold in THRESHOLDS[1:]:
        other = metrics[metrics.threshold == threshold].pivot(
            index="target", columns="explainer", values="acc_test"
        ).reindex(columns=METHODS)
        common = base.index.intersection(other.index)
        a = base.loc[common].to_numpy()
        b = other.loc[common].to_numpy()
        valid = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
        a = a[valid]
        b = b[valid]
        tau = kendalltau(
            rankdata(-np.median(a, axis=0)),
            rankdata(-np.median(b, axis=0)),
        ).statistic
        rank_rows.append(
            {
                "mcs_a": 50,
                "mcs_b": threshold,
                "n_targets": int(a.shape[0]),
                "n_methods": len(METHODS),
                "kendall_tau": float(tau),
            }
        )
    ranks = pd.DataFrame(rank_rows)
    ranks.to_csv(HERE / "mxc_raw_rank_stability.csv", index=False)

    pivot = metrics.pivot(
        index="target", columns=["explainer", "threshold"], values="acc_test"
    )
    columns = pd.MultiIndex.from_product(
        [METHODS, THRESHOLDS], names=["explainer", "threshold"]
    )
    pivot = pivot.reindex(columns=columns).dropna()
    cube = pivot.to_numpy().reshape(len(pivot), len(METHODS), len(THRESHOLDS))
    method_variance, definition_variance, ratio = fixed_effect_ratio(cube)

    archived_pivot = archived_cells.pivot(
        index="target", columns=["explainer", "threshold"], values="archived_acc_test"
    ).reindex(columns=columns)
    raw_complete = set(pivot.index)
    archived_complete = set(archived_pivot.dropna().index)
    finite_comparable = np.isfinite(
        cell_comparison[["acc_test", "archived_acc_test"]]
    ).all(axis=1)
    comparable_differences = cell_comparison.loc[
        finite_comparable, "absolute_difference"
    ]

    rank_value = float(ranks.kendall_tau.min())
    flip_pairs = int(flips.n_pairs.sum())
    flip_count = int(flips.n_flipped.sum())
    flip_rate = float(flip_count / flip_pairs)
    return {
        "rank_min_tau": rank_value,
        "rank_n_threshold_comparisons": len(ranks),
        "variance_ratio": ratio,
        "variance_n_complete_targets": len(pivot),
        "marginal_method_variance": method_variance,
        "marginal_definition_variance": definition_variance,
        "flip_rate": flip_rate,
        "flip_n_pairs": flip_pairs,
        "flip_n_flipped": flip_count,
        "flip_target_method_counts_exactly_match_archived": flip_target_method_exact,
        "n_data_targets": len(data_files),
        "n_metric_rows": len(metrics),
        "n_missing_mse_color_files": len(missing),
        "n_finite_comparable_raw_archived_cells": int(finite_comparable.sum()),
        "n_exact_raw_archived_cells": int((comparable_differences <= 1e-12).sum()),
        "n_within_0_02_raw_archived_cells": int(
            (comparable_differences <= 0.02).sum()
        ),
        "median_absolute_raw_archived_cell_difference": float(
            comparable_differences.median()
        ),
        "raw_complete_only_targets": sorted(raw_complete - archived_complete),
        "archived_complete_only_targets": sorted(archived_complete - raw_complete),
    }


def archived_gate_inputs():
    stability = pd.read_csv(ARCHIVED_RESULTS / "rank_stability.csv")
    stability = stability[
        (stability.source == "lin2")
        & (stability.loss == "MSE")
        & (stability.metric == "acc_test")
    ]
    variance = pd.read_csv(ARCHIVED_RESULTS / "variance_decomposition.csv")
    variance = variance[
        (variance.source == "lin2")
        & (variance.loss == "MSE")
        & (variance.metric == "acc_test")
    ].iloc[0]
    flips = pd.read_csv(ARCHIVED_RESULTS / "verdict_flips.csv")
    flips = flips[flips.loss == "MSE"]
    return {
        "rank_min_tau": float(stability.kendall_tau.min()),
        "variance_ratio": float(variance.definition_method_variance_ratio),
        "variance_n_complete_targets": int(variance.n_complete_targets),
        "flip_rate": float(flips.n_flipped.sum() / flips.n_pairs.sum()),
        "flip_n_pairs": int(flips.n_pairs.sum()),
        "flip_n_flipped": int(flips.n_flipped.sum()),
    }


def gate_reconciliation(raw, archived):
    rows = [
        {
            "gate_input": "minimum_Kendall_tau_MCS50_vs_55_95",
            "criterion_for_signal": "<0.80",
            "archived_value": archived["rank_min_tau"],
            "raw_recomputed_value": raw["rank_min_tau"],
            "absolute_difference": abs(
                archived["rank_min_tau"] - raw["rank_min_tau"]
            ),
            "raw_signal": bool(raw["rank_min_tau"] < 0.80),
            "raw_denominator_or_scope": (
                f"{raw['rank_n_threshold_comparisons']} threshold comparisons; "
                "pairwise-complete targets reported in mxc_raw_rank_stability.csv"
            ),
            "cause_status": "difference_recorded",
        },
        {
            "gate_input": "definition_to_method_marginal_variance_ratio",
            "criterion_for_signal": ">=1.0",
            "archived_value": archived["variance_ratio"],
            "raw_recomputed_value": raw["variance_ratio"],
            "absolute_difference": abs(
                archived["variance_ratio"] - raw["variance_ratio"]
            ),
            "raw_signal": bool(raw["variance_ratio"] >= 1.0),
            "raw_denominator_or_scope": (
                f"{raw['variance_n_complete_targets']} complete targets x 6 methods x 10 thresholds"
            ),
            "cause_status": "difference_recorded",
        },
        {
            "gate_input": "pooled_MSE_normalization_verdict_flip_rate",
            "criterion_for_signal": ">=0.20",
            "archived_value": archived["flip_rate"],
            "raw_recomputed_value": raw["flip_rate"],
            "absolute_difference": abs(archived["flip_rate"] - raw["flip_rate"]),
            "raw_signal": bool(raw["flip_rate"] >= 0.20),
            "raw_denominator_or_scope": (
                f"{raw['flip_n_flipped']} flipped / {raw['flip_n_pairs']} eligible pair-method rows"
            ),
            "cause_status": (
                "established_exact_reproduction"
                if raw["flip_target_method_counts_exactly_match_archived"]
                else "difference_detected"
            ),
        },
    ]
    frame = pd.DataFrame(rows)
    frame.to_csv(HERE / "mxc_gate_reconciliation.csv", index=False)
    return frame


def input_hashes(sentinel_data, sentinel_color):
    full_data = (
        FULL
        / "data"
        / SENTINEL_TARGET
        / f"{SENTINEL_TARGET}_seed_1337_test.pt"
    )
    full_color = (
        FULL
        / "colors"
        / SENTINEL_METHOD
        / SENTINEL_TARGET
        / f"{SENTINEL_TARGET}_seed_1337_nn_{SENTINEL_LOSS}_mean_1.0_{SENTINEL_METHOD}_test.pt"
    )
    paths = [
        ("sentinel_data_rawtest", sentinel_data),
        ("sentinel_data_full_archive", full_data),
        ("sentinel_color_rawtest", sentinel_color),
        ("sentinel_color_full_archive", full_color),
        ("archived_statistics", RELEASED / "mcs_attr_scores_350_stats_lin2.csv"),
        ("current_explain_color_code", REFERENCE / "molucn" / "evaluation" / "explain_color.py"),
        (
            "current_global_direction_code",
            REFERENCE / "molucn" / "evaluation" / "explain_direction_global.py",
        ),
        (
            "current_local_direction_code",
            REFERENCE / "molucn" / "evaluation" / "explain_direction_local.py",
        ),
        ("frozen_pilot_code", CASE / "src" / "run_pilot.py"),
    ]
    rows = []
    for role, path in paths:
        rows.append(
            {
                "role": role,
                "path_relative_to_project_root": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(HERE / "mxc_reconciliation_input_hashes.csv", index=False)
    if rows[0]["sha256"] != rows[1]["sha256"]:
        raise RuntimeError("Rawtest and full-archive sentinel data copies differ.")
    if rows[2]["sha256"] != rows[3]["sha256"]:
        raise RuntimeError("Rawtest and full-archive sentinel color copies differ.")
    return frame


def write_report(sentinel, gates, raw, archived, hashes):
    exact = int((sentinel.absolute_difference <= 1e-12).sum())
    within = int(sentinel.within_0_02.sum())
    exact_cells = ", ".join(
        f"{row.metric}@{row.threshold}"
        for row in sentinel[sentinel.absolute_difference <= 1e-12].itertuples(
            index=False
        )
    )
    within_nonzero_cells = ", ".join(
        f"{row.metric}@{row.threshold}"
        for row in sentinel[
            (sentinel.absolute_difference > 1e-12) & sentinel.within_0_02
        ].itertuples(index=False)
    )
    above_tolerance_cells = ", ".join(
        f"{row.metric}@{row.threshold}"
        for row in sentinel[~sentinel.within_0_02].itertuples(index=False)
    )
    raw_construct_signal = bool(gates.raw_signal.any())
    decision = "REVIEW_REQUIRED" if raw_construct_signal else "STOP"

    gate_lines = []
    for row in gates.itertuples(index=False):
        gate_lines.append(
            f"| {row.gate_input} | {row.criterion_for_signal} | "
            f"{row.archived_value:.12g} | {row.raw_recomputed_value:.12g} | "
            f"{row.absolute_difference:.12g} | {'是' if row.raw_signal else '否'} | "
            f"{row.raw_denominator_or_scope} |"
        )

    sentinel_lines = []
    for row in sentinel.itertuples(index=False):
        sentinel_lines.append(
            f"| {row.metric} | {row.threshold} | {row.archived_value:.12g} | "
            f"{row.recomputed_value:.12g} | {row.absolute_difference:.12g} | "
            f"{'是' if row.within_0_02 else '否'} | {row.cause_status} | "
            "基线 gate：集体进入；construct gate：不进入 |"
        )

    text = f"""# MXC sentinel reconciliation

## Summary

- The 40 sentinel cells include {exact} exact matches and {within} matches within an absolute tolerance of 0.02.
- All three construct-gate inputs were recomputed from the saved MSE attribution rows. The recorded thresholds yield **{decision}**.

## Sentinel differences

| Class | Count | Cells |
|---|---:|---|
| Exact match | {exact} | {exact_cells} |
| Nonzero difference ≤0.02 | {within - exact} | {within_nonzero_cells} |
| Difference >0.02 | {len(sentinel) - within} | {above_tolerance_cells} |

## Recomputed gate inputs

| gate 输入 | signal 阈值 | archived | raw-recomputed | 绝对差 | raw signal | 分母/范围 |
|---|---:|---:|---:|---:|---|---|
{chr(10).join(gate_lines)}

The ranking input is target-median method rank based on `MSE / acc_test` across MCS 50 and 55–95. The variance input is the marginal definition-to-method variance ratio across the target × method × threshold cube. The flip input is the pooled pair-method proportion whose global-direction verdict changes after L1 or signed-rank normalization.

The archived values are minimum tau {archived['rank_min_tau']:.12g}, variance ratio {archived['variance_ratio']:.12g}, and flip rate {archived['flip_n_flipped']}/{archived['flip_n_pairs']} = {archived['flip_rate']:.12g}. The recomputed marginal method variance is {raw['marginal_method_variance']:.12g} and the marginal definition variance is {raw['marginal_definition_variance']:.12g}.

Across {raw['n_finite_comparable_raw_archived_cells']} comparable primary-accuracy cells, {raw['n_exact_raw_archived_cells']} match exactly, {raw['n_within_0_02_raw_archived_cells']} are within 0.02, and the median absolute difference is {raw['median_absolute_raw_archived_cell_difference']:.12g}.

## Sentinel cells

| metric | threshold | archived | recomputed | 绝对差 | ≤0.02 | cause status | gate 使用 |
|---|---:|---:|---:|---:|---|---|---|
{chr(10).join(sentinel_lines)}

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
"""
    (HERE / "mxc_sentinel_reconciliation.md").write_text(text, encoding="utf-8")
    return decision


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    sentinel, sentinel_data, sentinel_color = sentinel_reconciliation()
    raw = recompute_raw_gate_inputs()
    archived = archived_gate_inputs()
    gates = gate_reconciliation(raw, archived)
    hashes = input_hashes(sentinel_data, sentinel_color)
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "sklearn": sklearn.__version__,
        "torch": torch.__version__,
        "torch_geometric": torch_geometric.__version__,
        "dill": dill.__version__,
        "numpy_seed_per_target_method_accuracy": 1337,
    }
    (HERE / "runtime_environment.json").write_text(
        json.dumps(environment, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    decision = write_report(sentinel, gates, raw, archived, hashes)
    summary = {
        "sentinel_cells": len(sentinel),
        "sentinel_exact": int((sentinel.absolute_difference <= 1e-12).sum()),
        "sentinel_within_0_02": int(sentinel.within_0_02.sum()),
        "raw_gate_inputs": raw,
        "archived_gate_inputs": archived,
        "raw_recomputed_decision": decision,
    }
    (HERE / "reconciliation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
