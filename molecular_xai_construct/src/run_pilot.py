import json
import os
import sys
from pathlib import Path

import dill
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, rankdata


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "source_data" / "molucn_rawtest"
RELEASED = ROOT / "source_data" / "molucn" / "results"
OUT = ROOT / "results"
METHODS = ["cam", "diff", "gradcam", "gradinput", "ig", "shap"]
RAW_METHOD = {"cam": "cam", "diff": "diff", "gradcam": "gradcam", "gradinput": "gradinput", "ig": "ig", "shap": "graphsvx"}
LOSSES = ["MSE", "MSE+AC", "MSE+UCN"]
METRICS = ["acc_test", "f1_test", "global_dir_test", "local_dir_test"]
THRESHOLDS = list(range(50, 100, 5))
TRANSFORMS = ["raw", "l1", "signed_rank"]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def load_pickle(path):
    with path.open("rb") as handle:
        return dill.load(handle)


def transform(values, mode):
    values = np.asarray(values, dtype=float)
    if mode == "raw":
        return values
    if mode == "l1":
        denom = np.abs(values).sum()
        return values if denom == 0 else values / denom
    ranks = rankdata(np.abs(values), method="average") / max(len(values), 1)
    out = np.sign(values) * ranks
    out[values == 0] = 0
    return out


def global_scores(dataset, colors, mode):
    scores = []
    for pair, predicted in zip(dataset, colors):
        data_i, data_j = pair["data_i"], pair["data_j"]
        pred_i, pred_j = transform(predicted[0], mode), transform(predicted[1], mode)
        mask_i, mask_j = np.asarray(data_i.mask), np.asarray(data_j.mask)
        idx_i, idx_j = np.flatnonzero(mask_i != 0), np.flatnonzero(mask_j != 0)
        mean_i = pred_i[idx_i].mean() if len(idx_i) else 0.0
        mean_j = pred_j[idx_j].mean() if len(idx_j) else 0.0
        diff = float(np.asarray(data_i.y).ravel()[0] - np.asarray(data_j.y).ravel()[0])
        higher = mean_i > mean_j
        scores.append(float((diff > 0 and higher) or (diff < 0 and not higher)))
    return np.asarray(scores)


def target_accuracy(dataset, colors, threshold_index=0):
    values = []
    for pair, predicted in zip(dataset, colors):
        if not bool(pair["mcs"][threshold_index]):
            continue
        for key, pred in zip(["data_i", "data_j"], predicted):
            mask = np.asarray(pair[key].mask)
            idx = np.flatnonzero(mask != 0)
            if len(idx):
                values.append(np.mean(mask[idx] == np.sign(np.asarray(pred)[idx])))
    return float(np.mean(values)) if values else np.nan


def random_accuracy(dataset, seed, threshold_index=0):
    rng = np.random.default_rng(seed)
    colors = []
    for pair in dataset:
        colors.append([rng.normal(size=len(pair["data_i"].mask)), rng.normal(size=len(pair["data_j"].mask))])
    return target_accuracy(dataset, colors, threshold_index)


def baseline_reproduction():
    sys.path.insert(0, str(ROOT / "reference" / "molucn"))
    from molucn.evaluation.explain_color import get_scores
    from molucn.evaluation.explain_direction_global import get_global_directions
    from molucn.evaluation.explain_direction_local import get_local_directions

    target, loss, method = "1D3G-BRE", "MSE", "diff"
    data = load_pickle(RAW / "data" / target / f"{target}_seed_1337_test.pt")
    colors = load_pickle(RAW / "colors" / method / target / f"{target}_seed_1337_nn_{loss}_mean_1.0_{method}_test.pt")
    np.random.seed(1337)
    acc, f1 = get_scores(data, colors, set="test")
    global_dir = get_global_directions(data, colors, set="test")
    local_dir = get_local_directions(data, colors, set="test")
    calculated = pd.DataFrame(
        {
            "mcs": THRESHOLDS,
            "acc_test": np.nanmean(acc, axis=0),
            "f1_test": np.nanmean(f1, axis=0),
            "global_dir_test": np.nanmean(global_dir, axis=0),
            "local_dir_test": np.nanmean(local_dir, axis=0),
        }
    )
    published = pd.read_csv(RELEASED / "mcs_attr_scores_350_stats_lin2.csv")
    published = published[(published.target == target) & (published.loss == loss) & (published.explainer == method)]
    rows = []
    for metric in METRICS:
        merged = calculated[["mcs", metric]].merge(published[["mcs", metric]], on="mcs", suffixes=("_calculated", "_published"))
        for row in merged.itertuples(index=False):
            diff = abs(row[1] - row[2])
            rows.append(
                {
                    "target": target,
                    "loss": loss,
                    "explainer": method,
                    "mcs": int(row.mcs),
                    "metric": metric,
                    "calculated": float(row[1]),
                    "published": float(row[2]),
                    "absolute_difference": float(diff),
                    "within_0_02": bool(diff <= 0.02),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "baseline_reproduction.csv", index=False)
    return result


def bootstrap_tau(matrix_a, matrix_b, n_boot=1000):
    point = kendalltau(rankdata(-np.median(matrix_a, axis=0)), rankdata(-np.median(matrix_b, axis=0))).statistic
    rng = np.random.default_rng(20260809)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, matrix_a.shape[0], matrix_a.shape[0])
        value = kendalltau(rankdata(-np.median(matrix_a[idx], axis=0)), rankdata(-np.median(matrix_b[idx], axis=0))).statistic
        if np.isfinite(value):
            values.append(value)
    return float(point), float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def fixed_effect_decomposition(cube):
    grand = cube.mean()
    n_target, n_method, n_definition = cube.shape
    ss_total = np.square(cube - grand).sum()
    target_mean = cube.mean(axis=(1, 2))
    method_mean = cube.mean(axis=(0, 2))
    definition_mean = cube.mean(axis=(0, 1))
    method_definition_mean = cube.mean(axis=0)
    ss_target = n_method * n_definition * np.square(target_mean - grand).sum()
    ss_method = n_target * n_definition * np.square(method_mean - grand).sum()
    ss_definition = n_target * n_method * np.square(definition_mean - grand).sum()
    interaction = method_definition_mean - method_mean[:, None] - definition_mean[None, :] + grand
    ss_interaction = n_target * np.square(interaction).sum()
    ss_residual = max(0.0, ss_total - ss_target - ss_method - ss_definition - ss_interaction)
    return {
        "ss_total": float(ss_total),
        "eta2_target": float(ss_target / ss_total),
        "eta2_method": float(ss_method / ss_total),
        "eta2_definition": float(ss_definition / ss_total),
        "eta2_method_definition": float(ss_interaction / ss_total),
        "eta2_residual": float(ss_residual / ss_total),
        "marginal_method_variance": float(np.var(method_mean)),
        "marginal_definition_variance": float(np.var(definition_mean)),
        "definition_method_variance_ratio": float(np.var(definition_mean) / np.var(method_mean)) if np.var(method_mean) else np.inf,
    }


def released_analysis():
    sources = {
        "lin1": RELEASED / "mcs_attr_scores_350_stats_lin1.csv",
        "lin2": RELEASED / "mcs_attr_scores_350_stats_lin2.csv",
    }
    ranking_rows, stability_rows, variance_rows = [], [], []
    for source, path in sources.items():
        frame = pd.read_csv(path)
        frame = frame[frame.explainer.isin(METHODS)]
        for loss in LOSSES:
            part_loss = frame[frame.loss == loss]
            for metric in METRICS:
                for threshold in THRESHOLDS:
                    scores = part_loss[part_loss.mcs == threshold].groupby("explainer")[metric].median().reindex(METHODS)
                    ranks = rankdata(-scores.to_numpy(), method="average")
                    for method, score, rank in zip(METHODS, scores, ranks):
                        ranking_rows.append({"source": source, "loss": loss, "metric": metric, "mcs": threshold, "explainer": method, "target_median": score, "rank": rank})
                base = part_loss[part_loss.mcs == 50].pivot(index="target", columns="explainer", values=metric).reindex(columns=METHODS)
                for threshold in THRESHOLDS[1:]:
                    other = part_loss[part_loss.mcs == threshold].pivot(index="target", columns="explainer", values=metric).reindex(columns=METHODS)
                    common = base.index.intersection(other.index)
                    a, b = base.loc[common].to_numpy(), other.loc[common].to_numpy()
                    valid = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
                    a, b = a[valid], b[valid]
                    if a.shape[0] and a.shape[1] >= 4:
                        tau, low, high = bootstrap_tau(a, b)
                        stability_rows.append({"source": source, "loss": loss, "metric": metric, "mcs_a": 50, "mcs_b": threshold, "n_targets": a.shape[0], "n_methods": a.shape[1], "kendall_tau": tau, "bootstrap_ci_low": low, "bootstrap_ci_high": high})
                pivot = part_loss.pivot(index="target", columns=["explainer", "mcs"], values=metric)
                columns = pd.MultiIndex.from_product([METHODS, THRESHOLDS], names=["explainer", "mcs"])
                pivot = pivot.reindex(columns=columns).dropna()
                if len(pivot):
                    cube = pivot.to_numpy().reshape(len(pivot), len(METHODS), len(THRESHOLDS))
                    row = {"source": source, "loss": loss, "metric": metric, "n_complete_targets": len(pivot)}
                    row.update(fixed_effect_decomposition(cube))
                    variance_rows.append(row)
    pd.DataFrame(ranking_rows).to_csv(OUT / "method_rankings.csv", index=False)
    stability = pd.DataFrame(stability_rows)
    variance = pd.DataFrame(variance_rows)
    stability.to_csv(OUT / "rank_stability.csv", index=False)
    variance.to_csv(OUT / "variance_decomposition.csv", index=False)
    return stability, variance


def raw_analysis():
    normalization_rows, flip_rows, real_rows, random_rows, missing = [], [], [], [], []
    data_files = sorted((RAW / "data").glob("*/*_test.pt"))
    total_pairs, total_pairs_mcs50 = 0, 0
    for index, data_path in enumerate(data_files, start=1):
        target = data_path.parent.name
        dataset = load_pickle(data_path)
        mcs = np.asarray([pair["mcs"] for pair in dataset], dtype=bool)
        total_pairs += len(dataset)
        total_pairs_mcs50 += int(mcs[:, 0].sum())
        for seed in [20260809, 20260810, 20260811]:
            random_rows.append({"target": target, "seed": seed, "acc_test_mcs50": random_accuracy(dataset, seed)})
        for method in METHODS:
            raw_method = RAW_METHOD[method]
            for loss in LOSSES:
                color_path = RAW / "colors" / raw_method / target / f"{target}_seed_1337_nn_{loss}_mean_1.0_{raw_method}_test.pt"
                if not color_path.exists():
                    missing.append({"target": target, "loss": loss, "explainer": method})
                    continue
                colors = load_pickle(color_path)
                if len(colors) != len(dataset):
                    raise RuntimeError(f"Length mismatch: {color_path}")
                if loss == "MSE":
                    real_rows.append({"target": target, "explainer": method, "acc_test_mcs50": target_accuracy(dataset, colors)})
                mode_scores = {}
                for mode in TRANSFORMS:
                    scores = global_scores(dataset, colors, mode)
                    mode_scores[mode] = scores
                    for idx, threshold in enumerate(THRESHOLDS):
                        eligible = mcs[:, idx]
                        normalization_rows.append({"target": target, "loss": loss, "explainer": method, "normalization": mode, "mcs": threshold, "global_dir_test": float(scores[eligible].mean()) if eligible.any() else np.nan, "n_pairs": int(eligible.sum())})
                eligible = mcs[:, 0]
                stacked = np.vstack([mode_scores[mode][eligible] for mode in TRANSFORMS]).T
                flipped = np.any(stacked != stacked[:, [0]], axis=1)
                flip_rows.append({"target": target, "loss": loss, "explainer": method, "n_pairs": int(eligible.sum()), "n_flipped": int(flipped.sum()), "flip_rate": float(flipped.mean()) if len(flipped) else np.nan})
        if index % 25 == 0 or index == len(data_files):
            print(f"RAW_PROGRESS {index}/{len(data_files)}", flush=True)
    normalization = pd.DataFrame(normalization_rows)
    flips = pd.DataFrame(flip_rows)
    real = pd.DataFrame(real_rows)
    random = pd.DataFrame(random_rows)
    normalization.to_csv(OUT / "normalization_sensitivity.csv", index=False)
    flips.to_csv(OUT / "verdict_flips.csv", index=False)
    real.to_csv(OUT / "real_accuracy_mcs50.csv", index=False)
    random.to_csv(OUT / "random_accuracy_mcs50.csv", index=False)
    pd.DataFrame(missing).to_csv(OUT / "missing_color_files.csv", index=False)
    random_mean = random.groupby("target", as_index=False).acc_test_mcs50.mean().rename(columns={"acc_test_mcs50": "random_mean"})
    negative_rows = []
    rng = np.random.default_rng(20260809)
    for method in METHODS:
        merged = real[real.explainer == method].merge(random_mean, on="target")
        diff = (merged.acc_test_mcs50 - merged.random_mean).to_numpy()
        boot = np.asarray([diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(2000)])
        low, high = np.quantile(boot, [0.025, 0.975])
        negative_rows.append({"explainer": method, "n_targets": len(diff), "mean_real_minus_random": float(diff.mean()), "median_real_minus_random": float(np.median(diff)), "bootstrap_ci_low": float(low), "bootstrap_ci_high": float(high), "absolute_difference_ge_0_02": bool(abs(diff.mean()) >= 0.02), "ci_excludes_zero": bool(low > 0 or high < 0)})
    negative = pd.DataFrame(negative_rows)
    negative.to_csv(OUT / "negative_control.csv", index=False)
    return {
        "n_test_targets": len(data_files),
        "n_test_pairs": total_pairs,
        "n_test_pairs_mcs50": total_pairs_mcs50,
        "n_color_files_missing": len(missing),
        "normalization_rows": len(normalization),
    }, flips, negative


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("BASELINE_START", flush=True)
    baseline = baseline_reproduction()
    print("RELEASED_ANALYSIS_START", flush=True)
    stability, variance = released_analysis()
    print("RAW_ANALYSIS_START", flush=True)
    raw_qc, flips, negative = raw_analysis()
    lin2 = pd.read_csv(RELEASED / "mcs_attr_scores_350_stats_lin2.csv")
    result_qc = {
        "released_rows": len(lin2),
        "released_targets": int(lin2.target.nunique()),
        "released_explainers": sorted(lin2.explainer.dropna().unique().tolist()),
        "released_thresholds": sorted(lin2.mcs.dropna().astype(int).unique().tolist()),
        "released_missing_fraction": {metric: float(lin2[metric].isna().mean()) for metric in METRICS},
    }
    result_qc.update(raw_qc)
    write_json(OUT / "qc_summary.json", result_qc)
    baseline_gate = bool(baseline.within_0_02.any())
    data_gate = bool(raw_qc["n_test_pairs_mcs50"] >= 1000 and raw_qc["n_test_targets"] >= 100 and len(METHODS) >= 4 and len(THRESHOLDS) >= 4)
    negative_gate = bool(((negative.absolute_difference_ge_0_02) & (negative.ci_excludes_zero)).any())
    primary_stability = stability[(stability.source == "lin2") & (stability.loss == "MSE") & (stability.metric == "acc_test")]
    rank_signal = bool((primary_stability.kendall_tau < 0.80).any())
    primary_variance = variance[(variance.source == "lin2") & (variance.loss == "MSE") & (variance.metric == "acc_test")]
    variance_signal = bool((primary_variance.definition_method_variance_ratio >= 1.0).any())
    primary_flips = flips[flips.loss == "MSE"]
    flip_rate = float(primary_flips.n_flipped.sum() / primary_flips.n_pairs.sum())
    flip_signal = bool(flip_rate >= 0.20)
    construct_gate = rank_signal or variance_signal or flip_signal
    gates = {
        "baseline_reproduction": baseline_gate,
        "data_completeness": data_gate,
        "negative_control": negative_gate,
        "construct_sensitivity": construct_gate,
        "non_redundancy": None,
    }
    preliminary = "INVALID" if not baseline_gate or not data_gate else ("GO_PENDING_NOVELTY_REVIEW" if negative_gate and construct_gate else "STOP")
    decision = {
        "run_id": "jc3_construct_validity_pilot_20260809",
        "preliminary_decision": preliminary,
        "gates": gates,
        "construct_signals": {
            "primary_rank_tau_below_0_80": rank_signal,
            "primary_definition_method_variance_ratio_ge_1": variance_signal,
            "mse_normalization_flip_rate": flip_rate,
            "mse_normalization_flip_rate_ge_0_20": flip_signal,
        },
        "baseline_status": "verified" if baseline_gate else "failed",
        "verification_plan": "independent verification and novelty review",
    }
    write_json(OUT / "go_stop.json", decision)
    print(json.dumps(decision, indent=2), flush=True)


if __name__ == "__main__":
    main()
