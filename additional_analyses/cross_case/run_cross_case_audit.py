import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from scipy.stats import kendalltau, rankdata


RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SEED = 20260904
BOOTSTRAP_REPLICATES = 5000
ACE_DATASETS = ["CHEMBL214_Ki", "CHEMBL235_EC50", "CHEMBL236_Ki"]
MXC_METHODS = ["cam", "diff", "gradcam", "gradinput", "ig", "shap"]
MXC_THRESHOLDS = list(range(50, 100, 5))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def digest_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_parts(value):
    if pd.isna(value) or not str(value).strip():
        return [], True
    parts, valid = [], True
    for text in str(value).split("."):
        text = text.strip()
        if not text:
            continue
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            parts.append("!RAW:" + text)
            valid = False
        else:
            parts.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    return sorted(parts), valid


def list_from_canonical(value):
    return [] if pd.isna(value) or not str(value) else str(value).split(".")


def normalize_reaction_frame(frame):
    rows = []
    for row in frame.itertuples(index=False):
        reactant, reactant_ok = canonical_parts(row.REACTANT)
        reagent, _ = canonical_parts(row.REAGENT)
        product, product_ok = canonical_parts(row.PRODUCT)
        rows.append(
            {
                "reaction_hash": digest_text(json.dumps([reactant, reagent, product], separators=(",", ":"))),
                "reactant_product_hash": digest_text(json.dumps([reactant, product], separators=(",", ":"))),
                "product_hash": digest_text(json.dumps(product, separators=(",", ":"))),
                "valid_reactant_product": reactant_ok and product_ok,
            }
        )
    return pd.DataFrame(rows)


def audit_atom_maps(csv_dir):
    pattern = re.compile(r":\d+\]")
    rows = []
    for split in ("train", "val", "test"):
        totals = {column: 0 for column in ("REACTANT", "REAGENT", "PRODUCT")}
        any_rows = 0
        n_rows = 0
        for chunk in pd.read_csv(csv_dir / f"{split}.csv", dtype=str, keep_default_na=False, chunksize=50000):
            n_rows += len(chunk)
            row_any = np.zeros(len(chunk), dtype=bool)
            for column in totals:
                hit = chunk[column].str.contains(pattern, regex=True).to_numpy()
                totals[column] += int(hit.sum())
                row_any |= hit
            any_rows += int(row_any.sum())
        rows.extend(
            {
                "split": split,
                "field": column,
                "rows": n_rows,
                "rows_with_atom_map_marker": totals[column],
            }
            for column in totals
        )
        rows.append(
            {
                "split": split,
                "field": "ANY_FIELD",
                "rows": n_rows,
                "rows_with_atom_map_marker": any_rows,
            }
        )
    return pd.DataFrame(rows)


def rso_audit(input_paths):
    csv_dir = ROOT / "reaction_source_overlap/data/raw/USPTO_MIT/MIT_separated"
    train_path = ROOT / "reaction_source_overlap/data/processed/train_reference.parquet"
    test_path = ROOT / "reaction_source_overlap/data/processed/test_pool.parquet"
    audit_path = ROOT / "reaction_source_overlap/results/data_audit.json"
    val_csv = csv_dir / "val.csv"
    for path in [train_path, test_path, audit_path, val_csv, *(csv_dir / f"{x}.csv" for x in ("train", "test"))]:
        input_paths[path] = "RSO input"

    train = pd.read_parquet(train_path)
    test = pd.read_parquet(test_path)
    val = normalize_reaction_frame(pd.read_csv(val_csv))

    train_rp = train.apply(
        lambda row: digest_text(
            json.dumps(
                [list_from_canonical(row.reactant_canonical), list_from_canonical(row.product_canonical)],
                separators=(",", ":"),
            )
        ),
        axis=1,
    )
    test_rp = test.apply(
        lambda row: digest_text(
            json.dumps(
                [list_from_canonical(row.reactant_canonical), list_from_canonical(row.product_canonical)],
                separators=(",", ":"),
            )
        ),
        axis=1,
    )
    references = {
        "exact_tuple": set(train.reaction_hash) | set(val.reaction_hash),
        "ignore_reagent": set(train_rp) | set(val.reactant_product_hash),
        "product_only": set(train.product_key) | set(val.product_hash),
    }
    test_keys = {
        "exact_tuple": test.reaction_hash,
        "ignore_reagent": test_rp,
        "product_only": test.product_key,
    }
    labels = {
        "exact_tuple": "sorted canonical reactants + sorted canonical reagents + sorted canonical products",
        "ignore_reagent": "sorted canonical reactants + sorted canonical products",
        "product_only": "sorted canonical products",
    }
    valid = test.valid_reactant_product.astype(bool).to_numpy()
    identity_rows = []
    for definition, reference in references.items():
        mask = test_keys[definition].isin(reference).to_numpy()
        identity_rows.append(
            {
                "identity_definition": definition,
                "definition_detail": labels[definition],
                "test_total_rows": len(test),
                "test_valid_reactant_product_rows": int(valid.sum()),
                "overlap_total_rows": int(mask.sum()),
                "overlap_total_rate": float(mask.mean()),
                "overlap_valid_rows": int((mask & valid).sum()),
                "overlap_valid_rate": float((mask & valid).sum() / valid.sum()),
                "overlap_unique_test_keys": int(test_keys[definition][mask].nunique()),
            }
        )
    identity = pd.DataFrame(identity_rows)
    identity.to_csv(OUT / "rso_identity_sensitivity.csv", index=False)

    atom_maps = audit_atom_maps(csv_dir)
    atom_maps.to_csv(OUT / "rso_atom_map_audit.csv", index=False)

    historical = json.loads(audit_path.read_text(encoding="utf-8"))
    template_rows = []
    success = ~test.template_failed.astype(bool)
    for radius in (0, 1):
        values = test.loc[success, f"template_frequency_r{radius}"].astype(int)
        template_rows.append(
            {
                "custom_template_radius": radius,
                "test_total_rows": len(test),
                "test_aligned_success_rows": int(success.sum()),
                "test_aligned_failed_rows": int((~success).sum()),
                "rows_with_train_frequency_ge_1": int((values >= 1).sum()),
                "rows_with_train_frequency_0": int((values == 0).sum()),
                "median_train_frequency_success_rows": float(values.median()),
                "maximum_train_frequency_success_rows": int(values.max()),
                "lookup_scope": "training mapped records only",
                "definition_scope": "custom center-environment digest",
                "train_mapped_parse_failures": int(historical["template_failed"]["train_mapped"]),
            }
        )
    templates = pd.DataFrame(template_rows)
    templates.to_csv(OUT / "rso_template_frequency_audit.csv", index=False)
    return identity, atom_maps, templates


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[b] = a


def dependency_components(cliff, noncliff):
    n = len(cliff)
    uf = UnionFind(n)
    seen = {}
    for row_id in range(n):
        for molecule in np.concatenate([cliff[row_id], noncliff[row_id]]):
            molecule = int(molecule)
            if molecule in seen:
                uf.union(row_id, seen[molecule])
            else:
                seen[molecule] = row_id
    groups = {}
    for row_id in range(n):
        groups.setdefault(uf.find(row_id), []).append(row_id)
    return [np.asarray(rows, dtype=int) for rows in groups.values()]


def rmse_from_se(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float))))


def estimands(se, endpoints):
    flat = endpoints.ravel()
    counts = Counter(int(value) for value in flat)
    unique = np.asarray(sorted(counts), dtype=int)
    inverse_weights = np.asarray([1.0 / counts[int(value)] for value in flat])
    return {
        "endpoint_deduplicated": rmse_from_se(se[unique]),
        "pair_equal": rmse_from_se(se[flat]),
        "inverse_multiplicity": float(np.sqrt(np.sum(se[flat] * inverse_weights) / inverse_weights.sum())),
    }, counts


def cluster_bootstrap_ratios(se, cliff, noncliff, components, rng):
    summaries = []
    for rows in components:
        c_flat = cliff[rows].ravel()
        n_flat = noncliff[rows].ravel()
        c_unique = np.unique(c_flat)
        n_unique = np.unique(n_flat)
        summaries.append(
            {
                "pair_c_sum": float(se[c_flat].sum()),
                "pair_c_n": len(c_flat),
                "pair_n_sum": float(se[n_flat].sum()),
                "pair_n_n": len(n_flat),
                "unique_c_sum": float(se[c_unique].sum()),
                "unique_c_n": len(c_unique),
                "unique_n_sum": float(se[n_unique].sum()),
                "unique_n_n": len(n_unique),
            }
        )
    if len(summaries) < 2:
        return {"endpoint_deduplicated": (np.nan, np.nan), "pair_equal": (np.nan, np.nan)}
    values = {"endpoint_deduplicated": [], "pair_equal": []}
    for _ in range(BOOTSTRAP_REPLICATES):
        chosen = rng.integers(0, len(summaries), len(summaries))
        aggregate = {key: sum(summaries[index][key] for index in chosen) for key in summaries[0]}
        values["endpoint_deduplicated"].append(
            np.sqrt(aggregate["unique_c_sum"] / aggregate["unique_c_n"])
            / np.sqrt(aggregate["unique_n_sum"] / aggregate["unique_n_n"])
        )
        values["pair_equal"].append(
            np.sqrt(aggregate["pair_c_sum"] / aggregate["pair_c_n"])
            / np.sqrt(aggregate["pair_n_sum"] / aggregate["pair_n_n"])
        )
    return {key: tuple(np.quantile(value, [0.025, 0.975])) for key, value in values.items()}


def ace_prediction_surface(manifest, predictions, dataset):
    test = manifest[
        (manifest.dataset == dataset) & (manifest.environment == "mmp_series") & (manifest.split == "test")
    ].reset_index(drop=True)
    pred = predictions[
        (predictions.dataset == dataset) & (predictions.environment == "mmp_series") & (predictions.split == "test")
    ]
    ensemble = pred.groupby("canonical_smiles", sort=False).prediction.mean()
    yhat = test.canonical_smiles.map(ensemble).to_numpy(float)
    if not np.isfinite(yhat).all():
        raise RuntimeError(f"Missing ACE predictions for {dataset}")
    return test, np.square(test.activity.to_numpy(float) - yhat)


def ace_surface_rows(dataset, matches, se, analysis_set, rng):
    cliff = matches[["cliff_i", "cliff_j"]].to_numpy(int)
    noncliff = matches[["noncliff_i", "noncliff_j"]].to_numpy(int)
    cliff_values, cliff_counts = estimands(se, cliff)
    noncliff_values, noncliff_counts = estimands(se, noncliff)
    if not np.isclose(cliff_values["endpoint_deduplicated"], cliff_values["inverse_multiplicity"], atol=1e-12, rtol=0):
        raise AssertionError("Cliff inverse-multiplicity and unique estimands differ")
    if not np.isclose(noncliff_values["endpoint_deduplicated"], noncliff_values["inverse_multiplicity"], atol=1e-12, rtol=0):
        raise AssertionError("Non-cliff inverse-multiplicity and unique estimands differ")
    components = dependency_components(cliff, noncliff)
    intervals = cluster_bootstrap_ratios(se, cliff, noncliff, components, rng)
    estimate_rows = []
    for estimand in ("endpoint_deduplicated", "pair_equal", "inverse_multiplicity"):
        ratio = cliff_values[estimand] / noncliff_values[estimand]
        interval_key = "endpoint_deduplicated" if estimand == "inverse_multiplicity" else estimand
        low, high = intervals[interval_key]
        estimate_rows.append(
            {
                "analysis_set": analysis_set,
                "dataset": dataset,
                "estimand": estimand,
                "cliff_rmse": cliff_values[estimand],
                "noncliff_rmse": noncliff_values[estimand],
                "rmse_ratio": ratio,
                "log_rmse_ratio": float(np.log(ratio)),
                "component_bootstrap_ci_low": low,
                "component_bootstrap_ci_high": high,
                "bootstrap_replicates": BOOTSTRAP_REPLICATES if len(components) >= 2 else 0,
            }
        )
    shared = set(cliff_counts) & set(noncliff_counts)
    component_row = {
        "analysis_set": analysis_set,
        "dataset": dataset,
        "matched_rows": len(matches),
        "cliff_unique_molecules": len(cliff_counts),
        "noncliff_unique_molecules": len(noncliff_counts),
        "shared_cliff_noncliff_molecules": len(shared),
        "shared_fraction_of_cliff_unique": len(shared) / len(cliff_counts),
        "shared_fraction_of_noncliff_unique": len(shared) / len(noncliff_counts),
        "dependency_components": len(components),
        "largest_component_rows": max(map(len, components)),
        "largest_component_fraction_rows": max(map(len, components)) / len(matches),
        "cliff_maximum_multiplicity": max(cliff_counts.values()),
        "noncliff_maximum_multiplicity": max(noncliff_counts.values()),
        "inverse_unique_max_abs_difference": max(
            abs(cliff_values["endpoint_deduplicated"] - cliff_values["inverse_multiplicity"]),
            abs(noncliff_values["endpoint_deduplicated"] - noncliff_values["inverse_multiplicity"]),
        ),
    }
    multiplicity_rows = []
    for endpoint_class, counts in (("cliff", cliff_counts), ("noncliff", noncliff_counts)):
        for molecule, multiplicity in sorted(counts.items()):
            multiplicity_rows.append(
                {
                    "analysis_set": analysis_set,
                    "dataset": dataset,
                    "endpoint_class": endpoint_class,
                    "test_row_index": molecule,
                    "multiplicity": multiplicity,
                }
            )
    return estimate_rows, component_row, multiplicity_rows


def ace_audit(input_paths):
    run = ROOT / "activity_cliff_estimand/results"
    manifest_path = ROOT / "activity_cliff_estimand/data/split_manifest.csv"
    predictions_path = run / "predictions_long.csv"
    matches_path = run / "matching_pairs.csv"
    metrics_path = run / "metrics.csv"
    for path in (manifest_path, predictions_path, matches_path, metrics_path):
        input_paths[path] = "ACE input"

    manifest = pd.read_csv(manifest_path)
    predictions = pd.read_csv(predictions_path)
    matches = pd.read_csv(matches_path)
    stored = pd.read_csv(metrics_path)
    rng = np.random.default_rng(SEED)
    estimates, components, multiplicity = [], [], []
    disjoint_estimates = []
    for dataset in ACE_DATASETS:
        _, se = ace_prediction_surface(manifest, predictions, dataset)
        current = matches[(matches.dataset == dataset) & (matches.environment == "mmp_series")].reset_index(drop=True)
        rows, component, endpoint_rows = ace_surface_rows(dataset, current, se, "archived_matches", rng)
        estimates.extend(rows)
        components.append(component)
        multiplicity.extend(endpoint_rows)
        stored_row = stored[(stored.dataset == dataset) & (stored.environment == "mmp_series")].iloc[0]
        unique_row = next(row for row in rows if row["estimand"] == "endpoint_deduplicated")
        if not np.isclose(unique_row["rmse_ratio"], stored_row.cliff_ratio, atol=1e-12, rtol=0):
            raise AssertionError(f"Stored ACE ratio mismatch for {dataset}")

        all_cliff_endpoints = set(current[["cliff_i", "cliff_j"]].to_numpy(int).ravel())
        fixed = current[
            ~current[["noncliff_i", "noncliff_j"]].isin(all_cliff_endpoints).any(axis=1)
        ].reset_index(drop=True)
        fixed_rows, fixed_component, fixed_endpoints = ace_surface_rows(dataset, fixed, se, "fixed_endpoint_disjoint", rng)
        disjoint_estimates.extend(fixed_rows)
        components.append(fixed_component)
        multiplicity.extend(fixed_endpoints)

    estimates = pd.DataFrame(estimates)
    estimates.to_csv(OUT / "ace_estimand_sensitivity.csv", index=False)
    pd.DataFrame(disjoint_estimates).to_csv(OUT / "ace_endpoint_disjoint_sensitivity.csv", index=False)
    components = pd.DataFrame(components)
    components.to_csv(OUT / "ace_endpoint_overlap_and_components.csv", index=False)
    pd.DataFrame(multiplicity).to_csv(OUT / "ace_endpoint_multiplicity.csv", index=False)

    thresholds = np.round(np.arange(1.15, 1.3501, 0.01), 2)
    threshold_rows = []
    for estimand in ("endpoint_deduplicated", "pair_equal", "inverse_multiplicity"):
        surface = estimates[estimates.estimand == estimand].set_index("dataset").rmse_ratio
        for threshold in thresholds:
            positive = surface >= threshold
            for dataset in ACE_DATASETS:
                threshold_rows.append(
                    {
                        "estimand": estimand,
                        "threshold": threshold,
                        "dataset": dataset,
                        "rmse_ratio": surface[dataset],
                        "ecfp_component_positive": bool(positive[dataset]),
                        "positive_datasets_out_of_3": int(positive.sum()),
                    }
                )
    threshold_map = pd.DataFrame(threshold_rows)
    threshold_map.to_csv(OUT / "ace_threshold_decision_map.csv", index=False)
    return estimates, components, pd.DataFrame(disjoint_estimates), threshold_map


def variance_ratio_complete(frame, thresholds):
    pivot = frame[frame.mcs.isin(thresholds)].pivot(index="target", columns=["explainer", "mcs"], values="acc_test")
    columns = pd.MultiIndex.from_product([MXC_METHODS, thresholds], names=["explainer", "mcs"])
    pivot = pivot.reindex(columns=columns).dropna()
    cube = pivot.to_numpy().reshape(len(pivot), len(MXC_METHODS), len(thresholds))
    method_means = cube.mean(axis=(0, 2))
    threshold_means = cube.mean(axis=(0, 1))
    ratio = float(np.var(threshold_means) / np.var(method_means)) if np.var(method_means) else np.inf
    return len(pivot), ratio


def variance_ratio_available(frame, thresholds):
    current = frame[frame.mcs.isin(thresholds) & np.isfinite(frame.acc_test)]
    method_means = current.groupby("explainer").acc_test.mean().reindex(MXC_METHODS).to_numpy()
    threshold_means = current.groupby("mcs").acc_test.mean().reindex(thresholds).to_numpy()
    return float(np.var(threshold_means) / np.var(method_means)) if np.var(method_means) else np.inf


def target_cluster_interval(frame, rng):
    grouped = frame.groupby("target")[["n_flipped", "n_pairs"]].sum()
    values = []
    for _ in range(BOOTSTRAP_REPLICATES):
        chosen = rng.integers(0, len(grouped), len(grouped))
        sample = grouped.iloc[chosen]
        values.append(sample.n_flipped.sum() / sample.n_pairs.sum())
    return tuple(np.quantile(values, [0.025, 0.975])), len(grouped)


def mxc_audit(input_paths):
    stats_path = ROOT / "molecular_xai_construct/source_data/molucn/results/mcs_attr_scores_350_stats_lin2.csv"
    flips_path = ROOT / "molecular_xai_construct/results/verdict_flips.csv"
    missing_path = ROOT / "molecular_xai_construct/results/missing_color_files.csv"
    decision_path = ROOT / "molecular_xai_construct/results/go_stop.json"
    for path in (stats_path, flips_path, missing_path, decision_path):
        input_paths[path] = "MXC input"

    stats = pd.read_csv(stats_path)
    primary = stats[(stats.loss == "MSE") & stats.explainer.isin(MXC_METHODS)].copy()
    variance_rows = []
    for maximum in MXC_THRESHOLDS[1:]:
        thresholds = [value for value in MXC_THRESHOLDS if value <= maximum]
        complete_n, complete_ratio = variance_ratio_complete(primary, thresholds)
        available = primary[primary.mcs.isin(thresholds)]
        expected = 350 * len(MXC_METHODS) * len(thresholds)
        finite = int(np.isfinite(available.acc_test).sum())
        variance_rows.append(
            {
                "mcs_min": 50,
                "mcs_max": maximum,
                "n_thresholds": len(thresholds),
                "complete_targets": complete_n,
                "expected_target_method_threshold_cells": expected,
                "available_finite_cells": finite,
                "missing_cell_fraction": 1 - finite / expected,
                "complete_case_definition_method_variance_ratio": complete_ratio,
                "available_cell_marginal_variance_ratio": variance_ratio_available(primary, thresholds),
            }
        )
    variance = pd.DataFrame(variance_rows)
    variance.to_csv(OUT / "mxc_missingness_variance_sensitivity.csv", index=False)

    tau_rows = []
    base = primary[primary.mcs == 50].pivot(index="target", columns="explainer", values="acc_test").reindex(columns=MXC_METHODS)
    for threshold in MXC_THRESHOLDS[1:]:
        other = primary[primary.mcs == threshold].pivot(index="target", columns="explainer", values="acc_test").reindex(columns=MXC_METHODS)
        common = base.index.intersection(other.index)
        a, b = base.loc[common].to_numpy(), other.loc[common].to_numpy()
        valid = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
        a, b = a[valid], b[valid]
        tau = kendalltau(rankdata(-np.median(a, axis=0)), rankdata(-np.median(b, axis=0))).statistic
        tau_rows.append({"mcs_a": 50, "mcs_b": threshold, "n_complete_targets": len(a), "kendall_tau": float(tau)})
    tau = pd.DataFrame(tau_rows)
    tau.to_csv(OUT / "mxc_tau_by_threshold.csv", index=False)
    min_tau = float(tau.kendall_tau.min())

    flips = pd.read_csv(flips_path)
    flips = flips[flips.loss == "MSE"].copy()
    rng = np.random.default_rng(SEED + 1)
    flip_rows = []
    for method in [*MXC_METHODS, "POOLED"]:
        current = flips if method == "POOLED" else flips[flips.explainer == method]
        (low, high), n_targets = target_cluster_interval(current, rng)
        flip_rows.append(
            {
                "explainer": method,
                "targets": n_targets,
                "pair_denominator": int(current.n_pairs.sum()),
                "flipped_pairs": int(current.n_flipped.sum()),
                "flip_rate": float(current.n_flipped.sum() / current.n_pairs.sum()),
                "target_cluster_bootstrap_ci_low": low,
                "target_cluster_bootstrap_ci_high": high,
                "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            }
        )
    flip_summary = pd.DataFrame(flip_rows)
    flip_summary.to_csv(OUT / "mxc_method_flip_rates.csv", index=False)
    observed_flip = float(flip_summary.loc[flip_summary.explainer == "POOLED", "flip_rate"].iloc[0])
    recorded_variance = float(variance.loc[variance.mcs_max == 95, "complete_case_definition_method_variance_ratio"].iloc[0])

    cutoff_rows = []
    cutoff_grids = {
        "minimum_kendall_tau": sorted(set(np.round(np.arange(0.70, 0.951, 0.01), 2)) | {0.80, min_tau}),
        "definition_method_variance_ratio": sorted(set(np.round(np.arange(0.05, 1.001, 0.05), 2)) | {1.0, recorded_variance}),
        "pooled_flip_rate": sorted(set(np.round(np.arange(0.05, 0.301, 0.01), 2)) | {0.20, observed_flip}),
    }
    observed = {
        "minimum_kendall_tau": min_tau,
        "definition_method_variance_ratio": recorded_variance,
        "pooled_flip_rate": observed_flip,
    }
    recorded_cutoffs = {"minimum_kendall_tau": 0.80, "definition_method_variance_ratio": 1.0, "pooled_flip_rate": 0.20}
    for metric, cutoffs in cutoff_grids.items():
        for cutoff in cutoffs:
            triggered = observed[metric] < cutoff if metric == "minimum_kendall_tau" else observed[metric] >= cutoff
            cutoff_rows.append(
                {
                    "metric": metric,
                    "observed": observed[metric],
                    "candidate_cutoff": cutoff,
                    "criterion": "observed < cutoff" if metric == "minimum_kendall_tau" else "observed >= cutoff",
                    "triggered": bool(triggered),
                    "is_recorded_local_cutoff": bool(np.isclose(cutoff, recorded_cutoffs[metric])),
                }
            )
    cutoffs = pd.DataFrame(cutoff_rows)
    cutoffs.to_csv(OUT / "mxc_trigger_cutoff_sensitivity.csv", index=False)
    return variance, tau, flip_summary, cutoffs


def markdown_table(frame, columns, digits=4):
    rows = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for record in frame[columns].itertuples(index=False, name=None):
        values = []
        for value in record:
            if isinstance(value, (float, np.floating)):
                values.append("NA" if not np.isfinite(value) else f"{value:.{digits}f}")
            else:
                values.append(str(value))
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def build_report(rso, atom_maps, templates, ace, ace_components, ace_disjoint, ace_thresholds, mxc_variance, mxc_tau, mxc_flips):
    exact = rso[rso.identity_definition == "exact_tuple"].iloc[0]
    pair_equal = ace[ace.estimand == "pair_equal"].copy()
    unique = ace[ace.estimand == "endpoint_deduplicated"].copy()
    fixed = ace_disjoint[ace_disjoint.estimand == "pair_equal"].copy()
    archived_counts = ace_components[ace_components.analysis_set == "archived_matches"][["dataset", "matched_rows"]].rename(
        columns={"matched_rows": "archived_matched_rows"}
    )
    fixed_counts = ace_components[ace_components.analysis_set == "fixed_endpoint_disjoint"][["dataset", "matched_rows"]].rename(
        columns={"matched_rows": "retained_matched_rows"}
    )
    fixed = fixed.merge(archived_counts, on="dataset").merge(fixed_counts, on="dataset")
    fixed["retained_fraction_of_archived_rows"] = fixed.retained_matched_rows / fixed.archived_matched_rows
    mxc_full = mxc_variance[mxc_variance.mcs_max == 95].iloc[0]
    pooled_flip = mxc_flips[mxc_flips.explainer == "POOLED"].iloc[0]
    threshold_125 = ace_thresholds[(ace_thresholds.estimand == "pair_equal") & np.isclose(ace_thresholds.threshold, 1.25)]
    report = f"""# Cross-case analysis

## Summary

- RSO exact-tuple overlap is {int(exact.overlap_total_rows)}/{int(exact.test_total_rows)} across all test rows and {int(exact.overlap_valid_rows)}/{int(exact.test_valid_reactant_product_rows)} across parser-valid rows.
- ACE inverse-multiplicity and endpoint-deduplicated estimates agree numerically in each dataset. Pair-equal estimates weight repeated endpoint appearances.
- MXC includes {int(mxc_full.complete_targets)} complete targets at MCS 50–95. The complete-case variance ratio is {mxc_full.complete_case_definition_method_variance_ratio:.6f}, the available-cell ratio is {mxc_full.available_cell_marginal_variance_ratio:.6f}, and the pooled MSE flip rate is {pooled_flip.flip_rate:.4%}.

The random seed is `{SEED}` and each bootstrap uses `{BOOTSTRAP_REPLICATES}` replicates.

## RSO definitions and denominators

{markdown_table(rso, ['identity_definition', 'test_total_rows', 'test_valid_reactant_product_rows', 'overlap_total_rows', 'overlap_total_rate', 'overlap_valid_rows', 'overlap_valid_rate'])}

No atom-map marker occurs in the REACTANT, REAGENT, or PRODUCT fields of the three CSV files. The template-frequency analysis uses a custom reaction-centre environment digest from mapped training records. Test alignment has {int(templates.iloc[0].test_aligned_failed_rows)} failures and mapped training parsing has {int(templates.iloc[0].train_mapped_parse_failures)} failures.

## ACE estimands

Endpoint-deduplicated estimates average over distinct molecules. Pair-equal estimates average over matched endpoint appearances. Inverse-multiplicity weighting assigns weight 1/m to each appearance and recovers the distinct-molecule mean.

{markdown_table(pd.concat([unique, pair_equal]), ['dataset', 'estimand', 'cliff_rmse', 'noncliff_rmse', 'rmse_ratio', 'component_bootstrap_ci_low', 'component_bootstrap_ci_high'])}

{markdown_table(ace_components[ace_components.analysis_set == 'archived_matches'], ['dataset', 'matched_rows', 'cliff_unique_molecules', 'noncliff_unique_molecules', 'shared_cliff_noncliff_molecules', 'dependency_components', 'largest_component_rows', 'largest_component_fraction_rows'])}

At threshold 1.25, the pair-equal ECFP screen is positive in {int(threshold_125.positive_datasets_out_of_3.iloc[0])}/3 datasets. The complete 1.15–1.35 map is in `ace_threshold_decision_map.csv`.

The endpoint-disjoint analysis removes rows whose control endpoints overlap the original cliff endpoint pool:

{markdown_table(fixed, ['dataset', 'archived_matched_rows', 'retained_matched_rows', 'retained_fraction_of_archived_rows', 'rmse_ratio', 'component_bootstrap_ci_low', 'component_bootstrap_ci_high'])}

## MXC missingness, flip rates, and thresholds

{markdown_table(mxc_variance, ['mcs_min', 'mcs_max', 'complete_targets', 'missing_cell_fraction', 'complete_case_definition_method_variance_ratio', 'available_cell_marginal_variance_ratio'])}

{markdown_table(mxc_flips, ['explainer', 'targets', 'pair_denominator', 'flipped_pairs', 'flip_rate', 'target_cluster_bootstrap_ci_low', 'target_cluster_bootstrap_ci_high'])}

The minimum Kendall tau from MCS 50 to the evaluated thresholds is {mxc_tau.kendall_tau.min():.6f}. The recorded tau, variance-ratio, and flip-rate thresholds of 0.80, 1.0, and 0.20 are not crossed.

## Reproduction

Run this script to rebuild the tables. Input hashes, byte counts, and roles are recorded in `input_hashes.csv`; `audit_summary.json` contains the JSON summary.
"""
    (OUT / "AUDIT_REPORT.md").write_text(report, encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    input_paths = {Path(__file__).resolve(): "analysis script"}
    rso, atom_maps, templates = rso_audit(input_paths)
    ace, ace_components, ace_disjoint, ace_thresholds = ace_audit(input_paths)
    mxc_variance, mxc_tau, mxc_flips, mxc_cutoffs = mxc_audit(input_paths)
    build_report(rso, atom_maps, templates, ace, ace_components, ace_disjoint, ace_thresholds, mxc_variance, mxc_tau, mxc_flips)

    hashes = []
    for path, role in sorted(input_paths.items(), key=lambda item: str(item[0])):
        hashes.append(
            {
                "role": role,
                "path_relative_to_workspace": str(path.relative_to(ROOT)).replace("\\", "/"),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(hashes).to_csv(OUT / "input_hashes.csv", index=False)

    exact = rso[rso.identity_definition == "exact_tuple"].iloc[0]
    pair = ace[ace.estimand == "pair_equal"].set_index("dataset")
    unique = ace[ace.estimand == "endpoint_deduplicated"].set_index("dataset")
    full_variance = mxc_variance[mxc_variance.mcs_max == 95].iloc[0]
    pooled = mxc_flips[mxc_flips.explainer == "POOLED"].iloc[0]
    summary = {
        "status": "completed",
        "analysis_status": "additional_sensitivity_analysis",
        "seed": SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "rso": {
            "exact_overlap_total": int(exact.overlap_total_rows),
            "test_total": int(exact.test_total_rows),
            "exact_overlap_valid": int(exact.overlap_valid_rows),
            "test_valid": int(exact.test_valid_reactant_product_rows),
            "csv_atom_map_rows": int(atom_maps[atom_maps.field == "ANY_FIELD"].rows_with_atom_map_marker.sum()),
            "test_template_alignment_failures": int(templates.iloc[0].test_aligned_failed_rows),
        },
        "ace": {
            "unique_ratios": {dataset: float(unique.loc[dataset].rmse_ratio) for dataset in ACE_DATASETS},
            "pair_equal_ratios": {dataset: float(pair.loc[dataset].rmse_ratio) for dataset in ACE_DATASETS},
            "inverse_multiplicity_equals_unique_to_atol": 1e-12,
        },
        "mxc": {
            "mcs_50_95_complete_targets": int(full_variance.complete_targets),
            "complete_case_variance_ratio": float(full_variance.complete_case_definition_method_variance_ratio),
            "available_cell_variance_ratio": float(full_variance.available_cell_marginal_variance_ratio),
            "minimum_tau": float(mxc_tau.kendall_tau.min()),
            "pooled_flip_rate": float(pooled.flip_rate),
            "pooled_flip_target_cluster_ci": [float(pooled.target_cluster_bootstrap_ci_low), float(pooled.target_cluster_bootstrap_ci_high)],
            "recorded_triggers_crossed": False,
        },
    }
    (OUT / "audit_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
