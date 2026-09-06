import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_pilot import pair_table, risk_metrics, rmse


def close(a, b, tolerance=1e-10):
    return np.isclose(float(a), float(b), rtol=tolerance, atol=tolerance)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--baseline-dir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    data_dir = Path(args.data_dir)
    baseline_dir = Path(args.baseline_dir)
    manifest = pd.read_csv(data_dir / "split_manifest.csv")
    split_qc = pd.read_csv(data_dir / "split_qc.csv")
    metrics = pd.read_csv(baseline_dir / "metrics.csv")
    predictions = pd.read_csv(baseline_dir / "predictions_long.csv")
    uq = pd.read_csv(baseline_dir / "uq_scores_long.csv")
    matches = pd.read_csv(baseline_dir / "matching_pairs.csv")
    errors = []
    balance_rows = []
    recomputed = []
    if len(metrics) != 6:
        errors.append(f"expected 6 metric rows, found {len(metrics)}")
    if predictions.groupby(["dataset", "environment", "split"])["seed"].nunique().ne(3).any():
        errors.append("not all baseline surfaces have three seeds")
    for _, reported in metrics.iterrows():
        dataset, environment = reported["dataset"], reported["environment"]
        frame = manifest[(manifest["dataset"] == dataset) & (manifest["environment"] == environment)].reset_index(drop=True)
        test = frame[frame["split"] == "test"].reset_index(drop=True)
        group = predictions[
            (predictions["dataset"] == dataset)
            & (predictions["environment"] == environment)
            & (predictions["split"] == "test")
        ]
        pivot = group.pivot(index="canonical_smiles", columns="seed", values="prediction")
        prediction = pivot.mean(axis=1).reindex(test["canonical_smiles"]).to_numpy()
        y = test["activity"].to_numpy(float)
        distance = uq[
            (uq["dataset"] == dataset)
            & (uq["environment"] == environment)
            & (uq["method"] == "ecfp_distance")
        ].set_index("canonical_smiles")["uncertainty"].reindex(test["canonical_smiles"]).to_numpy()
        variance = uq[
            (uq["dataset"] == dataset)
            & (uq["environment"] == environment)
            & (uq["method"] == "ensemble_variance")
        ].set_index("canonical_smiles")["uncertainty"].reindex(test["canonical_smiles"]).to_numpy()
        pairs = pair_table(test, distance)
        cliff_pairs = pairs[pairs["pair_type"] == "cliff"]
        match = matches[(matches["dataset"] == dataset) & (matches["environment"] == environment)]
        cliff_indices = match[["cliff_i", "cliff_j"]].to_numpy().ravel().astype(int)
        noncliff_indices = match[["noncliff_i", "noncliff_j"]].to_numpy().ravel().astype(int)
        cliff_rmse = rmse(y[sorted(set(cliff_indices))], prediction[sorted(set(cliff_indices))])
        noncliff_rmse = rmse(y[sorted(set(noncliff_indices))], prediction[sorted(set(noncliff_indices))])
        ratio = cliff_rmse / noncliff_rmse
        checks = {
            "rmse": rmse(y, prediction),
            "mae": float(np.mean(np.abs(y - prediction))),
            "cliff_rmse": cliff_rmse,
            "matched_noncliff_rmse": noncliff_rmse,
            "cliff_ratio": ratio,
        }
        for method, values in (("ensemble_variance", variance), ("ecfp_distance", distance)):
            _, _, aurc, sr80, rho = risk_metrics(y, prediction, values)
            checks[f"{method}_aurc"] = aurc
            checks[f"{method}_sr80"] = sr80
            checks[f"{method}_spearman_abs_error"] = rho
        for key, value in checks.items():
            if not close(value, reported[key]):
                errors.append(f"{dataset}/{environment} mismatch for {key}: {value} vs {reported[key]}")
        expected_pairs = int(
            split_qc[(split_qc["dataset"] == dataset) & (split_qc["environment"] == environment)][
                "test_cliff_pairs"
            ].iloc[0]
        )
        if len(cliff_pairs) != expected_pairs or len(match) != expected_pairs:
            errors.append(f"{dataset}/{environment} cliff or matching denominator mismatch")
        pair_lookup = {(int(row.i), int(row.j)): row for row in pairs.itertuples()}
        columns = ["max_similarity", "avg_mw", "mw_diff", "train_distance_mean", "mmp_edit_size"]
        for column in columns:
            cliff_values, noncliff_values = [], []
            for row in match.itertuples():
                cliff_value = getattr(pair_lookup[(row.cliff_i, row.cliff_j)], column)
                noncliff_value = getattr(pair_lookup[(row.noncliff_i, row.noncliff_j)], column)
                if np.isfinite(cliff_value) and np.isfinite(noncliff_value):
                    cliff_values.append(cliff_value)
                    noncliff_values.append(noncliff_value)
            pooled = np.sqrt((np.var(cliff_values) + np.var(noncliff_values)) / 2) if cliff_values else np.nan
            smd = abs(np.mean(cliff_values) - np.mean(noncliff_values)) / pooled if pooled > 0 else 0.0
            balance_rows.append(
                {
                    "dataset": dataset,
                    "environment": environment,
                    "covariate": column,
                    "matched_pairs_with_value": len(cliff_values),
                    "absolute_smd": smd,
                }
            )
        recomputed.append({"dataset": dataset, "environment": environment, **checks})
        print(f"VERIFIED {dataset} {environment} ratio={ratio:.6f}", flush=True)
    balance = pd.DataFrame(balance_rows)
    balance.to_csv(baseline_dir / "matching_balance.csv", index=False)
    pd.DataFrame(recomputed).to_csv(baseline_dir / "recomputed_metrics.csv", index=False)
    stress = pd.DataFrame(recomputed)
    stress = stress[stress["environment"] == "mmp_series"]
    report = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "verified_surfaces": int(len(recomputed)),
        "maximum_absolute_smd": float(balance["absolute_smd"].max()),
        "mmp_series_cliff_ratios": dict(zip(stress["dataset"], stress["cliff_ratio"])),
        "ecfp_stressor_positive_datasets": int((stress["cliff_ratio"] >= 1.25).sum()),
    }
    Path(args.report).write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if not errors else 1)


if __name__ == "__main__":
    main()
