import argparse
from pathlib import Path

import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(description="Recompute ACE estimands from external molecule-level inputs.")
parser.add_argument("--split-manifest", required=True, type=Path)
parser.add_argument("--run-dir", required=True, type=Path,
                    help="Directory containing predictions_long.csv, matching_pairs.csv and metrics.csv")
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
RUN = args.run_dir
MANIFEST = args.split_manifest
OUT = args.output


def rmse(values):
    return float(np.sqrt(np.mean(values)))


manifest = pd.read_csv(MANIFEST)
predictions = pd.read_csv(RUN / "predictions_long.csv")
matches = pd.read_csv(RUN / "matching_pairs.csv")
stored = pd.read_csv(RUN / "metrics.csv")
rows = []

for (dataset, environment), group in matches.groupby(["dataset", "environment"], sort=True):
    test = manifest[
        (manifest["dataset"] == dataset)
        & (manifest["environment"] == environment)
        & (manifest["split"] == "test")
    ].reset_index(drop=True)
    pred = predictions[
        (predictions["dataset"] == dataset)
        & (predictions["environment"] == environment)
        & (predictions["split"] == "test")
    ]
    ensemble = pred.groupby("canonical_smiles", sort=False)["prediction"].mean()
    yhat = test["canonical_smiles"].map(ensemble).to_numpy(float)
    squared_error = (test["activity"].to_numpy(float) - yhat) ** 2

    cliff_endpoints = group[["cliff_i", "cliff_j"]].to_numpy(int)
    noncliff_endpoints = group[["noncliff_i", "noncliff_j"]].to_numpy(int)
    cliff_unique = np.unique(cliff_endpoints)
    noncliff_unique = np.unique(noncliff_endpoints)

    unique_cliff_rmse = rmse(squared_error[cliff_unique])
    unique_noncliff_rmse = rmse(squared_error[noncliff_unique])
    pair_cliff_rmse = rmse(squared_error[cliff_endpoints].mean(axis=1))
    pair_noncliff_rmse = rmse(squared_error[noncliff_endpoints].mean(axis=1))
    stored_row = stored[(stored["dataset"] == dataset) & (stored["environment"] == environment)].iloc[0]
    assert np.isclose(unique_cliff_rmse, stored_row["matched_cliff_rmse"], rtol=0, atol=1e-12)
    assert np.isclose(unique_noncliff_rmse, stored_row["matched_noncliff_rmse"], rtol=0, atol=1e-12)
    assert np.isclose(unique_cliff_rmse / unique_noncliff_rmse, stored_row["cliff_ratio"], rtol=0, atol=1e-12)

    rows.append(
        {
            "dataset": dataset,
            "environment": environment,
            "matched_pair_rows_per_group": len(group),
            "cliff_endpoint_appearances": cliff_endpoints.size,
            "noncliff_endpoint_appearances": noncliff_endpoints.size,
            "unique_cliff_molecules": len(cliff_unique),
            "unique_noncliff_molecules": len(noncliff_unique),
            "stored_unique_cliff_rmse": stored_row["matched_cliff_rmse"],
            "recomputed_unique_cliff_rmse": unique_cliff_rmse,
            "stored_unique_noncliff_rmse": stored_row["matched_noncliff_rmse"],
            "recomputed_unique_noncliff_rmse": unique_noncliff_rmse,
            "stored_unique_ratio": stored_row["cliff_ratio"],
            "recomputed_unique_ratio": unique_cliff_rmse / unique_noncliff_rmse,
            "pair_equal_cliff_rmse": pair_cliff_rmse,
            "pair_equal_noncliff_rmse": pair_noncliff_rmse,
            "pair_equal_ratio": pair_cliff_rmse / pair_noncliff_rmse,
            "gate_threshold": 1.25,
            "stored_gate_pass": bool(stored_row["cliff_ratio"] >= 1.25),
            "pair_equal_gate_pass": bool(pair_cliff_rmse / pair_noncliff_rmse >= 1.25),
        }
    )

out = pd.DataFrame(rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False, lineterminator="\n")
print(out.to_string(index=False))
print("stored unique-molecule metrics reproduced within absolute tolerance 1e-12")
print(f"saved: {OUT}")
