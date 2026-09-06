import json
from pathlib import Path

import pandas as pd


root = Path(__file__).resolve().parents[1]
results = root / "results"
baseline = pd.read_csv(results / "baseline_reproduction.csv")
stability = pd.read_csv(results / "rank_stability.csv")
variance = pd.read_csv(results / "variance_decomposition.csv")
flips = pd.read_csv(results / "verdict_flips.csv")
negative = pd.read_csv(results / "negative_control.csv")
normalization = pd.read_csv(results / "normalization_sensitivity.csv")
qc = json.loads((results / "qc_summary.json").read_text(encoding="utf-8"))
decision = json.loads((results / "go_stop.json").read_text(encoding="utf-8"))

primary_stability = stability[(stability.source == "lin2") & (stability.loss == "MSE") & (stability.metric == "acc_test")]
primary_variance = variance[(variance.source == "lin2") & (variance.loss == "MSE") & (variance.metric == "acc_test")].iloc[0]
primary_flips = flips[flips.loss == "MSE"]
flip_rate = primary_flips.n_flipped.sum() / primary_flips.n_pairs.sum()
checks = {
    "baseline_gate": bool(baseline.within_0_02.any()),
    "baseline_passed_cells": int(baseline.within_0_02.sum()),
    "data_gate": bool(qc["n_test_pairs_mcs50"] >= 1000 and qc["n_test_targets"] >= 100),
    "negative_control_gate": bool(((negative.absolute_difference_ge_0_02) & negative.ci_excludes_zero).any()),
    "primary_min_kendall_tau": float(primary_stability.kendall_tau.min()),
    "primary_rank_signal": bool((primary_stability.kendall_tau < 0.8).any()),
    "primary_definition_method_variance_ratio": float(primary_variance.definition_method_variance_ratio),
    "primary_variance_signal": bool(primary_variance.definition_method_variance_ratio >= 1),
    "mse_normalization_flip_rate": float(flip_rate),
    "primary_flip_signal": bool(flip_rate >= 0.2),
    "normalization_rows": int(len(normalization)),
    "flip_rows": int(len(flips)),
}
checks["construct_gate"] = checks["primary_rank_signal"] or checks["primary_variance_signal"] or checks["primary_flip_signal"]
checks["recomputed_decision"] = "STOP" if checks["baseline_gate"] and checks["data_gate"] and checks["negative_control_gate"] and not checks["construct_gate"] else "REVIEW"
checks["matches_run_decision"] = checks["recomputed_decision"] == decision["preliminary_decision"]
checks["status"] = "PASS" if checks["matches_run_decision"] else "FAIL"
(results / "archived_result_verification.json").write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8", newline="\n")
print(json.dumps(checks, indent=2))
if checks["status"] != "PASS":
    raise SystemExit(1)
