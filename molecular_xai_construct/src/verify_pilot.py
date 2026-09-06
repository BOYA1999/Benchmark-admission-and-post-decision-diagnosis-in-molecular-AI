import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def digest(path, algorithm):
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main():
    baseline = pd.read_csv(RESULTS / "baseline_reproduction.csv")
    stability = pd.read_csv(RESULTS / "rank_stability.csv")
    variance = pd.read_csv(RESULTS / "variance_decomposition.csv")
    flips = pd.read_csv(RESULTS / "verdict_flips.csv")
    negative = pd.read_csv(RESULTS / "negative_control.csv")
    normalization = pd.read_csv(RESULTS / "normalization_sensitivity.csv")
    qc = json.loads((RESULTS / "qc_summary.json").read_text(encoding="utf-8"))
    decision = json.loads((RESULTS / "go_stop.json").read_text(encoding="utf-8"))
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
        "all_key_metrics_finite": bool(primary_stability[["kendall_tau", "bootstrap_ci_low", "bootstrap_ci_high"]].notna().all().all()),
    }
    checks["construct_gate"] = checks["primary_rank_signal"] or checks["primary_variance_signal"] or checks["primary_flip_signal"]
    checks["recomputed_decision"] = "STOP" if checks["baseline_gate"] and checks["data_gate"] and checks["negative_control_gate"] and not checks["construct_gate"] else "REVIEW"
    checks["matches_run_decision"] = checks["recomputed_decision"] == decision["preliminary_decision"]
    checks["archive_md5"] = {
        "data": digest(ROOT / "downloads" / "molucn_data.tar.gz", "md5"),
        "colors": digest(ROOT / "downloads" / "molucn_colors.tar.gz", "md5"),
        "results": digest(ROOT / "downloads" / "molucn_results.tar.gz", "md5"),
    }
    checks["analysis_sha256"] = digest(ROOT / "src" / "run_pilot.py", "sha256")
    (RESULTS / "verification.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
    print(json.dumps(checks, indent=2))
    if not checks["matches_run_decision"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
