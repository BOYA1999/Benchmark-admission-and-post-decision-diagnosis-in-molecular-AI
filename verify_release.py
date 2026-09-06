import csv
import hashlib
import math
import os
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "sha256_manifest.csv"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest():
    with MANIFEST.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    names = [row["relative_path"] for row in rows]
    if len(names) != len(set(names)):
        raise SystemExit("FAIL: duplicate manifest paths")
    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path != MANIFEST
    }
    if set(names) != actual:
        raise SystemExit("FAIL: manifest file set differs from archive contents")
    for row in rows:
        path = ROOT / Path(row["relative_path"])
        if path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"].lower():
            raise SystemExit(f"FAIL: manifest mismatch for {row['relative_path']}")
    return len(rows)


def verify_ace_aggregate_consistency():
    directory = ROOT / "activity_cliff_estimand" / "results"
    tables = {}
    for name in ("activity_cliff_estimand_audit", "metrics", "matching_pairs"):
        with (directory / f"{name}.csv").open(encoding="utf-8-sig", newline="") as handle:
            tables[name] = list(csv.DictReader(handle))
    audit = {(row["dataset"], row["environment"]): row for row in tables["activity_cliff_estimand_audit"]}
    metrics = {(row["dataset"], row["environment"]): row for row in tables["metrics"]}
    pairs = {}
    for row in tables["matching_pairs"]:
        pairs.setdefault((row["dataset"], row["environment"]), []).append(row)
    if len(tables["activity_cliff_estimand_audit"]) != 6 or len(tables["metrics"]) != 6 or len(tables["matching_pairs"]) != 495 or set(audit) != set(metrics) or set(audit) != set(pairs):
        raise SystemExit("FAIL: ACE aggregate table groups or counts differ")
    for key, row in audit.items():
        metric = metrics[key]
        group = pairs[key]
        checks = [
            len(group) == int(row["matched_pair_rows_per_group"]) == int(metric["n_matched_pairs"]),
            all(pair["representation"] == metric["representation"] == "ecfp" for pair in group),
            len({pair["match_id"] for pair in group}) == len(group),
        ]
        for label in ("cliff", "noncliff"):
            endpoints = {pair[f"{label}_{side}"] for pair in group for side in ("i", "j")}
            checks.extend([
                int(row[f"{label}_endpoint_appearances"]) == 2 * len(group),
                int(row[f"unique_{label}_molecules"]) == len(endpoints),
                math.isclose(float(row[f"stored_unique_{label}_rmse"]), float(row[f"recomputed_unique_{label}_rmse"]), rel_tol=0, abs_tol=1e-12),
                math.isclose(float(row[f"stored_unique_{label}_rmse"]), float(metric[f"matched_{label}_rmse"]), rel_tol=0, abs_tol=1e-12),
            ])
        for prefix in ("stored_unique", "recomputed_unique", "pair_equal"):
            ratio = float(row[f"{prefix}_cliff_rmse"]) / float(row[f"{prefix}_noncliff_rmse"])
            checks.append(math.isclose(float(row[f"{prefix}_ratio"]), ratio, rel_tol=0, abs_tol=1e-12))
        checks.append(math.isclose(float(row["stored_unique_ratio"]), float(metric["cliff_ratio"]), rel_tol=0, abs_tol=1e-12))
        for flag, ratio in (("stored_gate_pass", "stored_unique_ratio"), ("pair_equal_gate_pass", "pair_equal_ratio")):
            checks.append(row[flag].lower() == str(float(row[ratio]) >= float(row["gate_threshold"])).lower())
        if not all(checks):
            raise SystemExit(f"FAIL: ACE aggregate consistency for {key}")
    print("PASS: ACE aggregate consistency; 6 groups and 495 match rows; endpoint-level prediction recomputation requires external inputs")


def verify_mxc_reconciliation():
    directory = ROOT / "additional_analyses" / "mxc_sentinel_reconciliation"
    with (directory / "reconciliation_summary.json").open(encoding="utf-8-sig") as handle:
        import json
        summary = json.load(handle)
    raw = summary.get("raw_gate_inputs", {})
    if (
        summary.get("sentinel_cells") != 40
        or summary.get("sentinel_exact") != 6
        or summary.get("sentinel_within_0_02") != 15
        or summary.get("raw_recomputed_decision") != "STOP"
        or raw.get("flip_n_flipped") != 96402
        or raw.get("flip_n_pairs") != 647446
        or raw.get("flip_target_method_counts_exactly_match_archived") is not True
    ):
        raise SystemExit("FAIL: MXC sentinel reconciliation summary is incomplete or inconsistent")


count = verify_manifest()
verify_ace_aggregate_consistency()
verify_mxc_reconciliation()
environment = os.environ.copy()
environment["PYTHONDONTWRITEBYTECODE"] = "1"
commands = (
    [sys.executable, str(ROOT / "molecular_xai_construct" / "src" / "verify_archived_results.py")],
    [sys.executable, str(ROOT / "operational_framework" / "audit_decision.py"), "--self-test"],
)
for command in commands:
    subprocess.run(command, cwd=ROOT, env=environment, check=True)
verify_manifest()

privacy_patterns = (
    re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]"),
    re.compile(r"/(?:home|Users)/[^/\s]+/", re.I),
    re.compile(r"(?<![A-Za-z0-9_.+-])[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}(?![A-Za-z0-9_.-])"),
    re.compile("BEGIN " + "PRIVATE KEY"),
    re.compile("gh" + r"[pousr]_[A-Za-z0-9]{20,}"),
    re.compile("AKIA" + r"[A-Z0-9]{16}"),
    re.compile("sk" + r"-[A-Za-z0-9_-]{20,}"),
)
for path in ROOT.rglob("*"):
    if path.is_file() and path.suffix.lower() in {".md", ".txt", ".csv", ".json", ".py", ".yml"}:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
        if any(pattern.search(text) for pattern in privacy_patterns):
            raise SystemExit(f"FAIL: privacy or credential pattern in {path.relative_to(ROOT).as_posix()}")

print(f"PASS: {count} manifested files; ACE aggregate consistency, MXC, sentinel reconciliation, operational semantics, integrity, and privacy scan verified")
