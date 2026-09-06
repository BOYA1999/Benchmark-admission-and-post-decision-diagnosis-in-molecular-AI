import argparse
import hashlib
import json
import operator
from pathlib import Path
OPS = {"==": operator.eq, "!=": operator.ne, "<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge}
COMPONENTS = ("population", "estimand", "implementation", "contract", "test")
TO_STATE = {"GO": "ADMITTED", "STOP": "STOPPED", "REVISE": "REVISION_REQUIRED"}
def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
def digest(value):
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest().upper()
def require_keys(value, keys, label):
    missing = [key for key in keys if key not in value]
    if missing:
        raise ValueError(f"{label} missing: {', '.join(missing)}")

def decide(contract, observation):
    require_keys(contract, ["contract_id", "analysis_role", "audit_layer", "claim", "gate", "initial_state", "decision_summary"], "contract")
    require_keys(observation, ["observation_id", "contract_id", "recorded_at", "status", "metric", "value", "prior_state", "component_status"], "observation")
    if observation["contract_id"] != contract["contract_id"]:
        raise ValueError("contract_id mismatch")
    if observation["prior_state"] != contract["initial_state"]:
        raise ValueError("prior_state does not match contract initial_state")
    gate = contract["gate"]
    if observation["metric"] != gate["metric"]:
        raise ValueError("observation metric does not match gate metric")
    statuses = observation["component_status"]
    if set(statuses) != set(COMPONENTS) or any(statuses[key] not in {"retained", "revised", "replaced"} for key in COMPONENTS):
        raise ValueError("component_status must classify all five components exactly once")
    if observation["status"] not in {"complete", "incomplete"}:
        raise ValueError("status must be complete or incomplete")
    if observation["status"] == "incomplete":
        passed, decision = None, "REVISE"
        reason = "Observation is incomplete; the gate was not evaluated."
    else:
        passed = OPS[gate["operator"]](observation["value"], gate["threshold"])
        decision = "GO" if passed else gate["fail_decision"]
        relation = "satisfies" if passed else "does not satisfy"
        reason = f"{gate['metric']}={observation['value']} {relation} {gate['operator']} {gate['threshold']}."
    if contract["analysis_role"] == "post_stop_diagnostic":
        decision, reason = "REVISE", reason + " Post-decision diagnostic evidence cannot alter the prior decision."
    grouped = {status: [key for key in COMPONENTS if statuses[key] == status] for status in ("retained", "revised", "replaced")}
    return {
        "schema_version": "1.0",
        "transition_id": f"{contract['contract_id']}:{observation['observation_id']}",
        "contract_id": contract["contract_id"],
        "observation_id": observation["observation_id"],
        "recorded_at": observation["recorded_at"],
        "analysis_role": contract["analysis_role"],
        "audit_layer": contract["audit_layer"],
        "claim_id": contract["claim"]["claim_id"],
        "from_state": observation["prior_state"],
        "to_state": observation["prior_state"] if contract["analysis_role"] == "post_stop_diagnostic" else TO_STATE[decision],
        "decision": decision,
        "gate_evaluation": {"metric": gate["metric"], "observed": observation["value"], "operator": gate["operator"], "threshold": gate["threshold"], "passed": passed},
        "reason": reason,
        "decision_summary": contract["decision_summary"][decision],
        "components": grouped,
        "input_sha256": {"contract": digest(contract), "observation": digest(observation)},
    }

def check_expected(record, expected):
    for key in ("decision", "to_state", "decision_summary", "components"):
        assert record[key] == expected[key], f"unexpected {key}"

def self_test(root):
    contract = load(root / "fixtures/minimal_contract.json")
    observation = load(root / "fixtures/minimal_observation.json")
    check_expected(decide(contract, observation), load(root / "fixtures/minimal_expected.json"))
    stopped = dict(observation, value=1, observation_id="fixture-stop")
    assert decide(contract, stopped)["decision"] == "STOP"
    incomplete = dict(observation, status="incomplete", value=None, observation_id="fixture-revise")
    assert decide(contract, incomplete)["decision"] == "REVISE"
    assert decide(dict(contract, analysis_role="post_stop_diagnostic"), observation)["to_state"] == "STOPPED"
    rso_contract = load(root / "examples/rso_audit_contract.json")
    rso_observation = load(root / "examples/rso_observation.json")
    check_expected(decide(rso_contract, rso_observation), load(root / "examples/rso_expected.json"))
    print("self-test: PASS (GO, STOP, incomplete REVISE, diagnostic preservation, RSO STOP)")

def main():
    parser = argparse.ArgumentParser(description="Evaluate one audit contract against one observation.")
    parser.add_argument("--contract")
    parser.add_argument("--observation")
    parser.add_argument("--output")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if args.self_test:
        self_test(root)
        return
    if not args.contract or not args.observation:
        parser.error("--contract and --observation are required unless --self-test is used")
    result = decide(load(args.contract), load(args.observation))
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    Path(args.output).write_text(text, encoding="utf-8") if args.output else print(text, end="")

if __name__ == "__main__":
    main()
