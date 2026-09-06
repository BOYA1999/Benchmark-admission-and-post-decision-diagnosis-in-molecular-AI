# Major 9 operational artefact audit

## Outcome

This directory implements the evidence-state logic as an inspectable decision layer. It maps a declared claim and estimand to a numeric gate, evaluates the observation, and emits a `GO`, `STOP`, or `REVISE` transition with retained, revised, and replaced component lists.

The RSO demonstration applies the decision layer to the reported source-overlap case.

## Files and roles

| File | Role |
|---|---|
| `audit_contract.schema.json` | Declares the source-to-re-entry audit layer, claim, estimand, gate, decision summaries and five component baselines. |
| `benchmark_manifest.schema.json` | Declares benchmark provenance, assets, split counts and the identity rule. |
| `state_transition.schema.json` | Constrains the emitted decision record and component-change lists. |
| `audit_decision.py` | Standard-library-only evaluator; 100 physical lines. |
| `programme_record_sanitized.json` | Four-direction sampling record with archive-completeness eligibility. |
| `fixtures/` | Four-label fixture whose set intersection can be checked by inspection. |
| `examples/` | RSO contract, benchmark manifest, observation and expected decision. |
| `outputs/` | Generated state-transition records for the fixture and RSO example. |

## Decision semantics

The observation payload must contain `status`, the contract's named `metric`, its numeric `value`, the prior state, and a `component_status` entry for each of `population`, `estimand`, `implementation`, `contract`, and `test`.

- A complete observation satisfying the gate produces `GO` and `ADMITTED`.
- A complete observation failing the gate produces the contract's declared `STOP` or `REVISE` action.
- An incomplete observation produces `REVISE` before gate evaluation.
- A diagnostic observation produces `REVISE` while preserving the prior `STOPPED` state; a separate decision record handles re-entry.
- Every component is classified exactly once as `retained`, `revised`, or `replaced`.
- The output copies the matching decision statement and adds canonical JSON SHA-256 digests of both inputs.

`STOPPED` records the contract-specific branch action and its decision statement. The executable example evaluates exact source separation.

A new decision record states which components changed.

## Claim–estimand–gate–decision checks

### Hand-checkable fixture

The reference labels are `[A, B]` and the evaluation labels are `[C, D]`. Their exact, case-sensitive intersection is empty, so `overlap_count = 0`. The contract is `overlap_count == 0`; therefore the generated re-entry decision is `GO`, transitioning `STOPPED -> ADMITTED`. Population, estimand and contract are retained, implementation is revised, and the decision test is replaced. Changing the observed value to `1` yields `STOP`, while marking the observation incomplete yields `REVISE`.

### Existing RSO demonstration

The RSO contract uses all 40,000 distributed test rows as the decision population. The 39,995 parser-valid reactant/product rows are a separate validity statistic, not a replacement decision denominator. The exact normalized tuple audit found 739 test-row matches to train and 50 to validation, with 7 in both, hence a union of `739 + 50 - 7 = 782`. Because `782` does not satisfy the frozen `test_overlap_reference_union == 0` gate, the output is `STOP`.

The example evaluates exact normalized, reagent-preserving tuple overlap. The distributed strings contain no atom-map tokens.

## Reproduction

From this directory:

```text
python audit_decision.py --self-test
python audit_decision.py --contract fixtures/minimal_contract.json --observation fixtures/minimal_observation.json --output outputs/minimal_transition.json
python audit_decision.py --contract examples/rso_audit_contract.json --observation examples/rso_observation.json --output outputs/rso_transition.json
```

## Verification performed on 2026-09-04

- Draft 2020-12 meta-schema check: PASS for all three schemas.
- Instance validation: PASS for two contracts, two benchmark manifests and two emitted transitions.
- Manifest-reference SHA-256 checks: PASS for both examples.
- Program self-test: `PASS (GO, STOP, incomplete REVISE, diagnostic preservation, RSO STOP)`.
- RSO output: `STOP`, transition `AUDIT_PENDING -> STOPPED`, with all five components recorded as retained.

Runtime evaluation uses only the Python standard library. The independent QA step used `jsonschema` 4.26.0 to validate the schemas and instances.
