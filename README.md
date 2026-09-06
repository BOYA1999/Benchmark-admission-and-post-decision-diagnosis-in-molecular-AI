# Retrospective audit cases

This directory contains the three molecular-AI case studies and their executable verification materials:

- `reaction_source_overlap`: source-identity analysis for the reaction benchmark.
- `activity_cliff_estimand`: activity-cliff estimand analysis.
- `molecular_xai_construct`: molecular-XAI construct analysis.
- `operational_framework`: schemas, state-transition logic and a hand-checkable fixture.
- `additional_analyses`: additional cross-case analyses and MXC sentinel reconciliation.

From this directory, run:

`python verify_release.py`

The verifier checks the file manifest, checks ACE aggregate-table consistency, recomputes the packaged MXC gate statistics, runs the operational state-transition self-test, validates the MXC reconciliation summary and scans text files for private local paths.

The ACE estimand table can be rebuilt from externally held molecule-level inputs. From an external working directory, provide the prepared split manifest and the ECFP run directory containing predictions, matches and metrics:

`python <repository>/audit_cases/activity_cliff_estimand/scripts/audit_estimand.py --split-manifest data/split_manifest.csv --run-dir results/ecfp --output results/activity_cliff_estimand_audit.csv`

The distributed `activity_cliff_estimand/results/` contains aggregate summaries and integer-index matching records, not molecular structures, activities or individual predictions. The packaged consistency check does not replace this endpoint-level reconstruction.

The archived MXC gate calculations can be recomputed with:

`python molecular_xai_construct/src/verify_archived_results.py`

The additional analyses are stored under `additional_analyses/`. Their reports, tables and historical input hashes are included. In these provenance records, paths are public logical labels; byte counts and hashes refer to the archived inputs, not necessarily to the edited public code. Current package hashes are in `sha256_manifest.csv`.

Full cross-case and sentinel reruns require external upstream inputs. Use a separate working copy with the case directories named as above. RSO expects `data/raw/USPTO_MIT/MIT_separated/`, `data/processed/` and `results/data_audit.json`; ACE expects `data/split_manifest.csv` and the ECFP predictions, matches and metrics under `results/`; MXC expects its `source_data/`, `reference/molucn/` and archived `results/`. The contracts and input-hash tables identify the sources. Keep this working copy separate from the upload directory. Supplementary figure source tables and plotting scripts are under `supplementary_figures/`.

Key packaged results include the 40,000-row RSO identity audit, endpoint- and pair-weighted ACE summaries, MXC gate recomputation, 40-cell sentinel reconciliation and raw-row gate-input checks.

`environment.txt`, `environment.yml` and `requirements.txt` record the analysis environment and reconstruction dependencies. Repository-authored code and documentation use the root MIT License. Upstream resources and data-derived outputs retain the terms listed in `THIRD_PARTY_NOTICES.md`.
