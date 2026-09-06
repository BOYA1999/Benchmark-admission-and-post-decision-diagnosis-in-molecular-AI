# Benchmark admission and MolXAI-CRC

This repository provides code, frozen contracts, aggregate results, provenance records, and verification utilities for a molecular-AI benchmark study. The study combines retrospective audits of three completed cases with a separately specified conformal risk-control application.

The retrospective component covers Reaction-Source Overlap, Activity-Cliff Estimand, and Molecular-XAI Construct. The MolXAI-CRC component evaluates omission risk for non-empty task-derived atom reference masks together with retained atom cost. The Molecular-XAI Construct audit and MolXAI-CRC share a methodological setting but have different data, models, estimands, and empirical lineages.

## Repository structure

- `audit_cases/`: the three retrospective case archives, additional analyses, operational schemas, and verification code.
- `molxai_crc/`: MolXAI-CRC source code, public contracts, aggregate evidence, source provenance, and release checks.
- `article_results/`: the compact case-and-application summary used to cross-check reported values.
- `verify_public_package.py`: manifest, privacy, path-safety, and nested reproducibility checks.
- `SECURITY_AND_PRIVACY.md`: public-data and privacy handling.
- `THIRD_PARTY_NOTICES.md`: upstream attribution and licence terms.

## Verification

Create an environment from `requirements.txt`, then run:

```bash
python -B verify_public_package.py
```

The command verifies the complete file manifest, checks repository paths and text surfaces, runs the retrospective-audit verifier, and runs the MolXAI-CRC release verifier and core unit tests.
