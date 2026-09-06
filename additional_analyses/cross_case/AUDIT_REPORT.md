# Cross-case analysis

## Summary

- RSO exact-tuple overlap is 782/40,000 across all test rows and 782/39,995 across parser-valid rows.
- ACE endpoint-deduplicated ratios are 1.083082, 1.145647 and 1.140319; pair-equal ratios are 1.253834, 1.291421 and 1.291296.
- MXC has 233 complete targets across MCS 50–95. The complete-case variance ratio is 0.127583, the available-cell ratio is 0.134492, the minimum Kendall tau is 0.866667 and the pooled MSE flip rate is 14.8896%.

The random seed is `20260904`; bootstrap summaries use 5,000 replicates.

## RSO

`rso_identity_sensitivity.csv` reports exact tuple, reagent-removed and product-level identity definitions. `rso_atom_map_audit.csv` records atom-map counts, and `rso_template_frequency_audit.csv` records the custom reaction-centre environment digest analysis.

## ACE

`ace_estimand_sensitivity.csv` reports endpoint-deduplicated, inverse-multiplicity and pair-equal estimates. `ace_endpoint_multiplicity.csv` and `ace_endpoint_overlap_and_components.csv` record endpoint reuse and dependency components. `ace_endpoint_disjoint_sensitivity.csv` reports the endpoint-disjoint analysis. The 1.15–1.35 threshold map is in `ace_threshold_decision_map.csv`.

At threshold 1.25, the pair-equal ECFP screen is positive in all three datasets.

## MXC

`mxc_missingness_variance_sensitivity.csv` reports complete- and available-cell variance ratios. `mxc_tau_by_threshold.csv` reports rank stability, `mxc_method_flip_rates.csv` reports normalization flip rates, and `mxc_trigger_cutoff_sensitivity.csv` reports the threshold grid.

The recorded tau, variance-ratio and flip-rate thresholds of 0.80, 1.0 and 0.20 are not crossed.

## Reproduction

`run_cross_case_audit.py` rebuilds the tables from the inputs recorded in `input_hashes.csv`. `audit_summary.json` contains the machine-readable summary.
