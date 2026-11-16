# CLIF analysis — November 16, 2026

This fresh standalone copy retains scripts 01–08, the five-method analysis and local CCW/meta-analysis. The previous Analysis_2026_11_15 and its baseline-S/F sensitivity remain unchanged. No NCO or calibration analysis is included.

## Changes in this version

Chad's inclusive screening bounds replace the earlier bounds:

| Variable | Range |
|---|---|
| Temperature | 30–45 C |
| SpO2 | 70–100%; S/F additionally requires <=97% |
| MAP | 0–250 mmHg |
| WBC | 0–500 x10^3/uL |
| Serum glucose | 0–2000 mg/dL |
| Albumin | 0–15 g/dL |
| Platelets | 0–2000 x10^3/uL |
| Total bilirubin | 0–80 mg/dL |
| Creatinine | 0–30 mg/dL |
| GCS | 3–15 |
| FiO2 | 0.21–1; documented percent form normalized |
| S/F | 70–462, with SpO2 <=97% |

`variable_rules.R` is the source of these rules. Outside-range values become missing before temporal selection; source records are never overwritten. Zeros explicitly permitted by these ranges are retained. **No additional percentile winsorization** is performed in 07 or local 08. WBC/glucose retain log1p transformations. HFNO's plausibility lower bound changes from SpO2 80% to 70%, so the eligible cohort and time zero can change; rerun 01–04.

For all modeled measurement inputs, first retain the latest valid prior reading within the original lookback (12h vitals/support/FiO2; 24h WBC/glucose/GCS; 48h albumin/platelets/bilirubin/creatinine). For remaining inpatient gaps, use the nearest valid same-block measurement in either direction, with ties favoring earlier readings. No maximum fallback distance is imposed. Source measurements must fall between block admission and discharge, no later than the primary event or administrative cutoff. New adjacent fills are not made after discharge. No measurement from a different patient/block is borrowed.

**Future/post-treatment readings can supply these adjacent fills.** This implements the requested exploratory filling rule and is not a strictly prospective covariate definition. Read `adjacent_fill_audit_<variable>.csv`: counts, future-fill counts, distances, fills >24h and remaining gaps. Detailed gap/fill rows are sensitive. A block with no valid reading remains missing; numeric model gaps receive a median plus a missingness indicator. Respiratory support uses its explicit Unknown category. SOFA remains a five-component non-respiratory proxy and is missing when any component is unavailable after component filling.

S/F is a time-varying covariate in this full-cohort version. If the selected saturation exceeds 97%, S/F remains unavailable: a newer high saturation is not replaced by an older lower reading to manufacture an eligible ratio. SpO2 >97% does not exclude the patient. Remaining S/F gaps receive median imputation plus current/lagged missingness indicators. Arterial pH and PaCO2 remain excluded.

## Unchanged design

- Default grace: 24 hours; six-hour intervals; initiation within grace versus later/never initiation.
- Outcome: in-hospital death or hospice discharge by day 28 in the selected block. Other alive discharge is not censoring, and subsequent post-discharge deaths are excluded.
- Original age, ED admission, POA pneumonia, support-initiation and shock exclusions; block construction and one selected block per patient are retained.
- Categorical support: room air, conventional oxygen, HFNO, NIPPV/CPAP, IMV and unknown.
- 07 fits the same five methods with balance/weight/ESS diagnostics. Initiation weights freeze after grace; later initiation remains allowed in controls.
- 08 fits local CCW and pools contrasts by REML/Hartung–Knapp, also reporting fixed-effect results and site/six-method forests. Standard HK can give unusually narrow intervals when site estimates are very similar; inference is unchanged from the previous version.
- Influence-function uncertainty treats fitted nuisance models as fixed. There is no Cox model or NCO step.

## Run

Required packages: arrow, data.table, dplyr, ggplot2, survival, metafor.

```bash
bash /projects/emulate_icu/CLIF-Cohort/Analysis_2026_11_16/RUN_PIPELINE_COMMANDS.txt
```

Check the server `cd`, raw data path and administrative cutoff in the command file first. It defaults to `/projects/emulate_icu/rclif_2_1` and 2026-10-31. It runs 01–05, 07 and 08; 06 is sourced automatically. Each run gets its own directory:

```text
results/run_<timestamp>_grace24h/
  aggregate/eda/
  aggregate/analysis/
  aggregate/local_ccw_meta/
  sensitive/cohort/
  sensitive/source_profile/
  sensitive/steroid_candidates/
  sensitive/steroid_treatment/
  sensitive/eda_preparation/
  sensitive/analysis/
  sensitive/local_ccw_meta/
  logs/
```

Do not substitute the expanded history hospitalization table for the original index-encounter input without reviewing its effect on cohort/block construction. This version does not use the expanded table for history/NCOs; keep the intended index-cohort source files in the configured data directory.

Key files include screening_ranges.csv, raw_measurement_screening_counts.csv, adjacent_fill_audit_*.csv, model_relevant_missingness.csv, ipcw_model_specification.csv, five_method_results.csv, and six_method_results.csv. EDA/figures are unsuppressed; patient-level files remain in sensitive directories.

Synthetic checks: `Rscript tests/test_pipeline.R` and `Rscript tests/test_local_meta.R /path/to/synthetic/run`. They do not run on the actual EHR data.
