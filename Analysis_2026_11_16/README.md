# CLIF analysis — November 16, 2026

This standalone copy retains scripts 01–08, the five-method analysis and local CCW/meta-analysis. Optional script 09 diagnoses prior-year hospitalization associations; optional script 10 constructs a historical NCO panel and explores NN/MNN calibration. The previous Analysis_2026_11_15 and its baseline-S/F sensitivity remain unchanged.

## Optional historical NCO calibration (10)

Run after updated 07 and 09, using the same RUN_ROOT, expanded history, and
optional history-coverage start date:

```bash
Rscript 10_historical_nco_calibration.R "$RUN_ROOT" \
  /projects/emulate_icu/rclif_2_1/clif_hospitalization.parquet
# Add YYYY-MM-DD as argument 3 only if a system history-coverage start is known.
```

Chad has replaced that table with expanded encounters. Script 10 constructs
any admission, >=2 admissions, a stay >=7 days, ED-origin admission, and
SNF/rehab discharge in the year preceding index block admission. Historical
category fields are optional; missing fields/uncertain categories remain
unknown, never automatic negatives. Inspect the source-category mapping audit.
ICU history is not inferred without historical ADT.

Results live in `aggregate/nco_calibration/`; individual histories and model
objects stay in `sensitive/nco_calibration/`. Main files:

- `nco_calibration_inputs.csv`: final-grace-interval clipped-weight historical
  log prevalence ratios, SEs, counts, and inclusion/exclusion reasons.
- `nn_mnn_fit_diagnostics.csv` and `nn_mnn_bias_distributions.png`: NN/MNN
  parameters and fitted systematic-error densities; nonconverged MNN never
  silently becomes NN.
- `nco_patient_level_correlations.csv`, `nco_mnn_membership.csv`, and
  `leave_one_control_out.csv`: overlap, model membership, and stability checks.
- `exploratory_calibrated_ccw_rr.csv` and `exploratory_ccw_rr_forest.png`:
  primary CCW logRR minus fitted majority-component mean, exponentiated to RR.

`NCO_MIN_CELL_N` defaults to 10 observed positive and 10 negative patients in
each arm. Exact duplicate patient-level controls are excluded. NN requires
three usable controls; small-panel exploratory MNN requires four. These are
computational guards, not evidence that the panel is scientifically sufficient.
The controls are overlapping utilization measures, not independent draws.

**All calibration is an assumption-dependent sensitivity analysis.** Historical
prevalence contrasts in selected adherence/survival risk sets can reflect
selection as well as confounding. Transport to day-28 mortality log-RR bias,
and invariance of bias with true effect size, are not established. Consequently
`calibration_supported` remains FALSE. The majority component is designated
valid by assumption, not by proximity to zero or clinical validation.
Plug-in intervals combine primary sampling variance with between-NCO bias
dispersion, omitting calibration-parameter uncertainty, shared-patient
dependence, fitted-weight uncertainty, and transport uncertainty. No formal
calibrated p-values, RD or RMST correction are produced. No 07/08/09 results
are overwritten. Coverage is explicitly unverified if argument 3 is absent.

Synthetic checks: `Rscript tests/test_nco_calibration.R`; optional integration
argument is a synthetic `tests/test_pipeline.R` result folder, never real data.

## Changes in this version

Script 02 is a targeted intermittent-medication source check. It accepts only
the confirmed `steroid` group label, reports observed group
labels and counts, and summarizes routes and drug/route/dose-unit combinations.
Broad vocabulary discovery across all medication tables is no longer repeated.
The command arguments are unchanged; 03/04 do not depend on 02's profile CSVs.
Treatment eligibility remains defined by 03/04, including drug-specific doses
and routes in 04. Source-profile row counts are not treated-patient counts. The server source
check confirmed `med_group = steroid` with 204,291 administration rows; scripts
02--04 no longer accept the plural group label.

03 filters the confirmed `iv`/`enteral` routes and retains administration
actions for auditing. 04 requires `mar_action_category = given`; `not_given`
and unknown actions cannot establish treatment initiation. See
`steroid_administration_action_counts.csv` in 02 and
`steroid_administration_action_audit.csv` in 04. Drug/dose/route thresholds
are unchanged, and weight-based doses remain excluded. Rerun 02--10 after
this update using a fresh run with the existing 01 cohort outputs copied in,
or rerun the full pipeline in a fresh run. Old estimates have not been updated.
CPAP remains grouped with NIPPV for covariates but does not independently
qualify eligibility in 01; resolving acute versus chronic CPAP requires review.

Unmapped-site hospital blocks (including `hopkins_unspecified`, blank/missing IDs,
and unknown/unmapped/unspecified placeholder labels) are excluded in script 01,
after site resolution and before selecting one block per patient. Attrition
includes this exclusion. They do not enter analysis, site tables, or plots.
Scripts 05, 07, 08, and 09 reject stale inputs that still contain these labels.
After updating, run 01--08 in a fresh RUN_ROOT; rerun optional 09 afterward.
The previously downloaded 4,482-patient estimates have not been recomputed.

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

For other modeled measurement inputs (excluding SpO2/FiO2/SF), first retain the latest valid prior reading within the original lookback (12h vitals/support/FiO2; 24h WBC/glucose/GCS; 48h albumin/platelets/bilirubin/creatinine). For remaining inpatient gaps, use the nearest valid same-block measurement in either direction, with ties favoring earlier readings. No maximum fallback distance is imposed. Source measurements must fall between block admission and discharge, no later than the primary event or administrative cutoff. New adjacent fills are not made after discharge. No measurement from a different patient/block is borrowed.

**Future/post-treatment readings can supply these adjacent fills.** This implements the requested exploratory filling rule and is not a strictly prospective covariate definition. Read `adjacent_fill_audit_<variable>.csv`: counts, future-fill counts, distances, fills >24h and remaining gaps. Detailed gap/fill rows are sensitive. A block with no valid reading remains missing; numeric model gaps receive a median plus a missingness indicator. Respiratory support uses its explicit Unknown category. SOFA remains a five-component non-respiratory proxy and is missing when any component is unavailable after component filling.

S/F now uses strict chronological measurement pairing, not adjacent future fills. At each valid SpO2 timestamp, FiO2 is carried forward for at most four hours within the same continuous documented device episode. Device labels are forward-filled only. A nonblank device change starts a new episode, including a return to a previous device. Blank device rows inherit the last documented device; unknown initial device is not backfilled. At a simultaneous timestamp, the most intensive mapped device wins and FiO2 is taken only from that device (or a blank-device row if no device is specified at that timestamp).

At each six-hour decision time use the latest prior SpO2 observation no more than four hours old, and its already paired FiO2. A device change after that SpO2 invalidates the pair. SpO2 >97% remains undefined, rather than selecting an older lower reading. The FiO2 four-hour age is measured at SpO2 time; the SpO2 four-hour recency is a separate rule at decision time. Thus FiO2 can be older than four hours at the final decision if SpO2 itself has been carried forward. All timestamps are retained for audit. The existing script-01 HFNO cohort eligibility remains separate: it uses the original symmetric +/-4-hour clinical qualification window. This covariate update does not redefine time zero or require simultaneous SpO2 at HFNO initiation. Thus the no-future rule here governs IPCW covariates, not the unchanged cohort-eligibility window. A fully prospective redesign of HFNO eligibility would require a separately agreed time-zero definition.

The IPCW model uses S/F as a nominal categorical factor represented by separate indicators: **Undefined; >300 (reference); (250,300]; (200,250]; (120,200]; <=120**. Identifiers are not entered as an ordinal numeric score. No continuous S/F model, median S/F imputation, or separate sf_ratio_missing term is used; Undefined is its own indicator. Current and lagged category indicators enter initiation models. These are SOFA-2-inspired S/F bands, not the complete SOFA-2 respiratory score, which has additional support/ECMO rules.

Coverage outputs distinguish `sf_measurement_pairing_coverage.csv` (time/device pairing and eligible measured S/F), `sf_baseline_coverage.csv`, `sf_decision_time_coverage.csv`, and `sf_any_first_day_coverage.csv` (patients with any eligible first-day value, not a baseline confounder). `sf_measurement_pairs.parquet` in sensitive stores individual timestamps/episodes. Other variables retain the adjacent-fill rule and applicable missingness indicators; arterial pH and PaCO2 remain excluded.

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

## 09 — Primary historical NCO diagnostic

After 07 completes, run:

```bash
Rscript 09_prior_hospitalization_nco_diagnostic.R \
  "/path/to/completed/run" \
  "/path/to/expanded_clif_hospitalization.parquet"
```

An optional third argument gives a verified **system-level** history coverage start (`YYYY-MM-DD`). With that argument, diagnostic contrasts exclude patients whose entire 365-day lookback is not covered. Without it, zeros mean no observed prior within-system admission; complete history coverage is explicitly unverified. Do not use each patient's first encounter as the coverage-start date or require a previous admission to establish eligibility. Do not pass the original filtered index-only table as though it were expanded history.

The single NCO is **any distinct within-system hospital episode starting in the 365 days before index hospital-block admission**. The end of the window is index admission, not respiratory-failure time zero, so the current admission's pre-time-zero days cannot count as prior utilization. Encounters are deduplicated and overlap/<=6-hour transfer chains merged. Expanded episodes containing any selected index encounter are excluded, including connected prior transfer components. Only completed prior episodes count; potentially historical records with invalid dates make an otherwise negative NCO unknown. These records are reported rather than silently coded zero. Future admissions do not count. This definition counts hospitalizations represented by the delivered table; verify that it represents the intended admissions rather than an unfiltered mixture including outpatient/ED-only visits.

**This is not a replacement 28-day survival outcome.** The hospitalization already occurred before follow-up. 09 reuses the exact saved 07 clone risk sets, raw/clipped IPCW and confounder setup, then compares the historical binary NCO's weighted prevalence across the strategy-specific risk sets during grace. It also reports unweighted full-cohort observed-group comparisons and a saved-weight landmark comparison. The NCO is not included in the weight model. No weights or primary outcome estimates are refitted or calibrated.

CCW comparisons condition on evolving adherence/survival risk sets; the landmark comparison conditions on survival through grace. An imbalance can therefore reflect **selection as well as residual confounding**. These contrasts are diagnostics, not a formal causal NCO effect with a guaranteed zero null under every selection mechanism. Initial identical clones have equal baseline NCOs by construction, so initial equality is not validation. CIs combine both clones' contributions at the original-patient level, treating estimated weights as fixed. No bootstrap, multiple-NCO model or bias calibration is included.

Outputs:

- `aggregate/nco_prior_hospitalization/primary_nco_summary.csv`: pooled descriptive/last-grace-interval/landmark contrasts.
- `primary_nco_weighted_contrasts.csv`: site and pooled contrasts over grace, both raw/clipped weights, arm prevalences, prevalence differences and ESS.
- `primary_nco_ccw_diagnostic.png`: diagnostic differences over grace.
- `nco_availability_by_site.csv`, `historical_source_audit.csv`, `METHOD_NOTES.txt`: history availability/coverage and interpretation.
- `sensitive/nco_prior_hospitalization/`: the patient NCO and qualifying historical episode audit.

The existing 07/08 results remain unchanged. You can run 09 alone on a completed run. For optional execution after the complete command file, set `EXPANDED_HOSPITALIZATION_PATH` and, if known, `HISTORY_COVERAGE_START_DATE` before launching that file.

Software check: `Rscript tests/test_primary_nco.R` for boundary/coverage/clone-covariance tests. Add a synthetic completed run path for integration testing.
