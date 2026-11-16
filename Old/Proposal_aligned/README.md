# Standalone analysis with the proposal's covariates

This folder is independent of the parent analysis and CAPE_COD_aligned. All
scripts/helpers needed for construction, treatment classification, EDA, five
estimators and diagnostics are included locally. Existing code is unchanged.
The analysis uses the ORIGINAL adult pneumonia cohort and pooled systemic steroid
exposure definition, not the separate CAPE COD/hydrocortisone policy.

## Run on Discovery

From CLIF-Cohort:

```bash
export TREATMENT_GRACE_HOURS=24
export PROPOSAL_MAX_MISSING=0.70
export CLIF_DATA_CUTOFF_DATE=2026-10-31
export SUMMARY_SMALL_CELL_THRESHOLD=10
export CCW_MAX_PATIENTS=0
Rscript Proposal_aligned/run_pipeline.R ../rclif_2_1
```

With no data argument, the source defaults to ../../rclif_2_1 relative to this
folder, so it also works from another working directory. R packages: arrow,
data.table, dplyr, ggplot2. Output is always inside this folder at
`results/<timestamped-run-id>/`. Set PROPOSAL_RUN_ID to choose a name; existing
runs are refused. Per-stage logs identify failures. COMPLETE.txt is created only
after every stage finishes. Use a new run ID for another grace period or cutoff.

## Adjustment variables

Static clinical features:

- Age at admission; recorded sex plus missing/unknown-sex indicator.
- POA ICD-code proxies for chronic lung disease, uncontrolled diabetes,
  immunosuppression, metastatic solid tumor and hematologic malignancy.
- Glucose and albumin: closest VALID record at or before respiratory-failure
  time zero, within the selected hospital block. They are fixed at baseline.
  A record after time zero never fills a missing baseline value.

Time-varying clinical features considered:

- Temperature, WBC, available pH, arterial PaCO2, SpO2/FiO2 ratio.
- Respiratory support type, represented categorically: room air (reference),
  low-flow oxygen, HFNO, NIPPV/CPAP, IMV. Unknown type is missing, not room air.
- A nonrespiratory SOFA PROXY using platelets, bilirubin, MAP, creatinine and GCS.
  This is not full SOFA: cardiovascular vasopressor and renal urine-output rules
  are unavailable in this construction. All five proxy components must be
  observed (possibly by LOCF); otherwise the score is missing. Missing parts are
  never assumed to contribute zero points. The same 70% exclusion rule applies.

Structural features:

- Hospital type: JHH and BMC academic; HCG, SMH and SUB community (user supplied).
  Unmapped site identities stop analysis rather than default to community.
- Location at time zero: ED, Ward, IMC, ICU or Unknown. The latest ADT interval
  CONTAINING time zero is used. Explicit mappings replace loose substring
  matching; no future or nearest noncontaining location is used as a fallback.
  ED is the reference; unknown is explicit. Review time_zero_location_mapping.csv
  for source labels that need additional mapping.
- CLIF site defines local nuisance models (and pooled site-stratified identity
  control). Hospital type is constant inside each site, so its separate coefficient
  is not identifiable in site-specific models. Constant local terms are dropped
  and the site intercept absorbs structural differences. The type field remains
  in prepared data, the global specification and site-variation audit.
- Calendar year and season (sine/cosine encoding), derived from time zero.

NOT used as additional clinical adjustment variables: heart failure, CKD,
lactate, vasopressor indicators, historical/NCO proxies, weight, race or ethnicity.
Some source fields remain necessary for the original cohort's eligibility checks
(e.g., baseline shock exclusion) or descriptive data-quality audit. Their presence
in an input/cohort file does not mean they enter the adjustment models.

## LOCF, unavailable features and residual missingness

Only valid records at/preceding the six-hour decision are considered. Values
outside plausibility limits are removed BEFORE LOCF, so an invalid latest value
does not hide an older valid reading. Carry-forward is restricted to the selected
hospital block, starts at admission, has no expiry after a valid record, and never
backfills from a future measurement. Therefore a carried-forward observation can
be old; it is not a newly measured physiological state. Models remain vulnerable
to informative testing and unrecorded changes.

For S/F, first pair a SpO2 record with the most recent preceding FiO2 record no
more than four hours earlier. Carry that derived ratio forward; do not divide
independently stale SpO2 and FiO2 series. FiO2 percentages are normalized.
PaCO2 uses arterial categories only; venous PCO2 is not substituted. pH matching
attempts generic/arterial pH categories and excludes explicitly venous categories.
Validate source specimen labels/units; if unavailable, these features are dropped.

The user-selected missingness rule is **exclude if MORE THAN 70% missing**:

- Baseline glucose/albumin: one observation per selected patient.
- Time-varying features: untreated treatment-decision intervals within the grace
  period (including the interval in which initiation occurs). This avoids using
  the full post-treatment/post-discharge follow-up to decide availability.
- Missingness is measured AFTER baseline selection/LOCF and BEFORE imputation.
- Exactly 70% missing is retained; no observed values or a constant processed
  value is excluded. The corresponding missingness term is also excluded when
  the clinical feature is excluded. The rule is pooled across hospitals; inspect
  site-specific data availability before interpreting transportability.

For retained features, residual missingness means no valid prior reading exists.
These residual cells use median single imputation plus an explicit missingness
indicator. Medians/winsorization thresholds come from the appropriate baseline or
pre-initiation decision reference sample. Continuous values are winsorized at
0.5/99.5%; WBC/glucose use log1p; other continuous model variables are standardized.
Respiratory support is categorical and not winsorized as a continuous severity
score. Missingness indicators are processing terms for proposal variables, not
additional clinical risk factors. This is not multiple imputation. If every
acute time-varying feature is excluded, baseline-only estimation is permitted.

## Cohort, exposure and estimators

The copied cohort retains ED admission, POA J12-J18 pneumonia (including J17),
respiratory failure within seven days, baseline shock/code-status exclusions,
transfer merging within six hours and one randomly selected eligible block per
patient. The source's ED/ever-respiratory-support filters remain. Glucocorticoids
retain the original drug-specific per-administration dose/route thresholds and
intermittent MAR source. Prior steroid continuation and continuous steroid capture
are NOT newly addressed in this variable-only analysis; see the separate CAPE COD
folder for the different new-user/hydrocortisone design.

Target strategies are initiate qualifying systemic steroids within grace versus
no qualifying initiation through day 28. CCW, its pooled identity control,
landmark IPW, IPCW without cloning and unweighted per-protocol are reported.
The treatment clone retains observed deaths in the final grace interval and
censors surviving noninitiators after that interval; this corrects the legacy
final-interval boundary. The never-initiate control censors at any initiation,
including after grace. The landmark sensitivity uses baseline-only covariates,
not potentially treatment-affected landmark laboratory measurements; its balance
diagnostic uses the same specification. It remains conditional on survival to
grace end and is not the full-cohort CCW estimand.

Nuisance models are hospital-specific with shared numerator hazards. Constant
local terms are removed. Weights use the existing 1/99% truncation and fixed-weight
influence-function inference; no bootstrap is added. These are still parametric
models, and neither LOCF nor covariate inclusion guarantees confounding control.

## Outputs

- results/<run-id>/run_design.txt and source_versions.csv: settings/checksums.
- logs/: stage logs; potentially disclosive, keep on Discovery.
- sensitive/cohort/, steroid_candidates/, treatment/: patient-level inputs.
- sensitive/analysis/: prepared cohort/person-interval data and analysis objects.
- shareable/cohort_eda/: existing descriptive EDA and six-hour death/treatment plots.
- shareable/ccw_results/: estimates, survival curves, weight/balance diagnostics,
  covariate availability/processing and exact model-term specification.

Review in particular:

- covariate_availability_and_processing.csv: included/dropped features, reason,
  baseline versus time-varying role, missingness denominator, imputation values.
- baseline_feature_variation_by_site.csv: site-constant structural/categorical terms.
- time_zero_location_mapping.csv: observed location vocabulary and mappings.
- matched_lab_categories.csv: actual source pH/arterial PCO2 mappings.
- ccw_analysis_design.csv, ipcw_model_specification.csv, method_warnings.txt.

Patient-level data and logs remain on Discovery. Aggregate exports require the
project's usual disclosure review. EDA still describes the original cohort;
the covariate-processing report is authoritative for model feature availability.

## Tests

```bash
Rscript Proposal_aligned/tests/test_proposal_features.R
Rscript Proposal_aligned/tests/test_ccw_policy.R
Rscript Proposal_aligned/tests/test_end_to_end.R
```

Tests use synthetic data only. The end-to-end test copies the scripts into an
isolated temporary directory, constructs synthetic CLIF inputs and executes every
stage. It verifies missing-pH exclusion, fixed baseline labs, hospital classification,
absence of nonproposal clinical adjustments, and finite estimates. No clinical
results have been generated locally because rclif_2_1 is not available here.

## Resume after a native crash in 07

Vitals/labs now use single-threaded Arrow record batches with category filters
applied before materialization. LOCF/SF matching uses per-block binary searches
instead of non-equi joins and no longer copies the entire person-interval grid
for each variable. Native Arrow/data.table threading is limited to one. Logs
report package versions, table reads, batch progress and individual LOCF stages.
This reduces memory pressure and provides crash localization; it is not proof
of the cause of a native error on a different HPC environment.

After copying the updated folder to Discovery, rerun ONLY analysis using the
existing run directory, for example from CLIF-Cohort:

```bash
Rscript Proposal_aligned/resume_analysis.R \
  Proposal_aligned/results/run_20261106_123152
```

The resumer restores grace, the missingness cutoff and administrative cutoff
from the saved run; cohort and treatment files are reused. It writes a new
`07_resumed_<timestamp>.log`, reruns analysis outputs, records new source checksums
and marks completion only on success. An optional second argument overrides the
EHR source directory if it has moved. Batch size can be lowered with
`PROPOSAL_READ_BATCH_ROWS=16384` (default 65536). Small batches bound scan
materialization; selected records are still retained for analysis.

Additional tests:

```bash
Rscript Proposal_aligned/tests/test_large_alignment.R
Rscript Proposal_aligned/tests/test_batch_reader.R
```
