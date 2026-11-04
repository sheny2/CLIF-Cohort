# Historical-proxy sensitivity analysis (approach A)

This separate runner reuses the parent five-method pipeline without editing it.
It evaluates the parent's parsed R expressions in isolated environments, sets
its grace period from TREATMENT_GRACE_HOURS (default 96), and adds historical
proxies only in the augmented run. It records the parent source MD5 checksums.
It fails if expected insertion points are absent. Parent pipeline updates should
be reviewed for compatibility before subsequent runs.

This is measured baseline proxy adjustment, not empirical NN/MNN calibration
or a guarantee of eliminating unmeasured confounding. No NCO association is
subtracted from the mortality estimate.

## Historical definitions

All four binary indicators refer to separate completed hospital episodes
admitted in the 365 days preceding the index hospital-block admission:

- `prior_icu`: ICU location (`location_category == icu`).
- `prior_respiratory_hospitalization`: any ICD-10-CM J12–J18 or J96 diagnosis,
  or ICD-9-CM 480–486 / 518.81 / 518.82 / 518.84. This historical vulnerability
  definition includes pneumonia and respiratory failure, without a POA restriction.
- `prior_emergency_hospitalization`: an episode containing an ED admission
  (`admission_type_category == ed`).
- `prior_imv`: recorded invasive mechanical ventilation (`device_category == imv`).

Overlapping hospitalizations and gaps of six hours or less are merged using
running maximum discharge times. Episodes connected to the index admission,
future episodes, and episodes admitted before the lookback are excluded.
ICU/IMV timestamps must fall within their component hospitalization.

These are within-system histories. Hospitalization tables must include all
admissions, with linked diagnosis, ADT, and respiratory tables; a table filtered
to ED admissions or respiratory support cannot establish historical absence.
The runner requires an explicit source-completeness declaration and historical
coverage start. This is a data contract, not an inferred coverage test. Calendar
coverage does not establish that a patient always received care in this system.

For index dates without a full calendar year of source coverage, indicators are
stored as NA, then zero-filled for modeling with `prior_history_unavailable=1`.
All patients remain in both runs; the effects target the same analytic cohort.
This avoids coding unknown history as confirmed absence, but does not reconstruct
missing historical information. The runner stops if nobody has full coverage.
Missing admission/discharge or admission-type fields and duplicate encounter IDs
are rejected. Verify historical category mappings before interpreting results.

## Run on Discovery

First run existing scripts 01–04, with treatment assignment using the same grace
period. The new runner does not change cohort or treatment inputs, or rerun EDA.
Use an unfiltered historical CLIF directory, potentially distinct from the
analytic extract. Edit the historical source and coverage date below to actual
values; no historical date or completeness is assumed by the script.

```bash
cd /projects/emulate_icu/CLIF-Cohort
export TREATMENT_GRACE_HOURS=96
export HISTORY_COVERAGE_START=2020-01-01  # EXAMPLE: replace with actual coverage
export HISTORY_SOURCE_COMPLETE=yes      # only when the source meets the contract
export SUMMARY_SMALL_CELL_THRESHOLD=10
export CLIF_DATA_CUTOFF_DATE=2026-10-31
export CCW_MAX_PATIENTS=0

Rscript historical_proxy_sensitivity/run_analysis.R \
  /projects/emulate_icu/rclif_2_1 \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_treatment \
  /path/to/unfiltered_historical_clif \
  ./historical_proxy_sensitivity/output_96h_run1
```

The output root must be new; existing roots are refused. Both reference and
augmented runs use identical grace periods, cohorts, sampling settings,
covariate processing, and follow-up definitions. Constant proxies are dropped
by the parent's existing model availability checks. Proxies are baseline terms,
so they enter the treatment-initiation denominator and other models that use the
baseline list, including the censoring denominator and landmark model. The
initiation numerator specification remains the same.

## Outputs and interpretation

- `sensitive/historical_proxies.parquet`: original historical indicators and
  availability, one row per selected cohort block. Keep on Discovery.
- `sensitive/reference/` and `sensitive/augmented/`: prepared individual data and
  run provenance; keep on Discovery.
- `shareable/reference/` and `shareable/augmented/`: the usual aggregate results,
  model specifications, weight diagnostics, balance plots and processing summaries.
- `shareable/reference_vs_augmented.csv`: estimates, SEs and CIs side by side.
- `shareable/historical_proxy_summary.csv`: suppressed historical frequencies.

Review actual model terms, overlap/weight warnings, and proxy balance before
interpreting effect changes. The supplied uncertainty is the existing pipeline's
uncertainty; the comparison does not supply a paired bootstrap CI for the change.
The four indicators overlap and must not be interpreted as independent controls.
All aggregate exports remain subject to the project's usual disclosure review.

No patient data were used in local verification. Run the synthetic feature tests:

```bash
Rscript historical_proxy_sensitivity/test_history_helpers.R
```
