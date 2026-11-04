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

## Use the current rclif_2_1 extract

The default `HISTORY_MODE=observed_extract` implements option 2: proxies describe
**prior events observed in the filtered source**, not all prior hospitalizations.
Use the same rclif_2_1 directory for both the analytic and historical inputs.
`HISTORY_SOURCE_COMPLETE` is not required in this mode.

`HISTORY_COVERAGE_START` defaults to `2021-01-01`, following the documented project
extraction window; override it if the actual extract starts on a different date.
It describes calendar availability of this filtered extract, not complete patient
history. The minimum observed admission date is not used to infer completeness.

A zero means **no qualifying prior event found in the supplied extract/window**.
The code keeps observed positives even with less than 365 days of calendar
lookback. It adds `history_lookback_days` (capped at 365) and
`prior_history_short_lookback` to baseline adjustment and retains the same cohort
in both runs. No-history observations can reflect filtering or outside-system
care. They must not be interpreted as confirmed absence of prior hospitalization.

## Historical definitions

All four binary indicators refer to separate completed hospital episodes
admitted in the 365 days preceding the index hospital-block admission:

- `observed_prior_icu`: ICU location (`location_category == icu`).
- `observed_prior_respiratory_hospitalization`: any ICD-10-CM J12–J18 or J96 diagnosis,
  or ICD-9-CM 480–486 / 518.81 / 518.82 / 518.84. This historical vulnerability
  definition includes pneumonia and respiratory failure, without a POA restriction.
- `observed_prior_emergency_hospitalization`: an episode containing an ED admission
  (`admission_type_category == ed`).
- `observed_prior_imv`: recorded invasive mechanical ventilation (`device_category == imv`).

Overlapping hospitalizations and gaps of six hours or less are merged using
running maximum discharge times. Episodes connected to the index admission,
future episodes, and episodes admitted before the lookback are excluded.
ICU/IMV timestamps must fall within their component hospitalization.

These are within-system, extract-observed histories. The four proxies may overlap
strongly because the current extract already requires ED admission and respiratory
support. Constant columns are removed by the parent pipeline. They are baseline
vulnerability proxies, not independent negative-control experiments.

Optional future mode: set `HISTORY_MODE=complete_history`, supply unfiltered
historical tables, explicitly set their `HISTORY_COVERAGE_START`, and declare
`HISTORY_SOURCE_COMPLETE=yes`. In that mode, partial-year histories are stored
as NA and zero-filled only for modeling, with the short-lookback indicator. The
runner stops if nobody has a full calendar year. This mode still cannot capture
outside-system care. Missing admission/discharge or admission-type fields and
duplicate encounter IDs are rejected in both modes.

## Run on Discovery

First run existing scripts 01–04, with treatment assignment using the same grace
period. The new runner does not change cohort or treatment inputs, or rerun EDA.
Use the current filtered extract in both directory arguments below. The date
below follows the documented extraction window; change it if needed.

```bash
cd /projects/emulate_icu/CLIF-Cohort
export TREATMENT_GRACE_HOURS=96
export HISTORY_MODE=observed_extract
export HISTORY_COVERAGE_START=2021-01-01
export SUMMARY_SMALL_CELL_THRESHOLD=10
export CLIF_DATA_CUTOFF_DATE=2026-10-31
export CCW_MAX_PATIENTS=0

Rscript historical_proxy_sensitivity/run_analysis.R \
  /projects/emulate_icu/rclif_2_1 \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_treatment \
  /projects/emulate_icu/rclif_2_1 \
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

- `sensitive/historical_proxies.parquet`: extract-observed historical indicators, source mode and
  calendar lookback length, one row per selected cohort block. Keep on Discovery.
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
Rscript historical_proxy_sensitivity/test_runner_integration.R
```
