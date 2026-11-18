# CLIF pneumonia and glucocorticoid analysis

This folder contains a standalone analysis of early systemic glucocorticoid
initiation in adults with pneumonia and respiratory failure. It builds the
cohort from CLIF tables, prepares longitudinal covariates, estimates treatment
contrasts, and produces diagnostic summaries. Historical negative-control
outcome (NCO) analyses are optional.

**Start here:** review the study definitions below, then scripts **01 and 04**
for cohort and treatment construction. For methods, read **07 and 08** together
with the helpers they source. Use `RUN_PIPELINE_COMMANDS.txt` to run the pipeline.
No scripts from the parent folder are required.

## Study question and definitions

| Item | Current implementation |
|---|---|
| Comparison | Initiation within the grace period versus no initiation within that period; comparator patients may initiate later or never initiate|
| Grace period | 24 hours (by default) |
| Time zero | Qualifying respiratory-support initiation in the selected hospitalization block |
| Follow-up | Up to 28 days, represented in six-hour intervals |
| Primary outcome | In-hospital death or hospice discharge by day 28 in the selected block |
| Unit | One randomly selected eligible hospitalization block/episode per patient |
| Sites | Hospital IDs resolved from source records |

Overlapping encounters and encounters separated by no more than six hours are merged into an episode (that
can last days or weeks). 
Six hours is also the analysis-interval length. 


### Cohort eligibility: 01

The principal criteria are:

- Age at admission at least 18 and below 120 years.
- ED-origin admission (`admission_type_category = ed`).
- Qualifying ICD-10-CM pneumonia present on admission (poa), in any diagnosis
  position (`poa_any`): J12–J18, including J17. Influenza and the separately
  specified excluded pneumonia codes do not qualify; see `pneumonia_icd_helpers.R`.
- IMV or NIPPV initiation, or HFNO initiation with S/F <315 and FiO2 >0.50,
  within the first seven days of the block. HFNO qualification requires
  SpO2 ≤97%; there is no minimum future support-duration requirement though.
- No recorded death before time zero and no qualifying DNR/DNI-equivalent status.
- No concurrent shock proxy under the implemented rules: norepinephrine or
  epinephrine on NIPPV/HFNO; dose >0.25 mcg/kg/min on IMV; any vasopressin;
  or lactate ≥4 mmol/L. Pressor and lactate lookbacks are six hours.

Currently we have around ~4,500 patients. 
Review `adult_pneumonia_attrition.csv` and `adult_pneumonia_cohort_design.csv`
for the actual run. 

### Treatment definition: 03–04

Source: `clif_medication_admin_intermittent.parquet`, with `med_group = steroid`.
03 extracts IV/enteral candidates; 04 applies the following thresholds **per
administration**, not per day:

| Drug | Minimum dose | Eligible routes |
|---|---:|---|
| Dexamethasone | 6 mg | IV or enteral |
| Hydrocortisone | 50 mg | IV |
| Prednisone | 40 mg | Enteral |
| Methylprednisolone | 32 mg | IV or enteral |
| Prednisolone | 40 mg | IV or enteral |

**Only `mar_action_category = given` can establish treatment.** `not_given`
and unknown actions are excluded. Other drugs/routes and unsupported dose
units, including `mg/kg/dose`, do not qualify. 


### Outcome and discharge handling

`outcome_helpers.R` defines the endpoint consistently for the analyses:

- Death timestamp within the selected hospitalization block; an expired
  discharge supplies a fallback timestamp when needed.
- Hospice discharge at its discharge timestamp.
- The earliest qualifying component determines the event.
- Other alive discharges, including transfers, are not primary events and are
  not censoring at discharge; they remain event-free through day 28 unless administratively censored.
- Post-discharge death timestamps are excluded from this endpoint.

Thus the endpoint differs from 28-day all-cause mortality. Administrative
censoring uses the configured data cutoff. 
Again see `outcome_ascertainment_audit.csv`and `natural_censoring_counts.csv` for specific counts.

## Run on the data server

I prepared `RUN_PIPELINE_COMMANDS.txt` to run the full pipeline in one command. 
It is a bash script that sets up the environment and runs each numbered script in order. 

```bash
cd /projects/emulate_icu/CLIF-Cohort/Analysis_2026_11_16
bash RUN_PIPELINE_COMMANDS.txt
```

Before launching, review the exports in that command file:

| Setting | Command-file value |
|---|---|
| `CLIF_DATA_DIR` | `/projects/emulate_icu/rclif_2_1` via a relative-path calculation |
| `TREATMENT_GRACE_HOURS` | `24` |
| `CLIF_DATA_CUTOFF_DATE` | `2026-10-31` |
| `HOSPITAL_BLOCK_GAP_HOURS` | `6` |
| `COHORT_SELECTION_SEED` | `2026` |
| `PNEUMONIA_CASE_DEFINITION` | `poa_any` |
| `SUMMARY_SMALL_CELL_THRESHOLD` | `1` |

Edit the command-file grace-period export to try another grace period window. 
The timing EDA separately compares both 24 and 96 hours. 
The command file creates a fresh `RUN_ROOT` for each run.

We also now have `/projects/emulate_icu/rclif_2_1/clif_hospitalization.parquet` with expanded encounters before and after potential index encounters. 


## Script 

| Script | Purpose | Main output under the result folder |
|---|---|---|
| 01 | Construct eligible blocks, select one per patient | `sensitive/cohort/` |
| 02 | Check steroid labels, drugs, routes, units, and administration | `sensitive/source_profile/` |
| 03 | Extract candidate steroid rows linked to the cohort | `sensitive/steroid_candidates/` |
| 04 | Apply treatment rules and determine treatment timing | `sensitive/steroid_treatment/` |
| 05 | Cohort EDA, treatment/event timing, support trajectories, covariate availability | `aggregate/eda/` |
| 06 | Shared method library; sourced in 07 | Used by 07–10 |
| 07 | Prepare analysis data and run five methods with balance/weight diagnostics | `aggregate/analysis/` |
| 08 | Fit CCW independently within each site, then meta-analyze site contrasts | `aggregate/local_ccw_meta/` |
| 09, optional | Diagnose any prior-year hospitalization as a historical NCO | `aggregate/nco_prior_hospitalization/` |
| 10, optional | Build a historical NCO panel and explore NN/MNN bias calibration | `aggregate/nco_calibration/` |

07 reports federated CCW, pooled site-stratified CCW, landmark IPW, IPCW without
cloning, and unweighted per-protocol estimates. Outputs include risk differences,
risk ratios, odds ratios, and restricted mean event-free time differences.
Initiation/adherence weights freeze after grace; natural-censoring weights can
continue when applicable. Weight truncation uses the 1st/99th percentiles.
Landmark analysis conditions on reaching the end of grace without an event.

08 adds local CCW plus REML random-effects meta-analysis and writes site forest plots and comparisons of all six methods. 

09 and 10 are work in progress for me to explore NCOs.

## Covariates and measurement processing

`analysis_preparation.R` builds the analysis inputs. Baseline candidates include
age, sex, initial support/location, admission year/season, and coded
comorbidities. Time-varying candidates include temperature, WBC, glucose,
albumin, a non-respiratory SOFA proxy, categorical respiratory support, and categorical S/F.
The initiation denominator includes current and lagged time-varying terms.
Arterial pH and PaCO2 are excluded (too much missing). Constant or unavailable terms can be removed;
**`ipcw_model_specification.csv` lists the exact terms used in each run.**

### Screening

`variable_rules.R` contains the inclusive screening ranges:

| Variable | Range |
|---|---|
| Temperature | 30–45 °C |
| SpO2 | 70–100%; S/F additionally requires ≤97% |
| FiO2 | 0.21–1; documented percentage values are converted to fractions |
| S/F | 70–462 |
| MAP | 0–250 mmHg |
| WBC | 0–500 ×10³/µL |
| Glucose | 0–2000 mg/dL |
| Albumin | 0–15 g/dL |
| Platelets | 0–2000 ×10³/µL |
| Total bilirubin | 0–80 mg/dL |
| Creatinine | 0–30 mg/dL |
| GCS | 3–15 |

Outside-range measurements become missing; a patient is not removed solely for one invalid measurement. 
WBC and glucose use log(1 + value) transformations.

### Missingness and temporal selection

For non-S/F inputs, use a prior value within the variable-specific lookback,
then the nearest valid measurement in the same block to fill remaining inpatient
gaps. Ties favor earlier records; the fallback has no maximum distance. Source
records must be within the block and no later than the primary event or data
cutoff. **This can use future/post-treatment measurements.** Remaining numeric
gaps receive median imputation and missingness indicators. 
The non-respiratory SOFA proxy requires all five components and is not the complete SOFA score.

Respiratory support is represented as room air, conventional oxygen, HFNO, NIPPV/CPAP, IMV, or Unknown.


### S/F pairing

For analysis covariates, pair SpO2 with prior FiO2 documented within four hours
on the same continuous device episode. Device labels are forward-filled only;
a device change starts a new episode. At each decision time, use the latest
prior SpO2 no more than four hours old; an intervening device change invalidates
the pair. SpO2 >97% remains Undefined. 
FiO2 age is checked at SpO2 time, and SpO2 age is checked separately at decision time.

S/F categories are Undefined, >300 (reference), (250,300], (200,250],
(120,200], and ≤120, entered as separate indicators with their lags. S/F is not
median-imputed. The **01 HFNO eligibility window is separate**: it uses a
symmetric ±4-hour clinical qualification window, which can include later values.

## Where to find results

```text
results/run_<timestamp>_grace24h/
  aggregate/
    eda/
    analysis/
    local_ccw_meta/
    nco_prior_hospitalization/   # optional
    nco_calibration/            # optional
  sensitive/                   # patient-level inputs, histories, fits (not shared)
  logs/                        # one log per executed script
```

| FAQ | Files to inspect |
|---|---|
| Who entered the cohort? | `sensitive/cohort/adult_pneumonia_attrition.csv` and `adult_pneumonia_cohort_design.csv`; EDA also exports cohort summaries |
| Which administrations counted? | `sensitive/steroid_treatment/steroid_dose_route_rules.csv`, `steroid_administration_action_audit.csv`, and row-level eligibility parquet |
| What were the treatment estimates? | `aggregate/analysis/five_method_results.csv`; `aggregate/local_ccw_meta/six_method_results.csv` |
| Which confounders were fitted? | `aggregate/analysis/ipcw_model_specification.csv` |
| Where is balance poor? | `covariate_balance_by_method_and_interval.csv` and love-plot data in `aggregate/analysis/` |
| Were weights extreme? | Weight summaries, clipping fractions, ESS, and fit diagnostics in `aggregate/analysis/` |
| How much information was filled or missing? | `model_relevant_missingness.csv`, `adjacent_fill_audit_*.csv`, and `sf_*coverage.csv` in `aggregate/analysis/` |
| What happened within individual sites? | `aggregate/local_ccw_meta/local_ccw_site_results.csv` and `meta_analysis_results.csv` |

Aggregate figures/counts are unsuppressed by default, including cells below 10.
Patient-level files belong in `sensitive/` and should remain on the data server.
Review output contents before sharing; an aggregate folder name is not an
independent disclosure review. Sample sizes and estimates are run-specific and
are not hardcoded in this README.

## Optional historical NCO analyses

To run 09 and 10 on an existing completed run:

```bash
cd /projects/emulate_icu/CLIF-Cohort/Analysis_2026_11_16
export RUN_ROOT="$PWD/results/REPLACE_WITH_COMPLETED_RUN_FOLDER"
export EXPANDED_HOSPITALIZATION_PATH=/projects/emulate_icu/rclif_2_1/clif_hospitalization.parquet
mkdir -p "$RUN_ROOT/logs"
set -euo pipefail

Rscript 09_prior_hospitalization_nco_diagnostic.R \
  "$RUN_ROOT" "$EXPANDED_HOSPITALIZATION_PATH" 2>&1 | tee "$RUN_ROOT/logs/09.log"
Rscript 10_historical_nco_calibration.R \
  "$RUN_ROOT" "$EXPANDED_HOSPITALIZATION_PATH" 2>&1 | tee "$RUN_ROOT/logs/10.log"
```

Both accept an optional third argument: a verified **system-level** history
coverage-start date (`YYYY-MM-DD`). Without it, historical capture is explicitly
unverified. Do not use a patient's first observed encounter as that date.
To include 09/10 automatically in a full pipeline run, export
`EXPANDED_HOSPITALIZATION_PATH` before launching the command file; optionally
export `HISTORY_COVERAGE_START_DATE` if known.

09 measures any completed within-system episode starting in the 365 days before
index block admission. It merges connected encounters, excludes the entire
index-connected episode, and flags uncertain dates. It reuses saved 07 weights
and clone risk sets to compare historical prevalence during grace, plus
observed-group and landmark comparisons. It does not replace the primary
survival outcome or refit the clinical analysis.

10 adds recurrent admissions (≥2), stays ≥7 days, ED-origin admissions, and
SNF/rehab discharges. Missing historical category fields remain unknown; ICU
history is not inferred without historical location data. It fits NN
(Normal–Normal) and MNN (mixture Normal–Normal) systematic-error models to
final-grace-interval clipped-weight log prevalence ratios. MNN designates the
majority component as valid by assumption, not by proximity to zero.

Main outputs in `aggregate/nco_calibration/` are:

- `nco_calibration_inputs.csv`: estimates, corrected log-ratio SEs, and exclusions.
- `nn_mnn_fit_diagnostics.csv`, `nn_mnn_bias_distributions.png`: fits and densities.
- `nco_patient_level_correlations.csv`, `nco_mnn_membership.csv`,
  `leave_one_control_out.csv`: overlap, membership, and stability checks.
- `exploratory_calibrated_ccw_rr.csv`: primary log RR minus the fitted bias mean,
  exponentiated to RR, with plug-in sensitivity intervals.

**Calibration is exploratory.** These controls are correlated utilization
measures, and risk-set selection can generate associations without residual
confounding. Transfer to mortality bias is unverified; `calibration_supported`
therefore remains FALSE. Intervals omit calibration-parameter uncertainty,
shared-patient dependence, nuisance-model uncertainty, and transport uncertainty.
There is no formal calibrated p-value or RD/RMST correction. NN needs at least
three usable controls and exploratory MNN four; these are computational guards,
not evidence of an adequate panel. `NCO_MIN_CELL_N` defaults to 10 positive and
10 negative observations in each arm. Failed MNN fits are reported, not replaced
with NN under the MNN label.


Other tests cover hospital blocks, outcomes, adjacent filling, S/F pairing, and
local meta-analysis. Some accept a completed synthetic run for integration
checks; read their usage comments before supplying an argument.
