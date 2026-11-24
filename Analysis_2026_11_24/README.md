# CLIF pneumonia and glucocorticoid analysis

This folder contains a standalone analysis of early systemic glucocorticoid
initiation in adults with pneumonia and respiratory failure. It builds the
cohort from CLIF tables, prepares longitudinal covariates, estimates treatment
contrasts, and produces diagnostic summaries. This version includes scripts 01–08 through local CCW/meta-analysis.

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
| Follow-up | Up to 28 days, in configurable intervals (`HOSPITAL_BLOCK_HOURS`, default 6) |
| Primary outcome | In-hospital death or hospice discharge by day 28 in the selected block |
| Unit | One randomly selected eligible hospitalization block/episode per patient |
| Sites | Hospital IDs resolved from source records |

Overlapping encounters and encounters separated by at most `HOSPITAL_BLOCK_GAP_HOURS`
(default 6 hours) are merged into an episode that can last days or weeks.
`HOSPITAL_BLOCK_HOURS` separately controls the analysis-interval length. 


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
cd /projects/emulate_icu/CLIF-Cohort/Analysis_2026_11_24
bash RUN_PIPELINE_COMMANDS.txt
```

Before launching, review the exports in that command file:

| Setting | Command-file value |
|---|---|
| `CLIF_DATA_DIR` | `/projects/emulate_icu/rclif_2_1` via a relative-path calculation |
| `TREATMENT_GRACE_HOURS` | `24` |
| `CLIF_DATA_CUTOFF_DATE` | `2026-10-31` |
| `HOSPITAL_BLOCK_HOURS` | `6` (analysis interval length) |
| `HOSPITAL_BLOCK_GAP_HOURS` | `6` (encounter-merging gap) |
| `COHORT_SELECTION_SEED` | `2026` |
| `PNEUMONIA_CASE_DEFINITION` | `poa_any` |
| `SUMMARY_SMALL_CELL_THRESHOLD` | `1` |

Set values when launching the bash command file, or edit its defaults:

```bash
HOSPITAL_BLOCK_HOURS=3 bash RUN_PIPELINE_COMMANDS.txt
HOSPITAL_BLOCK_HOURS=12 bash RUN_PIPELINE_COMMANDS.txt
TREATMENT_GRACE_HOURS=96 HOSPITAL_BLOCK_HOURS=12 bash RUN_PIPELINE_COMMANDS.txt
# Separate sensitivity for merging hospital encounters:
HOSPITAL_BLOCK_GAP_HOURS=12 bash RUN_PIPELINE_COMMANDS.txt
```

`HOSPITAL_BLOCK_HOURS` must be a positive whole number of hours dividing both
`TREATMENT_GRACE_HOURS` and the 672-hour follow-up exactly. Grace must be at least
one interval and below 672 hours. With 24-hour grace, examples include 1, 2, 3,
4, 6, 8, 12, and 24 hours. The 3/6/12-hour settings yield 8/4/2 grace intervals
and 224/112/56 follow-up intervals. Invalid timing stops before cohort construction.
The run folder records grace, interval, and encounter-gap settings.

Changing the interval rebuilds decision-time covariates, event/treatment/censoring
intervals, and weights; covariate lags remain one analysis interval. Clinical
lookbacks (including four-hour S/F pairing and six-hour shock screening), treatment
thresholds, follow-up, and weight truncation are unchanged. The descriptive raw-event
timing histograms and `outcome_and_treatment_six_hour_counts.csv` deliberately retain
six-hour bins and the separate 24/96-hour comparisons; analysis-grid covariate EDA
uses the configured interval. Script 08 reads timing from the completed 07 design.
Each command-file invocation runs the whole 01–08 pipeline into a fresh folder.

We also now have `/projects/emulate_icu/rclif_2_1/clif_hospitalization.parquet` with expanded encounters before and after potential index encounters. 


## Script 

| Script | Purpose | Main output under the result folder |
|---|---|---|
| 01 | Construct eligible blocks, select one per patient | `sensitive/cohort/` |
| 02 | Check steroid labels, drugs, routes, units, and administration | `sensitive/source_profile/` |
| 03 | Extract candidate steroid rows linked to the cohort | `sensitive/steroid_candidates/` |
| 04 | Apply treatment rules and determine treatment timing | `sensitive/steroid_treatment/` |
| 05 | Cohort EDA, treatment/event timing, support trajectories, covariate availability | `aggregate/eda/` |
| 06 | Shared method library; sourced in 07 | Used by 07–08 |
| 07 | Prepare analysis data and run five methods with balance/weight diagnostics | `aggregate/analysis/` |
| 08 | Fit CCW independently within each site, then meta-analyze site contrasts | `aggregate/local_ccw_meta/` |

07 reports federated CCW, pooled site-stratified CCW, landmark IPW, IPCW without
cloning, and unweighted per-protocol estimates. Outputs include risk differences,
risk ratios, odds ratios, and restricted mean event-free time differences.
Initiation/adherence weights freeze after grace; natural-censoring weights can
continue when applicable. Weight truncation uses the 1st/99th percentiles.
Landmark analysis conditions on reaching the end of grace without an event.

08 adds local CCW plus REML random-effects meta-analysis and writes site forest plots and comparisons of all six methods. 



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
results/run_<timestamp>_grace24h_block6h_gap6h/
  aggregate/
    eda/
    analysis/
    local_ccw_meta/
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
