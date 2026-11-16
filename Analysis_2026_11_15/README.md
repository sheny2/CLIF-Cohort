# Standalone CLIF analysis — November 15, 2026

This folder contains a complete 01–07 pipeline. It does not source any script from the parent CLIF-Cohort folder and does not change the original, CAPE COD-aligned, or proposal-aligned analyses. Raw CLIF data remain in `../rclif_2_1` relative to CLIF-Cohort; from this analysis folder the source data path is `../../rclif_2_1`. The supplied commands use the renamed `Analysis_2026_11_15` folder. All results go into a new timestamped directory under this folder's `results/`, containing only `aggregate/`, `sensitive/`, and `logs/` at the top level.

## Run

Install R packages `arrow`, `data.table`, `dplyr`, `ggplot2`, and `survival`. Open `RUN_PIPELINE_COMMANDS.txt`, adjust the first server path and administrative cutoff if necessary, and execute its commands in bash/zsh. Script 06 is a library automatically sourced by 07. There is no Cox analysis in this version.

The default grace period is **24 hours**, with six-hour intervals and 28-day follow-up for in-hospital death or hospice discharge. Set `TREATMENT_GRACE_HOURS` to another positive multiple of six below 672 to change it. Results record the actual setting; figure names describe the grace period without assuming a fixed length. Timing EDA compares both 24 and 96 hours.

## Scripts

| Script | Purpose |
|---|---|
| 01 | Construct hospital blocks and the eligible adult ED-admitted pneumonia cohort; one eligible block per patient, seed 32284 |
| 02 | Inventory glucocorticoid source fields and routes |
| 03 | Extract steroid candidates for the selected blocks |
| 04 | Apply existing systemic drug/dose/route rules; retain qualifying administrations throughout follow-up |
| 05 | Cohort and treatment/primary-event timing summaries; screened measurement ranges, missingness, temperature-fill audit, categorical support proportions and transitions |
| 06 | Five-method library; early versus later/never initiation, with raw and clipped weights retained for diagnostics |
| 07 | Independently prepare covariates, fit all five methods, produce estimates, balance plots, explicit model specifications and weight diagnostics |

Helpers in this folder provide ICD definitions, follow-up conventions, measurement rules, shared preparation, EDA, and weight diagnostics. Keep them with the numbered scripts when copying this folder.

## Cohort and strategies

The cohort retains the original main-pipeline rules: completed hospital blocks joining overlapping encounters or gaps of at most six hours; age 18–119; ED admission; present-on-admission pneumonia J12–J18 including J17; exclusion of influenza, aspiration and ventilator-associated pneumonia; qualifying respiratory-support initiation within seven days of admission; no recorded death before time zero; and no DNR/DNI or specified concurrent shock proxy. The original norepinephrine/epinephrine threshold remains >0.25 mcg/kg/min on IMV; other shock criteria are unchanged. This is not a new CAPE COD-specific eligibility redesign.

**S/F eligibility changes:** HFNO qualification requires FiO2 >0.50, S/F <315, and SpO2 **80–97%**, retaining the original lower eligibility bound and pairing windows. This can change cohort membership/time zero compared with earlier runs. Rebuild 01 rather than reusing an old cohort. Every analytic S/F covariate also requires the selected SpO2 to be **50–97%** and FiO2 0.21–1. SpO2 98–100% remains a valid saturation reading but supplies no eligible S/F ratio; it is never capped at 97. Inputs are the latest valid prior readings within 12 hours, selected independently as in the original analysis. Their time separation can therefore be as large as 12 hours.

The two strategies are:

1. Initiate a qualifying systemic glucocorticoid within the configured grace period.
2. Do not initiate within the grace period; initiation later in follow-up is allowed.

Only initiation during grace artificially censors a control clone. Later initiation remains in the comparison strategy. Initiation denominator models use untreated decision rows **during grace only**. Treatment/adherence weights accumulate through grace, then freeze; administrative-censoring weights can continue if applicable. The stabilized initiation numerator is the pooled interval-specific initiation rate, broadcast to each site's model. Hospital-specific logistic denominator models and 1st/99th percentile site/strategy weight clipping are retained.

The methods are federated CCW, pooled site-stratified CCW, landmark IPW, IPCW without cloning, and unweighted per-protocol. Federated/pooled CCW are numerically identical by construction. The latter two observed-group analyses remain exploratory and vulnerable to treatment-timing bias. Landmark restricts to patients alive/observed at grace end and uses the inherited last grace-interval measurements; these may already be affected by early treatment.

The primary outcome throughout EDA and all five methods is **in-hospital death or hospice discharge by day 28**, measured from respiratory-failure start in the selected hospitalization block:

- A `death_dttm` within that block supplies the in-hospital death time. If discharge is `Expired` and no usable in-block death timestamp exists, use discharge time.
- A discharge category containing `hospice` supplies a hospice event at discharge time. If both components occur, use the earlier event, counting the patient once.
- Death after another alive discharge is excluded, even when a later death timestamp is recorded. Other alive discharges remain non-events for this endpoint; discharge is not censoring. Subsequent separate readmissions are not added to the selected block.
- Events after day 28 do not count. Follow-up is administratively censored at the supplied data cutoff when needed.

`outcome_definition.csv` and `outcome_ascertainment_audit.csv` document these rules, component counts, and ignored post-discharge death timestamps. Plot/table labels describe the composite outcome. The former survival curves now show **freedom from the primary outcome**, and RMST means **restricted mean event-free time**, not all-cause survival time.

Six-hour event/initiation discretization and clone boundary conventions are inherited, including censoring an untreated-by-deadline treatment clone before the final grace interval. Composite events during grace are not removed from the baseline cohort.

## Measurement screening and missingness

`variable_rules.R` is the single editable table of inclusive analytic screening ranges. These are broad plausibility screens, not clinical normal ranges and not proof that all excluded extremes are source errors. Outside-range/nonfinite measurements become missing **before** temporal selection; the patient/block is retained. Thus a zero temperature cannot override an older valid temperature. Raw sources are never overwritten. EDA saves `screening_ranges.csv`, `raw_measurement_screening_counts.csv`, and cleaned min/max/percentile summaries.

| Variable | Inclusive range |
|---|---|
| Temperature | 25–45 C |
| SpO2 | 50–100%; S/F additionally requires <=97% |
| MAP | 20–200 mmHg |
| WBC | 0.1–200 x10^3/uL |
| Serum glucose | 20–1000 mg/dL |
| Albumin | 0.5–7 g/dL |
| Platelets | 0–2000 x10^3/uL |
| Total bilirubin | 0–100 mg/dL |
| Creatinine | 0–30 mg/dL |
| GCS | 3–15 |
| FiO2 | 0.21–1; explicit percent-form values 21–100 are divided by 100 |
| S/F | 20–700, subject to the input requirements above |
| Non-respiratory SOFA proxy | 0–17 |

No automatic lab-unit conversion is inferred from extreme values. Lab aliases explicitly distinguish serum/blood measurements. **Arterial pH and arterial PaCO2 are excluded a priori**, including their missingness indicators and lags. Venous gases do not replace them.

Temperature first uses the latest valid prior reading within 12 hours. For remaining gaps **while the patient is within the selected hospital block**, use the nearest valid reading from the same block in either temporal direction; ties favor the earlier reading. No maximum fallback distance is imposed. Candidate readings must be within block admission/discharge, no later than the primary event or administrative cutoff. No other patient or block supplies a temperature, and post-discharge rows do not receive new nearest-temperature fills. If a block has no valid temperature, its gap remains missing and is handled with the numeric imputation below.

**Nearest temperature can use future/post-treatment measurements.** This is the requested adjacent-reading rule, not a strictly prospective covariate definition. `temperature_nearest_fill_audit.csv` reports counts of future fills, remaining gaps, distances and fills over 24 hours. The detailed row audit is sensitive. This rule may introduce bias; review those counts before interpreting estimates.

Other numeric covariates retain bounded prior-value carry-forward: vitals/support/FiO2 12 hours, WBC/glucose/GCS 24 hours, albumin/platelets/bilirubin/creatinine 48 hours. For numeric model values, remaining gaps receive the median and a missingness indicator. Winsorization at 0.5/99.5 percentiles, medians, and scaling are estimated from untreated decision rows during grace, avoiding domination by post-discharge rows. Raw screened values remain in the sensitive long table. Constant/unavailable model terms are omitted and the actual included terms are exported.

The non-respiratory SOFA **proxy** uses platelets, bilirubin, MAP, creatinine and GCS; it lacks vasopressor and urine-output components. Unlike the older preparation, a missing component is not silently scored as normal: all five must be available for the raw sum, otherwise the proxy is missing and receives its own imputation/indicator if usable. Diagnosis-based absence remains coded as no recorded comorbidity; sex/location coding follows the original primary implementation.

## Categorical respiratory support

| Category | Raw device categories |
|---|---|
| Room air | room_air |
| Conventional oxygen | nasal_cannula, face_mask, trach_collar |
| HFNO | high_flow_nc |
| NIPPV/CPAP | nippv, cpap |
| IMV | imv |
| Unknown | No mapped valid support observation within the prior 12 hours |

The model uses separate binary indicators for conventional oxygen, HFNO, NIPPV/CPAP, IMV and unknown; room air is the reference. Constant indicators are omitted. No numeric ordinal support score is entered. Baseline eligibility support retains separate IMV and NIPPV/CPAP indicators with HFNO as reference. Within an identical timestamp, the highest mapped intensity is selected deterministically; this ordering is only a tie-break and is not a numeric model effect. Unmapped/blank source rows cannot supply a known category; the latest mapped prior category may still be carried forward for up to 12 hours.

EDA counts **one patient per six-hour time point**. Full-cohort proportions include separate in-hospital-death, hospice-discharge, other-alive-discharge and administrative-censoring states; a second summary conditions on being alive and within the selected hospital block and includes unknown documentation. Other alive-discharge status means the patient has left this selected block without the primary event; later readmission outside the block is not assigned a support mode. Transitions summarize consecutive six-hour states, not every change between raw chart entries. Neither plot is a causal treatment-effect estimate.

## IPCW covariates and diagnostics

Candidate baseline variables: age, female/sex missingness, baseline support categories, initial ICU/ED location, admission year and season, chronic lung disease, uncontrolled diabetes, immunosuppression, metastatic solid tumor, hematologic malignancy, heart failure, and chronic kidney disease. Site defines separate nuisance models rather than another patient-level covariate.

Candidate current variables: temperature, WBC, serum glucose, albumin, eligible S/F, non-respiratory SOFA proxy, categorical support indicators, and applicable numeric missingness indicators. The initiation denominator includes current and one-interval lagged versions plus baseline terms and interval time. Initial lags are zero, following the original convention. The natural-censoring denominator includes time, baseline terms, and treatment history; its numerator includes time. The actual post-availability-check list is **`ipcw_model_specification.csv`**; `covariate_availability_and_processing.csv` records transformations.

Key outputs under `aggregate/analysis/`:

- `five_method_results.csv`: risks, RD, RR, OR, RMST and confidence intervals for all five methods; effect and event-free probability figures.
- `ipcw_model_specification.csv`: exact terms grouped as baseline, current, lagged, time and treatment.
- `covariate_balance_by_method_and_interval.csv` and Love plots: worst absolute SMD during grace; a single landmark assessment.
- `ccw_weight_summary_and_ess_by_interval.csv` and `ccw_pooled_ess_by_interval.csv`: raw/clipped tails, counts above 10/20/50/100, top-1% weight mass, ESS and ESS/N in the exact risk sets.
- `initiation_probability_summary.csv` and probability histograms: conditional six-hour initiation probabilities among untreated grace-period decisions, not one baseline PS.
- `weight_model_coefficients.csv`, `weight_model_fit_diagnostics.csv`, `method_warnings.txt`: large coefficients, rank/aliasing, convergence and warnings. Estimates are not automatically withheld or repaired when a fit warns; inspect these files before inference.
- `extreme_weight_covariate_profiles.csv`: encoded covariate means in the top raw-weight tail versus other person-intervals.
- `landmark_weight_summary_and_ess.csv`: raw/clipped landmark weight distribution and ESS.

Confidence intervals retain the existing influence-function approach treating estimated nuisance models as fixed; they do not propagate all uncertainty from fitting weights or temperature filling. Improved balance does not establish absence of unmeasured confounding.

## Data handling and validation

All EDA figures and aggregate tables are unsuppressed by default, including counts below 10. Patient-level data, individual extremes, temperature-fill rows, model objects and weight trajectories go only to sensitive directories. Aggregate output still needs the project's normal disclosure review before export.

Run the included synthetic integration check with `Rscript tests/test_pipeline.R`. It constructs synthetic CLIF tables, runs 01–05 and 07 (sourcing 06), and checks strict S/F eligibility, invalid-temperature removal, nearest future fill auditing, categorical model terms, support proportions, all five estimates, retention of late initiators as controls, and freezing of adherence weights after grace. This is software validation, not validation on the actual EHR cohort.

## Server crash workaround

Script 01 uses single-threaded data.table operations and base-R ordering/patient runs for hospital-block assignment. This avoids the grouped `:=`/`forderv` code path that crashed on the server. Overlap handling, the six-hour gap rule, block IDs and component ordering are unchanged. Runtime R/data.table/Arrow versions are printed in the log. This workaround does not establish the underlying cause of the native-library crash; a matching server environment may still require package maintenance if another operation crashes.

## Result directory layout

```text
results/run_<timestamp>_grace24h/
  aggregate/
    eda/                 # Script 05 summaries, figures and screening ranges
    analysis/            # Script 07 five-method results and weight diagnostics
  sensitive/
    cohort/              # Script 01 cohort/crosswalk
    source_profile/      # Script 02 source checks
    steroid_candidates/  # Script 03
    steroid_treatment/   # Script 04
    eda_preparation/     # Script 05 prepared individual data and fill audit
    analysis/            # Script 07 individual data, fits and weight trajectories
  logs/
```

Script 07's fourth argument is now the **run root**, not an `analysis/` subfolder. It creates `aggregate/analysis` and `sensitive/analysis` itself. The run command file supplies this correctly. Use a fresh run to avoid mixing the old all-cause endpoint/results with the new composite endpoint. No old result files are automatically moved or deleted.
