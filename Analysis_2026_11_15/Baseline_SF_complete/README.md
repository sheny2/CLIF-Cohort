# Baseline S/F-complete sensitivity analysis

This is a separate analysis within `Analysis_2026_11_15/Baseline_SF_complete`. All R code and helpers needed to run it are included here. No parent scripts are sourced, no original scripts are changed, and no existing result files are overwritten.

## Question and restriction

Among the original selected hospital-block cohort, retain patients with a **complete, eligible S/F ratio at the original time zero**. Use the latest valid prior SpO2 and FiO2 readings, independently selected within the preceding 12 hours and within the selected block:

- Valid SpO2 is 50–100%; the selected reading must then be **<=97%** for inclusion.
- FiO2 must be 0.21–1. Percent-form values 21–100 are normalized before screening.
- The resulting S/F must pass the existing 20–700 screening range.
- Both measurements must be at or before the original time zero. The code does not select a later time zero or use a future reading to establish eligibility.
- A newer valid SpO2 of 98–100% excludes the patient even if an older lower saturation exists in the lookback. High readings are not removed before temporal selection to manufacture eligibility.
- Later SpO2 >97%, missing later readings, treatment initiation, improvement or discharge do not remove additional patients or six-hour rows.

This restriction defines a new baseline population. The latest uploaded main run suggests approximately **2,162/4,483** would qualify, but the script recalculates the exact number from your raw data and selected cohort; it is not hardcoded.

## Model changes

S/F is a **baseline-only covariate**, named `baseline_sf_model`, constant within each patient. It is winsorized at the 0.5th/99.5th percentiles and standardized using one baseline value per retained patient; no baseline S/F imputation or missingness indicator is used. Numeric imputation and transformations for other covariates are re-estimated in the restricted cohort. No `sf_ratio_model`, `sf_ratio_missing`, or their lags enter the models.

Screened time-varying S/F may remain in the sensitive long table for audit, but is not an analysis covariate. The EDA shows a patient-level baseline S/F distribution instead of a time-varying S/F model/missingness summary.

The outcome remains **in-hospital death or hospice discharge by day 28** in the selected block, with the earliest qualifying event used. Recorded death within the block supplies the in-hospital death time; Expired-discharge time is the fallback. Post-discharge death after another alive discharge is excluded. Other alive discharge is not censoring. Administrative cutoff handling is unchanged.

The default strategies remain initiation by **24 hours** versus later/never initiation, using six-hour intervals. Temperature nearest filling, other lookback/missingness rules, categorical support, exclusion of arterial gases, cloning/censoring boundaries, weight clipping, and the five original estimators are unchanged. This is not a complete-case analysis of every covariate: missingness indicators remain for other incompletely measured variables.

08 independently fits CCW within each site and pools contrasts by REML/Hartung–Knapp meta-analysis, also reporting fixed-effect results. Local baseline S/F transformations are recomputed from one value per local patient. Site-constant terms can be omitted. The standard Hartung–Knapp implementation is retained from the parent analysis; this folder does not change meta-analysis inference. Its potentially narrow intervals with very similar site estimates still need sensitivity assessment.

## Run from existing cohort/treatment outputs

Required R packages: `arrow`, `data.table`, `dplyr`, `ggplot2`, `survival`, and `metafor`.

On the server:

```bash
cd /projects/emulate_icu/CLIF-Cohort/Analysis_2026_11_15/Baseline_SF_complete
bash RUN_PIPELINE_COMMANDS.txt /projects/emulate_icu/CLIF-Cohort/Analysis_2026_11_15/results/YOUR_MAIN_RUN
```

Replace `YOUR_MAIN_RUN` with the actual existing main run directory. Its `sensitive/cohort` and `sensitive/steroid_treatment` files are read-only inputs. Raw data default to `../../../rclif_2_1` from this folder, equivalent to the main project's `../rclif_2_1`. Set `CLIF_DATA_DIR` if needed. Grace defaults to 24 hours and the administrative cutoff to 2026-10-31.

The command runs **05, 07 and 08**. It reuses the existing cohort/treatment definitions rather than rerunning 01–04 and selecting a potentially different source population. Copies of 01–04 and the 06 method library are included for a self-contained code archive; 06 is sourced automatically. The main cohort remains unchanged.

For another grace period, use matching script-04 treatment assignment in the source run for EDA. Script 07/08 recompute initiation timing from all qualifying doses. The script defaults are intended for the current 24-hour main run.

## Outputs

```text
Baseline_SF_complete/results/run_<timestamp>_baselineSF_grace24h/
  aggregate/
    RUN_INPUTS.txt
    eda/
    analysis/
    local_ccw_meta/
  sensitive/
    eda_preparation/
    analysis/
    local_ccw_meta/
  logs/
```

Each output branch includes only the restricted cohort. The EDA applies selection before all cohort summaries, plots, treatment/outcome tables and timing summaries. Cohort attrition appends the baseline S/F restriction. Original cohort-construction metadata are retained as source metadata.

Important files:

- `baseline_sf_selection_summary.csv`: source cohort size, mutually exclusive exclusion reasons and retained count.
- `baseline_sf_selection_by_site.csv`: source/retained/excluded counts by site.
- `baseline_sf_sensitivity_definition.csv`: exact eligibility and model-role rules.
- `baseline_sf_processing.csv`: baseline winsorization/scaling; zero missing and no imputation.
- `baseline_sf_distribution.csv` and `plots/baseline_sf_distribution.png`: one observation per patient.
- `ipcw_model_specification.csv`: actual baseline/current/lagged terms. Only `baseline_sf_model` should appear for S/F.
- `five_method_results.csv`, `six_method_results.csv`, effect plots and site forests: same methods/outcome in the restricted population.
- Sensitive `baseline_sf_selection_rows.parquet`: patient-level selection audit, including excluded source patients; retain on server.

08 refuses to run on a main-analysis design without the `baseline_only` S/F marker. This prevents accidentally using the main run as its output target. Existing parent results are not moved, deleted or rewritten.

## Validation

Run `Rscript tests/test_pipeline.R` for the synthetic integration test, plus `Rscript tests/test_baseline_sf.R` for selection edge cases. The integration fixture includes high baseline saturation, missing baseline FiO2, first saturation after time zero, and later saturation above 97% in a retained patient. It verifies matching EDA/model populations, baseline S/F constancy/completeness, absence of time-varying S/F terms and 24h/96h behavior. The copied local-meta tests can be run on the resulting synthetic run root. Software validation does not substitute for reviewing the actual EHR sensitivity results.
