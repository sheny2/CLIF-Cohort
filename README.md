# Hopkins CLIF adult pneumonia cohort

This folder starts a cohort workflow for Hopkins EHR data stored in CLIF format.
It adapts the MIMIC-IV cohort logic to the Hopkins multisite CLIF.

The current Discovery HPC source at `/projects/emulate_icu/rclif_2_1` is a
limited CLIF 2.1 extract prepared with these upstream filters:

1. `admission_dttm` in the configured January 2021--December 2026 interval
   (the extract itself was created in early September 2026).
2. Admission from the emergency department.
3. Ever received IMV, NIPPV, or HFNO.

The cohort script intentionally reapplies its ED and respiratory-support rules
because the analytic criteria are narrower than these shipping filters. The
source delivery also contains `shipped_flow_chart.csv` for reconciliation with
the shipped table denominators.


## Current cohort definition

The extraction script constructs hospital blocks and retains one eligible block
per patient using the CAPE COD-aligned rules confirmed for this emulation:

1. Combine overlapping hospitalizations and consecutive hospitalizations whose gap is no more than six hours.
2. Require valid admission/discharge timestamps, a patient identifier, and age at least 18 (but less than 120 years).
3. Require the first component encounter to have `admission_type_category == "ed"`.
4. Require an ICD-10-CM code from `J12.*`, `J13-J16.*`, or `J18.*` in any diagnosis position with `poa_present == 1`. Influenza (`J09-J11`), pneumonia in diseases classified elsewhere (`J17.*`), aspiration (`J69.*`), and ventilator-associated pneumonia (`J95.851`) are exclusionary when present on admission. `U07.1` and organism/sepsis codes do not qualify without a qualifying pneumonia diagnosis; `J12.82` qualifies through `J12.*`.
5. Define time zero as the first eligible initiation of IMV or NIPPV, or HFNO initiation with a contemporaneous S/F ratio below 315 while FiO2 is greater than 0.50. SpO2 and HFNO FiO2 records are paired within four hours. No future respiratory-support duration is required.
6. Exclude DNR/DNI-equivalent code status (`and`, `dnr`, `dni`, `dnar`, `dnr/dni`, or `dnar/dni` after normalization) documented during the seven days before and through time zero.
7. Exclude the concurrent-shock proxy at time zero: any active norepinephrine/epinephrine on NIPPV or HFNO; norepinephrine/epinephrine at least 0.25 mcg/kg/min on IMV; any active vasopressin; or lactate at least 4 mmol/L. The most recent medication and lactate records in the preceding six hours are used. Missing or nonconvertible shock data are retained and flagged rather than treated as evidence of shock.
8. Exclude any patient whose recorded `death_dttm` occurs before time zero.
9. Randomly select one eligible block per patient.



## Required CLIF files

- `clif_patient.parquet`
- `clif_hospitalization.parquet`
- `clif_hospital_diagnosis.parquet`
- `clif_respiratory_support.parquet`
- `clif_vitals.parquet`
- `clif_labs.parquet`
- `clif_code_status.parquet`
- `clif_medication_admin_continuous.parquet`

For multisite attribution, `clif_adt.parquet` is strongly recommended. The
builder uses the `hospital_id` on the first chronological ADT record for each
component encounter, matching the collaborator example. If ADT is unavailable,
the builder can use a configured site column in `clif_hospitalization.parquet`
or a constant fallback for an extract already separated by site.

The Hopkins team confirmed that steroid administrations are stored in
`clif_medication_admin_intermittent.parquet`, with:

- `med_group %in% c("steroid", "steroids")` (both forms are accepted to
  accommodate the confirmed specification and possible stored vocabulary)
- administration time in `admin_dttm`
- candidate systemic routes defined by
  `med_route_category %in% c("enteral", "iv")`

Treatment is defined using these per-administration minimums:

| `med_category` | Minimum dose | Eligible `med_route_category` |
|---|---:|---|
| dexamethasone | 6 mg | enteral, iv |
| hydrocortisone | 50 mg | iv |
| prednisone | 40 mg | enteral |
| methylprednisolone | 32 mg | enteral, iv |
| prednisolone | 40 mg | enteral, iv |

The hydrocortisone rule is 50 mg per administration, reflecting typical q6h or
q8h dosing and the intended threshold of approximately more than 150 mg/day of
hydrocortisone-equivalent exposure. The treatment grace period remains a
separate target-trial parameter.

## 1. Validate the Hopkins CLIF extract

```bash
Rscript 00_check_clif_cohort_inputs.R \
  /path/to/hopkins/clif \
  ./pipeline_output/sensitive/input_checks
```

Review:

- `pipeline_output/sensitive/input_checks/required_schema_check.csv`
- `pipeline_output/sensitive/input_checks/sampled_category_values.csv`
- `pipeline_output/sensitive/input_checks/category_compatibility_check.csv`
- `pipeline_output/sensitive/input_checks/site_provenance_check.csv`
- `pipeline_output/sensitive/input_checks/optional_medication_table_check.csv`

The category report is important because CLIF column names can conform while
site-specific category values still differ.

## 2. Build the cohort

```bash
Rscript 01_build_adult_pneumonia_cohort.R \
  /path/to/hopkins/clif \
  ./pipeline_output/sensitive/cohort \
  poa_any \
  6 \
  32284 \
  hopkins_unspecified
```

The arguments are:

1. CLIF data directory (required unless `CLIF_DATA_DIR` is set)
2. Output directory
3. Pneumonia definition compatibility argument: `poa_any` only
4. Hospital-block gap in hours
5. One-block-per-patient random-selection seed
6. Fallback site identifier (or set `CLIF_SITE_ID`); use only when ADT/site
   provenance is missing or the input directory already represents one site
7. Optional site column in `clif_hospitalization.parquet` (or set
   `CLIF_SITE_COLUMN`); an explicit value overrides ADT-derived attribution

The same run can be configured with environment variables:

```bash
export CLIF_DATA_DIR=/path/to/hopkins/clif
export CLIF_SITE_ID=hopkins_unspecified
Rscript 01_build_adult_pneumonia_cohort.R
```

### Multisite 

By default, the script assigns each component encounter to its first ADT `hospital_id` and assigns a hospital block to the first component encounter's hospital. 
The crosswalk retains `component_site_id`, while the block-level output includes
`component_site_ids` and `crosses_site_boundary`.
This permits site-stratified summaries while keeping transfers within the same patient hospital block.
Because `clif_adt` can be very large, the script waits until clinical
eligibility is known and queries ADT only for component hospitalizations in
eligible blocks. Hospital attribution does not affect cohort eligibility.

If the data is already split into separate directories, run the script once per directory with a distinct `CLIF_SITE_ID` and output directory. If the site is stored directly on hospitalization, set it explicitly, for example:

```bash
export CLIF_SITE_COLUMN=admission_hospital_id
Rscript 01_build_adult_pneumonia_cohort.R \
  /path/to/hopkins/clif \
  ./pipeline_output/sensitive/cohort
```

## 3. Locate glucocorticoid information

Run the discovery script before defining treatment:

```bash
Rscript 02_profile_glucocorticoid_sources.R \
  /path/to/hopkins/clif \
  ./pipeline_output/sensitive/glucocorticoid_source_profile
```

Review:

- `pipeline_output/sensitive/glucocorticoid_source_profile/medication_table_inventory.csv`
- `pipeline_output/sensitive/glucocorticoid_source_profile/glucocorticoid_candidate_values.csv`
- `pipeline_output/sensitive/glucocorticoid_source_profile/glucocorticoid_candidate_routes.csv`
- `pipeline_output/sensitive/glucocorticoid_source_profile/glucocorticoid_search_dictionary.csv`

This step verifies the confirmed `med_group`, administration-time, dose, and
route fields while also profiling medication names and categories.

## 4. Extract candidate steroid administrations

After building the cohort, extract route-eligible steroid administrations:

```bash
Rscript 03_extract_candidate_steroid_administrations.R \
  /path/to/hopkins/clif \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_candidates
```

This links medication administrations to the selected cohort through component
hospitalization IDs and computes `hours_from_time_zero` when the cohort contains
`respiratory_failure_start`. It writes:

- `candidate_steroid_administrations_in_final_cohort.parquet`
- `candidate_steroid_attrition.csv`
- `candidate_steroid_definition.csv`

These are candidate administrations. Drug-specific treatment qualification is
performed in the next step.

## 5. Apply the steroid dose and route definition

Run without a grace period to audit qualifying administrations:

```bash
Rscript 04_define_steroid_treatment.R \
  ./pipeline_output/sensitive/steroid_candidates \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_treatment
```

Once the target-trial grace period is finalized, provide it in hours as the
fourth argument. For example, if the finalized protocol specifies 24 hours:

```bash
Rscript 04_define_steroid_treatment.R \
  ./pipeline_output/sensitive/steroid_candidates \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_treatment \
  24
```

The script normalizes mg, mcg/ug, and g to mg, applies the five drug-specific
thresholds, records a reason for every nonqualifying administration, and writes
a one-row-per-block treatment assignment when a grace period is supplied.

## Data governance

This is a limited data set governed by the approved Discovery HPC environment.
Do not move patient-level parquet files or other row-level data off the HPC.
Only appropriately reviewed aggregate results may be moved to the SAFE desktop;
only publication-level summary results may be placed on OneDrive. Generated
patient-level directories are git-ignored. Suppressed diagnostic CSVs default
to suppressing counts below 10; this does not replace project-specific IRB or
disclosure review.

## 6. Generate aggregate cohort summaries and EDA

Run this step only on Discovery HPC. It reads patient-level cohort and treatment
files in place but writes only suppressed aggregate CSVs and PNG figures:

```bash
Rscript 05_cohort_summary_eda.R \
  /projects/emulate_icu/rclif_2_1 \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_treatment \
  ./pipeline_output/shareable/cohort_summary_eda \
  10 \
  2026-10-31
```

The fifth argument is the minimum reportable cell size. The current study
cutoff is October 31, 2026. It is also the script default and can be overridden
by the sixth argument or `CLIF_DATA_CUTOFF_DATE`. Outputs cover attrition,
hospital/site counts,
treatment uptake and timing, baseline characteristics, missingness, discharge
outcomes, mortality, and follow-up completeness. Identifier columns are blocked
from all CSV exports.

## Output organization

All generated files now live under one root:

```text
pipeline_output/
├── sensitive/
│   ├── input_checks/
│   ├── glucocorticoid_source_profile/
│   ├── cohort/
│   ├── steroid_candidates/
│   ├── steroid_treatment/
│   └── ccw_analysis/
└── shareable/
    ├── cohort_summary_eda/
    └── ccw_results/
```

Everything under `sensitive/` must remain on Discovery. The `shareable/`
branch contains only identifier-free, small-cell-suppressed summaries and
figures, including copies of cohort attrition, design, and ICD-definition
metadata. These outputs still require the usual disclosure review before
being moved from the HPC.

## 7. Run the five-method target-trial analysis

`06_generalized_five_methods.R` contains the reusable estimator library.
`07_prepare_and_run_five_methods.R` prepares the Hopkins person-interval data
and compares:

1. Federated CCW with hospital-specific nuisance models and central aggregation
   of weighted risk/event counts.
2. Pooled site-stratified CCW, implemented as the exact centralized identity
   control for the federated estimator.
3. Federated landmark IPW among patients alive and observed at 24 hours.
4. IPCW without cloning.
5. Unweighted per-protocol analysis.

The target strategies are initiation of a qualifying systemic glucocorticoid
within 24 hours versus no qualifying initiation through day 28. Follow-up uses
six-hour intervals through 28 days. Run:

```bash
Rscript 07_prepare_and_run_five_methods.R \
  /projects/emulate_icu/rclif_2_1 \
  ./pipeline_output/sensitive/cohort \
  ./pipeline_output/sensitive/steroid_treatment \
  ./pipeline_output/sensitive/ccw_analysis \
  ./pipeline_output/shareable/ccw_results \
  0 \
  2026-10-31
```

The sixth argument is the maximum patient count; `0` uses the full cohort.
Patient/block-level person-period data, crosswalks, fitted components, and RDS
objects remain under `sensitive/ccw_analysis`; the interval-specific weighted
risk/event counts transmitted by hospitals also remain there because small
event cells can occur. Identifier-free method estimates, survival curves,
effect plots, model specifications, covariate-processing summaries, hospital
summaries, warnings, covariate-balance diagnostics, and the
federated-versus-pooled identity check are written to
`shareable/ccw_results`.

The balance diagnostics use the exact estimator-specific risk sets and
weights. They include interval-level SMDs, a compact balance summary, and
separate Love plots for CCW, landmark IPW, IPCW without cloning, and the
unweighted per-protocol analysis. Federated and pooled site-stratified CCW
share one CCW plot because their local nuisance models, weights, and risk sets
are identical by design. The Love plots show the worst absolute SMD across the
four six-hour treatment-decision intervals; the landmark plot is evaluated at
24 hours. For the unweighted per-protocol method, the before/after points are
intentionally identical.

Baseline adjustment includes age, sex and a missing-sex indicator, baseline
respiratory support, initial ICU/ED location, admission year and season, and
present-on-admission comorbidity proxies. Time-varying adjustment considers
temperature, WBC, glucose, albumin, pH, PaCO2, S/F ratio, respiratory-support
intensity, and a non-respiratory SOFA proxy, with explicit missingness
indicators. Availability checks remove variables that are absent or constant.

## Principal outputs

- `adult_pneumonia_cohort.parquet` and `.csv`: one eligible block per patient
- `adult_pneumonia_eligible_blocks_before_patient_selection.parquet`
- `adult_pneumonia_hospitalization_crosswalk.parquet`
- `hospital_block_crosswalk_all.parquet`
- `adult_pneumonia_diagnoses.parquet`
- `adult_pneumonia_respiratory_eligibility.parquet`
- `adult_pneumonia_eligible_sf_pairs.parquet`
- `adult_pneumonia_eligible_pf_pairs.parquet`
- `adult_pneumonia_attrition.csv`
- `adult_pneumonia_cohort_design.csv`
- `pneumonia_icd_definition.csv`



## Packages (if needed)

```r
install.packages(c("arrow", "data.table", "dplyr", "ggplot2"))
```

## Native-memory errors

If R exits with a libc message such as `corrupted size vs. prev_size`, record
the last progress message and run `sessionInfo()` in a fresh R session. This is
a native Arrow/data.table failure rather than a catchable R condition. The
current cohort builder avoids collecting the full ADT table and selects only
needed hospitalization/diagnosis columns to reduce that risk.
