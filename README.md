# Hopkins CLIF adult pneumonia cohort

This folder starts a site-native cohort workflow for Hopkins EHR data stored in
CLIF format. It adapts the validated MIMIC-in-CLIF cohort logic, but it does not
use a MIMIC path, construct pseudo-sites, or assume shifted MIMIC dates.

## Current cohort definition

The extraction script constructs hospital blocks and retains one eligible block
per patient using these rules:

1. Combine overlapping hospitalizations and consecutive hospitalizations whose
   gap is no more than six hours.
2. Require valid admission/discharge timestamps, a patient identifier, and age
   at least 18 but less than 120 years.
3. Require the first component encounter to have
   `admission_type_category == "ed"`.
4. Apply the expanded pneumonia coding definition: principal pneumonia, or
   secondary pneumonia with principal sepsis or respiratory failure in the
   same component hospitalization.
5. Require at least 24 continuous hours of IMV, NIPPV, or HFNC at least
   30 L/min, beginning within seven days of block admission.
6. Define time zero as the beginning of the first qualifying respiratory-
   support episode.
7. During the first 24 hours after time zero, require FiO2 documentation,
   exclude tracheostomy, and require an eligible P/F or S/F measurement paired
   to FiO2 within four hours. S/F requires SpO2 80--97%; P/F requires PaO2
   40--700 and excludes ratios of 600 or greater.
8. Exclude DNI-equivalent code status (`and`, `dnr/dni`, or `dnar/dni` after
   normalization) recorded from seven days before through 24 hours after time
   zero.
9. Randomly select one eligible block per patient using a fixed seed.

The initial script intentionally matches the previous primary cohort. It does
not yet impose a calendar restriction, a septic-shock exclusion, an S/F <=315
severity threshold, or a treatment definition. Those decisions should be made
and documented before the treatment/outcome analysis is added.

## Required CLIF files

- `clif_hospitalization.parquet`
- `clif_hospital_diagnosis.parquet`
- `clif_respiratory_support.parquet`
- `clif_vitals.parquet`
- `clif_labs.parquet`
- `clif_code_status.parquet`

## 1. Validate the Hopkins CLIF extract

```bash
Rscript 00_check_clif_cohort_inputs.R \
  /path/to/hopkins/clif \
  ./input_check
```

Review:

- `input_check/required_schema_check.csv`
- `input_check/sampled_category_values.csv`
- `input_check/category_compatibility_check.csv`

The category report is important because CLIF column names can conform while
site-specific category values still differ.

## 2. Build the cohort

```bash
Rscript 01_build_adult_pneumonia_cohort.R \
  /path/to/hopkins/clif \
  ./output \
  expanded \
  6 \
  32284 \
  hopkins
```

The arguments are:

1. CLIF data directory (required unless `CLIF_DATA_DIR` is set)
2. Output directory
3. Pneumonia definition: `expanded`, `principal`, or `any`
4. Hospital-block gap in hours
5. One-block-per-patient random-selection seed
6. Real source-site identifier (or set `CLIF_SITE_ID`)

The same run can be configured with environment variables:

```bash
export CLIF_DATA_DIR=/path/to/hopkins/clif
export CLIF_SITE_ID=hopkins
Rscript 01_build_adult_pneumonia_cohort.R
```

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

Every cohort and crosswalk row carries `source_site_id`. No pseudo-site is
created at this stage.

## Packages

```r
install.packages(c("arrow", "data.table", "dplyr"))
```
