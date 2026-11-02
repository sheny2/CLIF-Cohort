# Hopkins CLIF adult pneumonia cohort

This folder starts a cohort workflow for Hopkins EHR data stored in CLIF format.
It adapts the MIMIC-IV cohort logic to the Hopkins multisite CLIF.


## Current cohort definition

The extraction script constructs hospital blocks and retains one eligible block
per patient using these rules (largely based on previous codes):

1. Combine overlapping hospitalizations and consecutive hospitalizations whose gap is no more than six hours.
2. Require valid admission/discharge timestamps, a patient identifier, and age at least 18 (but less than 120 years).
3. Require the first component encounter to have `admission_type_category == "ed"`.
4. Apply the pneumonia coding definition: principal pneumonia, or secondary pneumonia with principal sepsis or respiratory failure in the same component hospitalization.
5. Require at least 24 continuous hours of IMV, NIPPV, or HFNC at least 30 L/min, beginning within seven days of block admission.
6. Define time zero as the beginning of the first qualifying respiratory-support episode.
7. During the first 24 hours after time zero, require FiO2 documentation, exclude tracheostomy, and require an eligible P/F or S/F measurement paired to FiO2 within four hours. S/F requires SpO2 80--97%; P/F requires PaO2 40--700 and excludes ratios of 600 or greater.
8. Exclude DNI-equivalent code status (`and`, `dnr/dni`, or `dnar/dni` after normalization) recorded from seven days before through 24 hours after time zero.
9. Randomly select one eligible block per patient.



## Required CLIF files

- `clif_hospitalization.parquet`
- `clif_hospital_diagnosis.parquet`
- `clif_respiratory_support.parquet`
- `clif_vitals.parquet`
- `clif_labs.parquet`
- `clif_code_status.parquet`

For multisite attribution, `clif_adt.parquet` is strongly recommended. The
builder uses the `hospital_id` on the first chronological ADT record for each
component encounter, matching the collaborator example. If ADT is unavailable,
the builder can use a configured site column in `clif_hospitalization.parquet`
or a constant fallback for an extract already separated by site.

The Hopkins team confirmed that steroid administrations are stored in
`clif_medication_admin_intermittent.parquet`, with:

- `med_group == "steroid"`
- administration time in `admin_dttm`
- candidate systemic routes defined by
  `med_route_category %in% c("enteral", "iv")`

The minimum qualifying dose and any glucocorticoid dose-equivalence rules are
still pending clinical specification.

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
- `input_check/site_provenance_check.csv`
- `input_check/optional_medication_table_check.csv`

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
  hopkins_unspecified
```

The arguments are:

1. CLIF data directory (required unless `CLIF_DATA_DIR` is set)
2. Output directory
3. Pneumonia definition: `expanded`, `principal`, or `any`
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
Rscript 01_build_adult_pneumonia_cohort.R /path/to/hopkins/clif ./output
```

## 3. Locate glucocorticoid information

Run the discovery script before defining treatment:

```bash
Rscript 02_profile_glucocorticoid_sources.R \
  /path/to/hopkins/clif \
  ./glucocorticoid_source_check
```

Review:

- `glucocorticoid_source_check/medication_table_inventory.csv`
- `glucocorticoid_source_check/glucocorticoid_candidate_values.csv`
- `glucocorticoid_source_check/glucocorticoid_candidate_routes.csv`
- `glucocorticoid_source_check/glucocorticoid_search_dictionary.csv`

This step verifies the confirmed `med_group`, administration-time, and route
fields while also profiling medication names and categories. It does **not** yet
classify treatment because the qualifying dose threshold and dose-equivalence
rules have not been finalized.

## 4. Extract candidate steroid administrations

After building the cohort, extract route-eligible steroid administrations:

```bash
Rscript 03_extract_candidate_steroid_administrations.R \
  /path/to/hopkins/clif \
  ./output \
  ./steroid_candidates
```

This links medication administrations to the selected cohort through component
hospitalization IDs and computes `hours_from_time_zero` when the cohort contains
`respiratory_failure_start`. It writes:

- `candidate_steroid_administrations_in_final_cohort.parquet`
- `candidate_steroid_attrition.csv`
- `candidate_steroid_definition.csv`

These are candidate administrations, not the final treatment indicator. The
definition file explicitly records the dose threshold and dose-equivalence
rules as pending.

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
install.packages(c("arrow", "data.table", "dplyr"))
```

## Native-memory errors

If R exits with a libc message such as `corrupted size vs. prev_size`, record
the last progress message and run `sessionInfo()` in a fresh R session. This is
a native Arrow/data.table failure rather than a catchable R condition. The
current cohort builder avoids collecting the full ADT table and selects only
needed hospitalization/diagnosis columns to reduce that risk.
