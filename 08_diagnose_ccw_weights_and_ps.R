#!/usr/bin/env Rscript
# Standalone weight/PS/ESS diagnostics after scripts 01--05; no 07/09 outputs
# or execution are required. Existing source data and results are not changed.
# Usage:
# Rscript 08_diagnose_ccw_weights_and_ps.R CLIF_DATA_DIR [COHORT_DIR] \
#   [TREATMENT_DIR] [NEW_OUTPUT_DIR] [DATA_CUTOFF_DATE]
# Uses the original site-specific logistic models only. No spline, ridge,
# model comparison, or cross-validation is run.
# Defaults: TREATMENT_GRACE_HOURS=24; DIAGNOSTIC_MAX_PATIENTS=0 (all).
# Raw CLIF covariates/death data and the outputs of 01/04 are required;
# 02/03 support 04 and 05 supplies optional EDA.
# Output: aggregate/ (unsuppressed summaries/plots) and sensitive/ (prepared
# data, row-level weights/probabilities, extreme trajectories; keep on server).
# Script 06 is sourced as a function library; it need not be executed separately.

suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
argv <- commandArgs(TRUE)
sf <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
script_dir <- dirname(normalizePath(gsub("~+~"," ",sf,fixed=TRUE)))
source(file.path(script_dir,"06_generalized_five_methods.R"))
arg <- function(i, default) if(length(argv)>=i && nzchar(argv[i])) argv[i] else default
data_dir <- arg(1,Sys.getenv("CLIF_DATA_DIR","/projects/emulate_icu/rclif_2_1"))
cohort_dir <- arg(2,Sys.getenv("COHORT_OUTPUT_DIR",file.path(script_dir,"pipeline_output/sensitive/cohort")))
treatment_dir <- arg(3,Sys.getenv("STEROID_TREATMENT_DIR",file.path(script_dir,"pipeline_output/sensitive/steroid_treatment")))
output_dir <- arg(4,Sys.getenv("DIAGNOSTIC_OUTPUT_DIR",file.path(script_dir,"pipeline_output",
  paste0("weight_diagnostics_",format(Sys.time(),"%Y%m%d_%H%M%S")))))
cutoff_text <- arg(5,Sys.getenv("CLIF_DATA_CUTOFF_DATE","2026-10-31"))
grace_hours <- suppressWarnings(as.numeric(Sys.getenv("TREATMENT_GRACE_HOURS","24")))
max_patients <- suppressWarnings(as.integer(Sys.getenv("DIAGNOSTIC_MAX_PATIENTS","0")))
if(!is.finite(grace_hours) || grace_hours<6 || grace_hours>=672 || grace_hours %% 6 != 0)
  stop("TREATMENT_GRACE_HOURS must be a positive multiple of 6 below 672")
if(!is.finite(max_patients) || max_patients<0 || (max_patients>0 && max_patients<50))
  stop("DIAGNOSTIC_MAX_PATIENTS must be 0 (all) or at least 50")
if(!dir.exists(data_dir)) stop("CLIF directory not found: ",data_dir)
if(normalizePath(output_dir,mustWork=FALSE) %in% vapply(c(data_dir,cohort_dir,treatment_dir),
  normalizePath,character(1),mustWork=FALSE)) stop("Use a separate output directory")
aggregate_dir <- file.path(output_dir,"aggregate");sensitive_dir <- file.path(output_dir,"sensitive")
dir.create(aggregate_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(sensitive_dir,recursive=TRUE,showWarnings=FALSE)
write_agg <- function(x,name) {
  if(any(c("patient_id","hospital_block_id",".uid",".id") %in% names(x))) stop("Identifier in aggregate output: ",name)
  fwrite(x,file.path(aggregate_dir,name),na="")
}

# ---- Standalone covariate/outcome preparation (same rules as main pipeline) --
prepare_diagnostic_inputs <- function(data_dir, cohort_dir, treatment_dir, sensitive_dir,
                               shareable_dir, grace_hours, cutoff_text, max_patients) {
  step_hours <- 6L; horizon_days <- 28L
  tau <- as.integer(grace_hours / step_hours)
  M <- as.integer(horizon_days * 24 / step_hours)
  sample_seed <- 20260814L
  data_cutoff_dttm <- as.POSIXct(paste0(as.Date(cutoff_text), " 23:59:59"), tz="UTC")
  if (is.na(data_cutoff_dttm)) stop("DATA_CUTOFF_DATE must be YYYY-MM-DD")
  cohort_path <- file.path(cohort_dir, "adult_pneumonia_cohort.parquet")
  crosswalk_path <- file.path(cohort_dir, "adult_pneumonia_hospitalization_crosswalk.parquet")
  steroid_path <- file.path(treatment_dir, "qualifying_steroid_administrations.parquet")
  for (f in c(cohort_path, crosswalk_path, steroid_path))
    if (!file.exists(f)) stop("Missing upstream output: ", f, "; run scripts 01--04 first.")
  write_shareable_csv <- function(x, filename) {
    forbidden <- intersect(c("patient_id", "hospitalization_id", "hospital_block_id"), names(x))
    if (length(forbidden)) stop("Identifier in aggregate preparation output")
    fwrite(as.data.table(x), file.path(shareable_dir, filename), na="")
  }
  write_sensitive_csv <- function(x, filename) fwrite(as.data.table(x), file.path(sensitive_dir, filename), na="")
required_file <- function(filename) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) stop("Required CLIF file not found: ", path)
  path
}

read_selected <- function(filename, columns) {
  path <- required_file(filename)
  dataset <- open_dataset(path)
  schema_names <- dataset$schema$names
  missing <- setdiff(columns, schema_names)
  if (length(missing)) {
    stop(filename, " is missing: ", paste(missing, collapse = ", "))
  }
  query <- dataset |> select(all_of(columns))
  if ("hospitalization_id" %in% columns && exists("crosswalk", inherits = TRUE)) {
    ids <- unique(get("crosswalk", inherits = TRUE)$hospitalization_id)
    query <- query |> filter(.data$hospitalization_id %in% ids)
  }
  if ("patient_id" %in% columns && !"hospitalization_id" %in% columns &&
      exists("cohort", inherits = TRUE)) {
    ids <- unique(get("cohort", inherits = TRUE)$patient_id)
    query <- query |> filter(.data$patient_id %in% ids)
  }
  query |> collect() |> as.data.table()
}

clean_code <- function(x) toupper(gsub("[^A-Z0-9]", "", as.character(x)))

to_interval <- function(hours) {
  result <- rep.int(M + 1L, length(hours))
  observed <- is.finite(hours)
  result[observed] <- pmin(
    M + 1L,
    pmax(1L, as.integer(ceiling(hours[observed] / step_hours)))
  )
  result
}

# ---- Hospitalization, demographics, outcome, comorbidities ------------------

message("Building baseline hospital-block cohort and 28-day mortality outcome...")
cohort <- as.data.table(read_parquet(cohort_path, as_data_frame = TRUE))
required_cohort_columns <- c(
  "hospital_block_id", "patient_id", "admission_dttm", "discharge_dttm",
  "age_at_admission", "respiratory_failure_start",
  "respiratory_support_category", "source_site_id"
)
missing_cohort_columns <- setdiff(required_cohort_columns, names(cohort))
if (length(missing_cohort_columns)) {
  stop(
    "The cohort is not the revised block-level output; missing: ",
    paste(missing_cohort_columns, collapse = ", "),
    "\nRerun 01_build_adult_pneumonia_cohort.R."
  )
}
cohort[, `:=`(
  hospital_block_id = as.character(hospital_block_id),
  patient_id = as.character(patient_id),
  admission_dttm = as.POSIXct(admission_dttm, tz = "UTC"),
  discharge_dttm = as.POSIXct(discharge_dttm, tz = "UTC"),
  respiratory_failure_start = as.POSIXct(respiratory_failure_start, tz = "UTC"),
  index_dttm = as.POSIXct(respiratory_failure_start, tz = "UTC")
)]
if (anyDuplicated(cohort$patient_id)) {
  stop("The revised cohort must contain exactly one block per patient.")
}

crosswalk <- as.data.table(read_parquet(crosswalk_path, as_data_frame = TRUE))
required_crosswalk_columns <- c(
  "patient_id", "hospital_block_id", "hospitalization_id"
)
missing_crosswalk_columns <- setdiff(required_crosswalk_columns, names(crosswalk))
if (length(missing_crosswalk_columns)) {
  stop(
    "Hospital-block crosswalk is missing: ",
    paste(missing_crosswalk_columns, collapse = ", ")
  )
}
crosswalk[, `:=`(
  patient_id = as.character(patient_id),
  hospital_block_id = as.character(hospital_block_id),
  hospitalization_id = as.character(hospitalization_id)
)]
crosswalk <- unique(crosswalk)
if (!all(crosswalk$hospital_block_id %in% cohort$hospital_block_id)) {
  stop("Crosswalk contains blocks absent from the selected cohort.")
}
if (crosswalk[, uniqueN(hospital_block_id), by = hospitalization_id][V1 > 1L, .N]) {
  stop("A hospitalization_id maps to more than one selected hospital block.")
}

# Sampling occurs before auxiliary CLIF tables are loaded or the longitudinal
# grid is built. Because 01 already selected one block per patient, the sampled
# patient, block, and independent analysis-unit counts are identical.
full_patients <- uniqueN(cohort$patient_id)
full_blocks <- uniqueN(cohort$hospital_block_id)
full_component_hospitalizations <- uniqueN(crosswalk$hospitalization_id)
if (max_patients > 0L && full_patients > max_patients) {
  set.seed(sample_seed)
  sampled_patient_ids <- sample(
    sort(unique(cohort$patient_id)),
    size = max_patients,
    replace = FALSE
  )
  cohort <- cohort[patient_id %in% sampled_patient_ids]
  crosswalk <- crosswalk[patient_id %in% sampled_patient_ids]
  message(
    "Random sample retained ", uniqueN(cohort$patient_id),
    " patients, ", uniqueN(cohort$hospital_block_id), " hospital blocks, and ",
    uniqueN(crosswalk$hospitalization_id), " component hospitalization IDs."
  )
} else {
  message(
    "Full cohort retained because it contains only ", full_patients,
    " patients."
  )
}

sampling_summary <- data.table(
  measure = c(
    "sample_seed", "maximum_patients_requested", "full_cohort_patients",
    "full_cohort_hospital_blocks", "full_component_hospitalizations",
    "sampled_patients", "sampled_hospital_blocks",
    "sampled_component_hospitalizations"
  ),
  value = c(
    sample_seed, max_patients, full_patients, full_blocks,
    full_component_hospitalizations, uniqueN(cohort$patient_id),
    uniqueN(cohort$hospital_block_id), uniqueN(crosswalk$hospitalization_id)
  )
)
write_shareable_csv(sampling_summary, "analysis_sampling_summary.csv")

map_to_selected_blocks <- function(x) {
  x[, hospitalization_id := as.character(hospitalization_id)]
  x <- x[hospitalization_id %in% crosswalk$hospitalization_id]
  merge(
    x,
    crosswalk[, .(hospitalization_id, hospital_block_id)],
    by = "hospitalization_id",
    all = FALSE,
    sort = FALSE
  )
}

patient <- read_selected(
  "clif_patient.parquet",
  c(
    "patient_id", "sex_category", "race_category", "ethnicity_category",
    "death_dttm"
  )
)
patient[, patient_id := as.character(patient_id)]
patient <- unique(patient, by = "patient_id")

new_patient_columns <- setdiff(names(patient), names(cohort))
cohort <- merge(
  cohort,
  patient[, c("patient_id", new_patient_columns), with = FALSE],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

diagnosis <- read_selected(
  "clif_hospital_diagnosis.parquet",
  c("hospitalization_id", "diagnosis_code", "diagnosis_code_format",
    "poa_present")
)
diagnosis[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  code = clean_code(diagnosis_code),
  format = toupper(gsub("[^A-Z0-9]", "", as.character(diagnosis_code_format)))
)]
diagnosis <- diagnosis[
  tolower(trimws(as.character(poa_present))) %chin%
    c("1", "true", "t", "yes", "y")
]
diagnosis <- map_to_selected_blocks(diagnosis)

# Claims-based proxies for the chronic conditions listed in the screenshot.
diagnosis[, chronic_lung := as.integer(
  (format == "ICD10CM" & grepl("^J(4[0-7]|6[0-7])", code)) |
    (format == "ICD9CM" & grepl("^(49[0-6]|50[0-5])", code))
)]
diagnosis[, uncontrolled_diabetes := as.integer(
  (format == "ICD10CM" & grepl("^E(08|09|10|11|13)(1|65)", code)) |
    (format == "ICD9CM" & grepl("^250[0-9][23]$", code))
)]
diagnosis[, immunosuppressed := as.integer(
  (format == "ICD10CM" & grepl("^(D8[0-4]|B20|Z94|Z796)", code)) |
    (format == "ICD9CM" & grepl("^(042|279|V42|V5865)", code))
)]
diagnosis[, metastatic_solid_tumor := as.integer(
  (format == "ICD10CM" & grepl("^C(7[7-9]|80)", code)) |
    (format == "ICD9CM" & grepl("^19[6-9]", code))
)]
diagnosis[, hematologic_malignancy := as.integer(
  (format == "ICD10CM" & grepl("^(C8[1-9]|C9[0-6]|D4[5-7])", code)) |
    (format == "ICD9CM" & grepl("^(20[0-8]|2387)", code))
)]
diagnosis[, heart_failure := as.integer(
  (format == "ICD10CM" & grepl("^I50", code)) |
    (format == "ICD9CM" & grepl("^428", code))
)]
diagnosis[, chronic_kidney_disease := as.integer(
  (format == "ICD10CM" & grepl("^N18", code)) |
    (format == "ICD9CM" & grepl("^585", code))
)]

comorbidities <- diagnosis[
  ,
  .(
    chronic_lung = max(chronic_lung),
    uncontrolled_diabetes = max(uncontrolled_diabetes),
    immunosuppressed = max(immunosuppressed),
    metastatic_solid_tumor = max(metastatic_solid_tumor),
    hematologic_malignancy = max(hematologic_malignancy),
    heart_failure = max(heart_failure),
    chronic_kidney_disease = max(chronic_kidney_disease)
  ),
  by = hospital_block_id
]
cohort <- merge(
  cohort,
  comorbidities,
  by = "hospital_block_id",
  all.x = TRUE,
  sort = FALSE
)
comorbidity_names <- c(
  "chronic_lung", "uncontrolled_diabetes", "immunosuppressed",
  "metastatic_solid_tumor", "hematologic_malignancy", "heart_failure",
  "chronic_kidney_disease"
)
for (variable in comorbidity_names) setnafill(cohort, "const", 0, cols = variable)

# Use patient death_dttm for 28-day death after respiratory-failure start,
# including deaths after discharge. Fall back to an expired-discharge timestamp
# only when death_dttm is absent. Administrative censoring uses the explicit
# study data cutoff.
cohort[, event_dttm := as.POSIXct(death_dttm, tz = "UTC")]
cohort[
  is.na(event_dttm) & tolower(discharge_category) == "expired",
  event_dttm := as.POSIXct(discharge_dttm, tz = "UTC")
]
cohort[, event_hours := as.numeric(difftime(
  event_dttm,
  index_dttm,
  units = "hours"
))]
cohort[
  !is.finite(event_hours) | event_hours < 0,
  `:=`(event_dttm = as.POSIXct(NA, tz = "UTC"), event_hours = NA_real_)
]
cohort[, censor_hours := as.numeric(difftime(
  data_cutoff_dttm, index_dttm, units = "hours"
))]
if (any(cohort$censor_hours < 0, na.rm = TRUE)) {
  stop("At least one cohort time zero occurs after the data cutoff.")
}
cohort[, delta := as.integer(
  is.finite(event_hours) & event_hours <= horizon_days * 24 &
    event_hours <= censor_hours
)]
cohort[, T_event := fifelse(delta == 1L, to_interval(event_hours), M + 1L)]
cohort[, T_censor := fifelse(
  censor_hours < horizon_days * 24,
  pmin(M, pmax(1L, as.integer(floor(censor_hours / step_hours)) + 1L)),
  M + 1L
)]
cohort[, T_obs := pmin(T_event, T_censor, M)]

cohort[, `:=`(
  age_z = (as.numeric(age_at_admission) - 65) / 10,
  female = fifelse(tolower(sex_category) == "female", 1L, 0L, 0L),
  sex_missing = as.integer(
    is.na(sex_category) | !tolower(sex_category) %in% c("female", "male")
  ),
  baseline_support_imv = as.integer(respiratory_support_category == "imv"),
  baseline_support_nippv = as.integer(respiratory_support_category == "nippv"),
  admission_year_z = (as.integer(format(index_dttm, "%Y")) - 2021) / 5,
  admission_season_sin = sin(2 * pi * as.integer(format(index_dttm, "%m")) / 12),
  admission_season_cos = cos(2 * pi * as.integer(format(index_dttm, "%m")) / 12)
)]

# Initial location proxies. Hospital identity is the federating/site-stratifying
# variable and therefore is not duplicated as a patient-level model covariate.
adt <- read_selected(
  "clif_adt.parquet",
  c("hospitalization_id", "in_dttm", "out_dttm", "location_category")
)
adt <- map_to_selected_blocks(adt)
adt <- merge(
  adt,
  cohort[, .(hospital_block_id, index_dttm)],
  by = "hospital_block_id",
  all.x = TRUE,
  sort = FALSE
)
adt[, `:=`(
  in_dttm = as.POSIXct(in_dttm, tz = "UTC"),
  out_dttm = as.POSIXct(out_dttm, tz = "UTC")
)]
adt[, `:=`(
  contains_index = as.integer(
    !is.na(in_dttm) & in_dttm <= index_dttm &
      (is.na(out_dttm) | out_dttm >= index_dttm)
  ),
  prior_to_index = as.integer(!is.na(in_dttm) & in_dttm <= index_dttm),
  distance = abs(as.numeric(difftime(in_dttm, index_dttm, units = "hours")))
)]
setorder(adt, hospital_block_id, -contains_index, -prior_to_index, distance)
initial_location <- adt[, .SD[1L], by = hospital_block_id]
initial_location[, location_text := tolower(as.character(location_category))]
initial_location[, `:=`(
  initial_location_icu = as.integer(grepl("icu|intensive", location_text)),
  initial_location_ed = as.integer(grepl("emerg|ed", location_text))
)]
cohort <- merge(
  cohort,
  initial_location[, .(
    hospital_block_id,
    initial_location_icu,
    initial_location_ed
  )],
  by = "hospital_block_id",
  all.x = TRUE,
  sort = FALSE
)
setnafill(
  cohort,
  "const",
  0,
  cols = c("initial_location_icu", "initial_location_ed")
)

# ---- Qualifying glucocorticoid initiation -----------------------------------

message("Reading qualifying systemic glucocorticoid administrations...")
steroid <- as.data.table(read_parquet(steroid_path, as_data_frame = TRUE))
required_steroid <- c("hospital_block_id", "admin_dttm")
missing_steroid <- setdiff(required_steroid, names(steroid))
if (length(missing_steroid)) {
  stop("Qualifying steroid file is missing: ",
       paste(missing_steroid, collapse = ", "))
}
steroid[, `:=`(
  hospital_block_id = as.character(hospital_block_id),
  admin_dttm = as.POSIXct(admin_dttm, tz = "UTC")
)]
steroid <- unique(steroid[
  hospital_block_id %in% cohort$hospital_block_id,
  .(hospital_block_id, admin_dttm)
])
steroid <- merge(
  steroid,
  cohort[, .(hospital_block_id, index_dttm, discharge_dttm, event_dttm)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
steroid[, treatment_hours := as.numeric(difftime(
  admin_dttm, index_dttm, units = "hours"
))]
steroid <- steroid[
  is.finite(treatment_hours) & treatment_hours >= 0 &
    treatment_hours <= horizon_days * 24 & admin_dttm <= discharge_dttm &
    admin_dttm <= data_cutoff_dttm &
    (is.na(event_dttm) | admin_dttm <= event_dttm)
]
first_steroid <- steroid[order(admin_dttm), .SD[1L], by = hospital_block_id]
cohort <- merge(
  cohort,
  first_steroid[, .(hospital_block_id, treatment_dttm = admin_dttm,
                     treatment_hours)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
cohort[, S := to_interval(treatment_hours)]
cohort[S > T_event | S > T_censor, S := M + 1L]
cohort[, A_tau := as.integer(S <= tau)]

# Raw outcome counts. These are descriptive cohort quantities, not adjusted
# causal estimates. In particular, treatment-group comparisons can be affected
# by treatment-confounding and by deaths that occur before the grace period ends.
deaths_within_grace <- cohort[
  delta == 1L & is.finite(event_hours) & event_hours <= grace_hours,
  .N
]
deaths_after_grace <- cohort[
  delta == 1L & is.finite(event_hours) & event_hours > grace_hours,
  .N
]
deaths_28d <- cohort[delta == 1L, .N]
cohort_mortality_summary <- data.table(
  measure = c(
    "final sampled cohort",
    paste0("death within ", grace_hours, " hours of time zero"),
    paste0("death after ", grace_hours, " hours through day ", horizon_days),
    paste0("all deaths within ", horizon_days, " days"),
    paste0("no observed death within ", horizon_days, " days")
  ),
  n = c(
    nrow(cohort),
    deaths_within_grace,
    deaths_after_grace,
    deaths_28d,
    nrow(cohort) - deaths_28d
  )
)
cohort_mortality_summary[, percent_of_cohort := 100 * n / nrow(cohort)]
write_shareable_csv(cohort_mortality_summary, "cohort_mortality_summary.csv")

mortality_by_observed_treatment <- cohort[
  ,
  .(
    hospital_blocks = .N,
    deaths_within_grace = sum(
      delta == 1L & is.finite(event_hours) & event_hours <= grace_hours
    ),
    deaths_28d = sum(delta == 1L)
  ),
  by = .(glucocorticoid_within_grace = A_tau)
]
mortality_by_observed_treatment[, `:=`(
  deaths_within_grace_percent = 100 * deaths_within_grace / hospital_blocks,
  deaths_28d_percent = 100 * deaths_28d / hospital_blocks,
  interpretation = paste(
    "Descriptive only; observed treatment groups are not baseline randomized",
    "groups and early death can prevent treatment initiation."
  )
)]
write_shareable_csv(
  mortality_by_observed_treatment,
  "cohort_mortality_by_observed_grace_period_treatment.csv"
)

# ---- Actual Hopkins hospitals as sites --------------------------------------

cohort[, site := as.character(source_site_id)]
if (anyNA(cohort$site) || any(!nzchar(cohort$site))) {
  stop("Every analysis unit must have a nonmissing source_site_id.")
}
if (uniqueN(cohort$site) < 2L) stop("At least two hospitals are required.")
patient_sites <- unique(cohort[, .(patient_id, site)])

# ---- Person-interval grid ----------------------------------------------------

grid <- cohort[
  ,
  .(m = seq_len(T_obs)),
  by = .(
    hospital_block_id,
    patient_id,
    site,
    index_dttm,
    S,
    A_tau,
    T_event,
    T_censor,
    T_obs,
    delta,
    age_z,
    female,
    sex_missing,
    admission_year_z,
    admission_season_sin,
    admission_season_cos,
    baseline_support_imv,
    baseline_support_nippv,
    initial_location_icu,
    initial_location_ed,
    chronic_lung,
    uncontrolled_diabetes,
    immunosuppressed,
    metastatic_solid_tumor,
    hematologic_malignancy,
    heart_failure,
    chronic_kidney_disease
  )
]
grid[, decision_dttm := index_dttm + (m - 1L) * step_hours * 3600]

rolling_last <- function(measurement, time_col, value_col, output_col,
                         lookback_hours = 24) {
  if (!nrow(measurement)) {
    answer <- unique(grid[, .(hospital_block_id, m)])
    answer[, (output_col) := NA_real_]
    return(answer)
  }
  x <- measurement[
    !is.na(get(time_col)) & is.finite(as.numeric(get(value_col))),
    .(
      hospital_block_id = as.character(hospital_block_id),
      recorded_dttm = as.POSIXct(get(time_col), tz = "UTC"),
      measurement_seconds = as.numeric(get(time_col)),
      value = as.numeric(get(value_col))
    )
  ]
  setkey(x, hospital_block_id, recorded_dttm)
  answer <- x[
    grid[, .(hospital_block_id, m, decision_dttm)],
    on = .(hospital_block_id, recorded_dttm <= decision_dttm),
    mult = "last",
    .(
      hospital_block_id = i.hospital_block_id,
      m = i.m,
      decision_dttm = i.decision_dttm,
      measurement_seconds = x.measurement_seconds,
      value = x.value
    )
  ]
  answer[, age_hours := as.numeric(decision_dttm) / 3600 - measurement_seconds / 3600]
  answer[!is.finite(age_hours) | age_hours < 0 | age_hours > lookback_hours,
         value := NA_real_]
  answer[, c("decision_dttm", "measurement_seconds", "age_hours") := NULL]
  setnames(answer, "value", output_col)
  answer
}

# ---- Vitals ------------------------------------------------------------------

message("Aligning temperature, SpO2, and MAP...")
vitals <- read_selected(
  "clif_vitals.parquet",
  c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value")
)
vitals <- map_to_selected_blocks(vitals)
vitals[, vital_category := gsub(
  "[^a-z0-9]+", "_", tolower(trimws(as.character(vital_category)))
)]
vitals <- vitals[
  vital_category %chin% c(
    "temperature", "temperature_c", "temp_c", "spo2", "map",
    "mean_arterial_pressure"
  )
]
temperature <- rolling_last(
  vitals[vital_category %chin% c("temperature", "temperature_c", "temp_c")],
  "recorded_dttm", "vital_value", "temperature", 12
)
spo2 <- rolling_last(
  vitals[vital_category == "spo2"],
  "recorded_dttm", "vital_value", "spo2", 12
)
map_value <- rolling_last(
  vitals[vital_category %chin% c("map", "mean_arterial_pressure")],
  "recorded_dttm", "vital_value", "map_value", 12
)
rm(vitals); gc(verbose = FALSE)

# ---- Labs --------------------------------------------------------------------

message("Aligning laboratory covariates...")
lab_path <- required_file("clif_labs.parquet")
lab_dataset <- open_dataset(lab_path)
lab_schema <- lab_dataset$schema$names
lab_time_column <- if ("lab_collect_dttm" %in% lab_schema) {
  "lab_collect_dttm"
} else if ("lab_result_dttm" %in% lab_schema) {
  "lab_result_dttm"
} else stop("clif_labs.parquet lacks lab_collect_dttm/lab_result_dttm.")
lab_categories <- lab_dataset |>
  distinct(lab_category) |>
  collect() |>
  pull(lab_category) |>
  as.character()

normalize_lab_category <- function(x) {
  gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(x))))
}
# Explicit blood/serum aliases: do not combine venous and arterial blood gases,
# bilirubin fractions, or urine and serum measurements through broad regexes.
lab_aliases <- list(
  wbc = c("wbc", "white_blood_cell_count", "white_blood_cells"),
  glucose = c("glucose_serum", "serum_glucose", "glucose", "blood_glucose"),
  albumin = c("albumin", "albumin_serum", "serum_albumin"),
  ph = c("ph_arterial", "arterial_ph", "ph_blood_arterial"),
  paco2 = c("pco2_arterial", "arterial_pco2", "paco2",
            "carbon_dioxide_arterial", "arterial_carbon_dioxide"),
  platelet = c("platelet_count", "platelets", "platelet"),
  bilirubin = c("bilirubin_total", "total_bilirubin", "bilirubin"),
  creatinine = c("creatinine", "creatinine_serum", "serum_creatinine")
)
lab_map <- lapply(lab_aliases, function(aliases) {
  lab_categories[normalize_lab_category(lab_categories) %in% aliases]
})
lab_mapping_audit <- rbindlist(lapply(names(lab_map), function(variable) {
  matched <- lab_map[[variable]]
  data.table(
    variable = variable,
    matched_lab_category = if (length(matched)) paste(matched, collapse = ";") else NA_character_,
    mapping_available = length(matched) > 0L,
    definition = if (variable %in% c("ph", "paco2")) "arterial only" else "blood/serum"
  )
}))
write_shareable_csv(lab_mapping_audit, "laboratory_category_mapping.csv")
selected_lab_categories <- unique(unlist(lab_map, use.names = FALSE))
if (length(selected_lab_categories)) {
  labs <- lab_dataset |>
    filter(
      .data$hospitalization_id %in% crosswalk$hospitalization_id,
      .data$lab_category %in% selected_lab_categories
    ) |>
    select(
      hospitalization_id,
      all_of(lab_time_column),
      lab_category,
      lab_value_numeric
    ) |>
    collect() |>
    as.data.table()
  labs[, hospitalization_id := as.character(hospitalization_id)]
  setnames(labs, lab_time_column, "lab_dttm")
  labs <- map_to_selected_blocks(labs)
} else {
  labs <- data.table(
    hospitalization_id = character(),
    hospital_block_id = character(),
    lab_dttm = as.POSIXct(character()),
    lab_category = character(),
    lab_value_numeric = numeric()
  )
}

lab_roll <- function(name, lookback = 24) {
  rolling_last(
    labs[lab_category %in% lab_map[[name]]],
    "lab_dttm",
    "lab_value_numeric",
    name,
    lookback
  )
}
wbc <- lab_roll("wbc")
glucose <- lab_roll("glucose")
albumin <- lab_roll("albumin", 48)
ph <- lab_roll("ph")
paco2 <- lab_roll("paco2")
platelet <- lab_roll("platelet", 48)
bilirubin <- lab_roll("bilirubin", 48)
creatinine <- lab_roll("creatinine", 48)
rm(labs); gc(verbose = FALSE)

# ---- Respiratory support and GCS --------------------------------------------

message("Aligning respiratory support and GCS...")
respiratory <- read_selected(
  "clif_respiratory_support.parquet",
  c("hospitalization_id", "recorded_dttm", "device_category", "fio2_set")
)
respiratory <- map_to_selected_blocks(respiratory)
respiratory[, device_category_clean := gsub(
  "[^a-z0-9]+", "_", tolower(trimws(as.character(device_category)))
)]
respiratory[, support_score := fcase(
  device_category_clean == "room_air", 0,
  device_category_clean %chin% c("nasal_cannula", "face_mask", "trach_collar"), 1,
  device_category_clean %chin% c("high_flow_nc", "nippv", "cpap"), 2,
  device_category_clean == "imv", 3,
  default = NA_real_
)]
fio2 <- rolling_last(respiratory, "recorded_dttm", "fio2_set", "fio2", 12)
support <- rolling_last(
  respiratory,
  "recorded_dttm",
  "support_score",
  "support_score",
  12
)
rm(respiratory); gc(verbose = FALSE)

assessment <- read_selected(
  "clif_patient_assessments.parquet",
  c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value")
)
assessment <- map_to_selected_blocks(assessment)
assessment <- assessment[
  assessment_category == "gcs_total"
]
gcs <- rolling_last(
  assessment,
  "recorded_dttm",
  "numerical_value",
  "gcs",
  24
)
rm(assessment); gc(verbose = FALSE)

# ---- Merge and construct analysis covariates --------------------------------

message("Constructing longitudinal covariates and SOFA proxy...")
long <- copy(grid)
for (object_name in c(
  "temperature", "spo2", "map_value", "wbc", "glucose", "albumin", "ph",
  "paco2", "platelet", "bilirubin", "creatinine", "fio2", "support", "gcs"
)) {
  long <- merge(
    long,
    get(object_name),
    by = c("hospital_block_id", "m"),
    all.x = TRUE,
    sort = FALSE
  )
}
setorder(long, site, hospital_block_id, m)

long[fio2 > 1.5, fio2 := fio2 / 100]
long[fio2 < 0.21 | fio2 > 1, fio2 := NA_real_]
long[spo2 < 50 | spo2 > 100, spo2 := NA_real_]
long[, sf_ratio := spo2 / fio2]

# Simplified non-respiratory SOFA: coagulation, liver, cardiovascular MAP,
# renal creatinine, and neurologic GCS. Vasopressor dose and urine-output
# criteria are not available in this first-pass construction.
long[, sofa_platelet := fcase(
  platelet < 20, 4, platelet < 50, 3, platelet < 100, 2,
  platelet < 150, 1, is.finite(platelet), 0, default = NA_real_
)]
long[, sofa_bilirubin := fcase(
  bilirubin >= 12, 4, bilirubin >= 6, 3, bilirubin >= 2, 2,
  bilirubin >= 1.2, 1, is.finite(bilirubin), 0, default = NA_real_
)]
long[, sofa_map := fcase(
  map_value < 70, 1, is.finite(map_value), 0, default = NA_real_
)]
long[, sofa_creatinine := fcase(
  creatinine >= 5, 4, creatinine >= 3.5, 3, creatinine >= 2, 2,
  creatinine >= 1.2, 1, is.finite(creatinine), 0, default = NA_real_
)]
long[, sofa_gcs := fcase(
  gcs < 6, 4, gcs < 10, 3, gcs < 13, 2, gcs < 15, 1,
  is.finite(gcs), 0, default = NA_real_
)]

sofa_parts <- c(
  "sofa_platelet", "sofa_bilirubin", "sofa_map",
  "sofa_creatinine", "sofa_gcs"
)
for (variable in sofa_parts) {
  long[is.na(get(variable)), (variable) := 0]
}
long[, nonresp_sofa_proxy := rowSums(.SD), .SDcols = sofa_parts]

raw_tv <- c(
  "temperature", "wbc", "glucose", "albumin", "ph", "paco2",
  "sf_ratio", "support_score", "nonresp_sofa_proxy"
)

plausible_limits <- list(
  temperature = c(25, 45),
  wbc = c(0.1, 200),
  glucose = c(20, 1000),
  albumin = c(0.5, 7),
  ph = c(6.5, 8.0),
  paco2 = c(5, 200),
  sf_ratio = c(20, 700),
  support_score = c(0, 3),
  nonresp_sofa_proxy = c(0, 17)
)

covariate_qa <- rbindlist(lapply(raw_tv, function(variable) {
  limits <- plausible_limits[[variable]]
  x <- as.numeric(long[[variable]])
  x[x < limits[1L] | x > limits[2L]] <- NA_real_
  observed <- is.finite(x)
  observed_fraction <- mean(observed)
  if (any(observed)) {
    winsor <- quantile(x[observed], c(0.005, 0.995), na.rm = TRUE, names = FALSE)
    x[observed] <- pmin(pmax(x[observed], winsor[1L]), winsor[2L])
    fill <- median(x[observed], na.rm = TRUE)
  } else {
    winsor <- c(NA_real_, NA_real_)
    fill <- 0
  }
  missing_name <- paste0(variable, "_missing")
  long[, (missing_name) := as.integer(!observed)]
  x[!observed] <- fill
  if (variable %in% c("wbc", "glucose")) x <- log1p(x)
  scale_value <- sd(x)
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
  center_value <- mean(x)
  model_name <- paste0(variable, "_model")
  if (variable %in% c("support_score", "nonresp_sofa_proxy")) {
    long[, (model_name) := x]
  } else {
    long[, (model_name) := (x - center_value) / scale_value]
  }
  data.table(
    variable = variable,
    observed_fraction = observed_fraction,
    imputation_value = fill,
    winsor_lower = winsor[1L],
    winsor_upper = winsor[2L],
    model_variable = model_name,
    usable = observed_fraction >= 0.01 && uniqueN(x) > 1L
  )
}))
write_shareable_csv(covariate_qa, "covariate_availability_and_processing.csv")
write_shareable_csv(
  data.table(
    target = names(lab_map),
    matched_clif_categories = vapply(
      lab_map,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  ),
  "matched_lab_categories.csv"
)

baseline_covariates <- c(
  "age_z", "female", "sex_missing", "baseline_support_imv",
  "baseline_support_nippv", "initial_location_icu", "initial_location_ed",
  "admission_year_z", "admission_season_sin", "admission_season_cos",
  comorbidity_names
)
baseline_covariates <- baseline_covariates[
  vapply(long[, ..baseline_covariates], function(x) uniqueN(x) > 1L, logical(1))
]
usable_tv <- covariate_qa[usable == TRUE, variable]
time_varying_covariates <- c(
  covariate_qa[usable == TRUE, model_variable],
  paste0(usable_tv, "_missing")
)
time_varying_covariates <- time_varying_covariates[
  vapply(long[, ..time_varying_covariates], function(x) uniqueN(x) > 1L,
         logical(1))
]

if (!length(time_varying_covariates)) {
  stop("No requested time-varying covariates passed availability checks.")
}

# Record the exact, post-availability-check terms used to estimate the
# stabilized IPCW components. This is more informative than the candidate list
# because constant baseline variables and poorly observed longitudinal
# variables have already been removed.
ipcw_terms <- list(
  initiation_denominator = c(
    ".time", baseline_covariates, time_varying_covariates,
    paste0("lag_", time_varying_covariates)
  ),
  initiation_numerator = ".time",
  censoring_denominator = c(".time", baseline_covariates, ".trt"),
  censoring_numerator = ".time"
)
ipcw_model_specification <- rbindlist(lapply(
  names(ipcw_terms),
  function(component) {
    data.table(
      model_component = component,
      term = ipcw_terms[[component]],
      term_type = fcase(
        ipcw_terms[[component]] == ".time", "interval time",
        ipcw_terms[[component]] == ".trt", "current treatment",
        ipcw_terms[[component]] %chin% baseline_covariates, "baseline",
        grepl("^lag_", ipcw_terms[[component]]), "lagged time-varying",
        default = "current time-varying"
      ),
      note = if (grepl("^censoring", component)) {
        paste(
          "Administrative censoring occurs only when the 28-day horizon",
          "extends beyond the prespecified data cutoff."
        )
      } else {
        paste(
          "The analysis uses the centrally aggregated interval-specific",
          "initiation numerator broadcast to all sites; denominator models",
          "are hospital-specific for both federated and pooled",
          "site-stratified CCW."
        )
      }
    )
  }
))
write_shareable_csv(ipcw_model_specification, "ipcw_model_specification.csv")

long[, init_m := as.integer(S == m & S <= M)]
long[, event_m := as.integer(delta == 1L & T_event == m)]
long[, censor_m := as.integer(
  T_censor <= M & T_censor == m & (delta == 0L | T_event > T_censor)
)]

required_analysis <- c(baseline_covariates, time_varying_covariates)
if (anyNA(long[, ..required_analysis]) ||
    any(!is.finite(as.matrix(long[, ..required_analysis])))) {
  stop("Non-finite values remain in analysis covariates.")
}
stopifnot(
  all(long$event_m + long$censor_m <= 1L),
  all(long[, .N, by = hospital_block_id]$N <= M),
  uniqueN(cohort$patient_id) == uniqueN(patient_sites$patient_id)
)

# ---- Save prepared data and diagnostics -------------------------------------

design <- list(
  source = "Hopkins multisite EHR data in CLIF 2.1 format",
  analysis_unit = "one eligible six-hour hospital block per patient",
  time_zero = "respiratory-failure start",
  treatment = "first qualifying systemic glucocorticoid administration",
  strategy_1 = "initiate qualifying glucocorticoid within the grace period",
  strategy_0 = "do not initiate qualifying glucocorticoid through 28 days",
  outcome = "death within 28 days of respiratory-failure start",
  alive_discharge = paste(
    "not censoring; patient death_dttm provides 28-day follow-up and",
    "discharge Expired is used only when death_dttm is absent"
  ),
  step_hours = step_hours,
  grace_hours = grace_hours,
  tau_intervals = tau,
  horizon_days = horizon_days,
  horizon_intervals = M,
  data_cutoff = format(data_cutoff_dttm, "%Y-%m-%d"),
  observed_hospitals = sort(unique(cohort$site)),
  maximum_sampled_patients = max_patients,
  sample_seed = sample_seed,
  full_cohort_patients = full_patients,
  full_cohort_hospital_blocks = full_blocks,
  full_component_hospitalizations = full_component_hospitalizations,
  sampled_patients = uniqueN(cohort$patient_id),
  sampled_hospital_blocks = uniqueN(cohort$hospital_block_id),
  sampled_component_hospitalizations = uniqueN(crosswalk$hospitalization_id),
  baseline_covariates = baseline_covariates,
  time_varying_covariates = time_varying_covariates,
  omitted_requested_features = c(
    "hospital type: not available as a patient-level varying field",
    "full non-respiratory SOFA: vasopressor dose and urine-output components omitted"
  )
)

  write_parquet(as.data.frame(cohort), file.path(sensitive_dir, "pneumonia_ccw_cohort.parquet"))
  write_parquet(as.data.frame(long), file.path(sensitive_dir, "pneumonia_ccw_long.parquet"))
  write_parquet(as.data.frame(crosswalk), file.path(sensitive_dir, "sampled_hospital_block_crosswalk.parquet"))
  saveRDS(design, file.path(sensitive_dir, "analysis_design.rds"))
  list(design=design, long=as.data.frame(long))
}

message("Preparing diagnostics directly from CLIF and scripts 01/04 outputs...")
message("Grace period: ",grace_hours," hours; administrative cutoff: ",cutoff_text)
inputs <- prepare_diagnostic_inputs(data_dir,cohort_dir,treatment_dir,sensitive_dir,
                                   aggregate_dir,grace_hours,cutoff_text,max_patients)
design <- inputs$design;long <- inputs$long
rm(inputs);invisible(gc())
M <- as.integer(design$horizon_intervals);tau <- as.integer(design$tau_intervals)
step <- as.numeric(design$step_hours)

# ---- Original logistic weight models and diagnostics ------------------------
prepared <- g6_prepare_data(long, "site", "hospital_block_id", "m", "event_m",
                           "censor_m", "init_m", NULL, design$baseline_covariates,
                           design$time_varying_covariates, M, lag_initial = 0)
d <- prepared$data
rm(long); invisible(gc())
terms <- list(
  initiation_denominator = c(".time", design$baseline_covariates,
    design$time_varying_covariates, prepared$lag_names),
  initiation_numerator = ".time",
  censoring_denominator = c(".time", design$baseline_covariates, ".trt"),
  censoring_numerator = ".time"
)
learner <- g6_default_binary_learner()
cfg <- list(tau = tau, horizon = M, terms = terms, trunc = c(0, 1),
            learners = list(initiation = learner, censoring = learner))
common_num <- g6_common_init_numerator(d, M)

# Extend a LOCAL COPY's return value to expose already-computed fitted hazards
# and fit objects. The source file is untouched; assert its expected interface.
# Setting truncation to [0,1] retains the full range of the raw weights.
diagnostic_components <- g6_local_components
b <- as.list(body(diagnostic_components)); last <- b[[length(b)]]
if (!is.call(last) || !identical(last[[1L]], as.name("list")))
  stop("The script-06 component interface changed; review the diagnostic adapter.")
b[[length(b)]] <- as.call(c(as.list(last), list(
  Hden = quote(Hden), Hnum = quote(Hnum), initiation_fit = quote(afit),
  censor_fit = quote(cfit), censor_numerator_fit = quote(cnfit))))
body(diagnostic_components) <- as.call(b)

warnings_log <- data.table(site = character(), stage = character(), message = character())
capture <- function(expr, site, stage) withCallingHandlers(expr, warning = function(w) {
  warnings_log <<- rbind(warnings_log, data.table(site = site, stage = stage,
                                                message = conditionMessage(w)))
  invokeRestart("muffleWarning")
})
q <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) as.numeric(quantile(x, p, names = FALSE)) else rep(NA_real_, length(p))
}
weight_summary <- function(w) {
  finite <- is.finite(w) & w > 0
  z <- w[finite]; n <- length(z)
  qs <- q(z, c(0.01, 0.5, 0.95, 0.99, 0.999))
  ess <- if (n) sum(z)^2 / sum(z^2) else NA_real_
  list(n = length(w), n_invalid = sum(!finite), minimum = if(n) min(z) else NA_real_,
       q01 = qs[1], median = qs[2], q95 = qs[3], q99 = qs[4], q999 = qs[5],
       maximum = if(n) max(z) else NA_real_, mean = if(n) mean(z) else NA_real_,
       cv = if(n > 1) sd(z) / mean(z) else NA_real_, ess = ess,
       ess_fraction = ess / n, n_gt_5 = sum(z > 5), n_gt_10 = sum(z > 10),
       n_gt_20 = sum(z > 20), n_gt_50 = sum(z > 50),
       top_1pct_weight_share = if(n) sum(sort(z, decreasing = TRUE)[seq_len(max(1, ceiling(n*.01)))])/sum(z) else NA_real_)
}
fit_summary <- function(fit, dataset, response, site, component) {
  model <- fit$model
  cf <- if (inherits(model, "glm")) coef(model) else numeric()
  data.table(site = site, component = component, rows = nrow(dataset),
    events = sum(dataset[[response]]), patients = uniqueN(dataset$.uid),
    constant_response = isTRUE(fit$constant),
    coefficients = length(cf), rank = if(length(cf)) model$rank else 0L,
    aliased_coefficients = sum(is.na(cf)),
    converged = if(length(cf)) isTRUE(model$converged) else TRUE,
    largest_abs_coefficient = if(any(is.finite(cf))) max(abs(cf[is.finite(cf)])) else NA_real_)
}
components <- list(); weight_rows <- list(); probability_rows <- list(); fit_rows <- list()
coefficient_rows <- list()
sites <- sort(unique(d$.site))
for (site in sites) {
  message("Reconstructing original models and untruncated weights: ", site)
  ds <- d[d$.site == site, , drop = FALSE]
  x <- capture(diagnostic_components(ds, cfg, common_num), site, "original_model_refit")
  if (any(!is.finite(x$W1[x$obs])) || any(!is.finite(x$W0[x$obs])))
    stop("Non-finite reconstructed weights at ", site, "; inspect model fit.")
  ad <- ds[ds$.time <= pmin(ds$.S, M), , drop = FALSE]
  cd <- ds[ds$.event == 0L, , drop = FALSE]
  fit_rows[[site]] <- rbindlist(list(
    fit_summary(x$initiation_fit, ad, ".init", site, "initiation"),
    fit_summary(x$censor_fit, cd, ".censor", site, "natural_censoring")))
  if (inherits(x$initiation_fit$model, "glm")) {
    co <- capture(summary(x$initiation_fit$model)$coefficients, site, "coefficient_summary")
    coefficient_rows[[site]] <- data.table(site = site, term = rownames(co),
      estimate = co[,1], standard_error = co[,2], z_value = co[,3], p_value = co[,4])
  }
  acell <- cbind(match(as.character(ad$.uid), x$uid), ad$.time)
  pd <- data.table(site = site, .uid = as.character(ad$.uid),
                   interval = ad$.time, initiated = ad$.init,
                   probability = x$Hden[acell], numerator = x$Hnum[acell])
  pd[, probability_before_clamp := capture(
    as.numeric(x$initiation_fit$learner$predict(x$initiation_fit$model, ad)),
    site, "original_probability_prediction")]
  pd[, observed_action_probability := fifelse(initiated == 1L, probability, 1-probability)]
  pd[, stabilized_increment := fifelse(initiated == 1L,
                                      numerator/probability, (1-numerator)/(1-probability))]
  pd[, within_grace := interval <= tau]
  probability_rows[[site]] <- pd
  ids <- as.character(ds$.uid)
  cell <- cbind(match(ids, x$uid), ds$.time)
  for (arm in c("treatment", "no_treatment")) {
    W <- if(arm == "treatment") x$W1 else x$W0
    R <- if(arm == "treatment") x$R1 else x$R0
    # Exactly match script 06: clipping thresholds use all observed cells,
    # including cells outside the strategy-specific risk set.
    lim <- q(W[x$obs], c(.01, .99))
    w <- W[cell]; used <- R[cell]
    weight_rows[[paste(site, arm)]] <- data.table(
      site = site, .uid = ids, interval = ds$.time, strategy = arm,
      in_risk_set = used, raw_weight = w,
      truncated_weight = pmin(pmax(w, lim[1]), lim[2]),
      lower_limit = lim[1], upper_limit = lim[2])
  }
  x[c("Hden", "Hnum", "initiation_fit", "censor_fit", "censor_numerator_fit")] <- NULL
  components[[site]] <- x
}
weights_all <- rbindlist(weight_rows)
probabilities <- rbindlist(probability_rows)
write_agg(rbindlist(fit_rows), "original_model_fit_diagnostics.csv")
if (length(coefficient_rows)) write_agg(rbindlist(coefficient_rows), "original_initiation_model_coefficients.csv")
weights <- weights_all[in_risk_set == TRUE]
write_parquet(weights_all, file.path(sensitive_dir, "weights_by_patient_interval.parquet"))
write_parquet(probabilities, file.path(sensitive_dir, "initiation_probabilities.parquet"))
weights_long <- melt(weights, id.vars = c("site", ".uid", "interval", "strategy"),
                     measure.vars = c("raw_weight", "truncated_weight"),
                     variable.name = "weight_version", value.name = "weight")
by_site_time <- weights_long[, weight_summary(weight), by = .(site, strategy, weight_version, interval)]
pooled_time <- weights_long[, weight_summary(weight), by = .(strategy, weight_version, interval)]
pooled_time[, site := "All sites"]
ess <- rbindlist(list(by_site_time, pooled_time), use.names = TRUE)
ess[, interval_end_hours := interval * step]
write_agg(ess, "weight_and_patient_ess_by_interval.csv")
overall <- weights_long[, weight_summary(weight), by = .(site, strategy, weight_version)]
overall[, unit := "person-intervals; not patient ESS"]
write_agg(overall, "weight_summary_all_person_intervals.csv")
clipping <- weights[, .(risk_set_person_intervals = .N,
  percent_clipped_low = 100*mean(raw_weight < lower_limit),
  percent_clipped_high = 100*mean(raw_weight > upper_limit),
  lower_limit = lower_limit[1], upper_limit = upper_limit[1]), by = .(site, strategy)]
write_agg(clipping, "truncation_thresholds_and_fraction_clipped.csv")

prob_summary <- function(p, action) list(rows = length(p), minimum = min(p),
  q01 = q(p,.01), median = median(p), q99 = q(p,.99), maximum = max(p),
  percent_p_lt_001 = 100*mean(p < .01), percent_p_gt_099 = 100*mean(p > .99),
  percent_observed_action_probability_lt_001 = 100*mean(action < .01))
write_agg(probabilities[, prob_summary(probability, observed_action_probability),
  by = .(site, interval, initiated)], "initiation_probability_summary.csv")
write_agg(probabilities[, .(patients_at_risk = uniqueN(.uid), initiations = sum(initiated),
  initiation_percent = 100*mean(initiated),
  percent_predictions_clamped = 100*mean(probability_before_clamp < 1e-6 |
                                         probability_before_clamp > 1-1e-6)),
  by = .(site, interval)], "initiation_risk_set_counts.csv")
missing_flags <- grep("_missing$", design$time_varying_covariates, value = TRUE)
if (length(missing_flags)) {
  decision_data <- as.data.table(d[d$.time <= pmin(d$.S,tau), , drop = FALSE])
  miss <- melt(decision_data[, c(".site", ".time", missing_flags), with = FALSE],
               id.vars = c(".site", ".time"), variable.name = "variable", value.name = "missing")
  miss <- miss[, .(risk_set_rows = .N, missing_rows = sum(missing),
                   missing_percent = 100*mean(missing)), by = .(.site, .time, variable)]
  setnames(miss, c(".site", ".time"), c("site", "interval"))
  write_agg(miss, "untreated_decision_risk_set_missingness.csv")
}
write_agg(probabilities[, .(rows = .N, largest_increment = max(stabilized_increment),
  q99_increment = q(stabilized_increment,.99),
  percent_increment_gt_5 = 100*mean(stabilized_increment > 5)),
  by = .(site, within_grace, initiated)], "initiation_weight_increment_summary.csv")

# Identify extreme trajectories and descriptive covariate associations. These
# associations are diagnostic, not causal explanations of treatment effects.
top <- weights[, .(peak_raw_weight = max(raw_weight)), by = .(site, strategy, .uid)]
setorder(top, site, strategy, -peak_raw_weight)
top <- top[, head(.SD, 20L), by = .(site, strategy)]
trajectory <- merge(weights, top, by = c("site", "strategy", ".uid"))
trace_vars <- unique(c(design$baseline_covariates, design$time_varying_covariates, prepared$lag_names))
covdata <- as.data.table(d)[, c(".uid", ".time", trace_vars), with = FALSE]
covdata[, .uid := as.character(.uid)]; setnames(covdata, ".time", "interval")
trajectory <- merge(trajectory, covdata, by = c(".uid", "interval"), all.x = TRUE)
trajectory <- merge(trajectory, probabilities, by = c("site", ".uid", "interval"), all.x = TRUE)
write_parquet(trajectory, file.path(sensitive_dir, "top_weight_patient_trajectories.parquet"))
profiles <- rbindlist(lapply(sites, function(site) {
  selected_site <- site
  rbindlist(lapply(c("treatment", "no_treatment"), function(arm) {
    z <- merge(weights[site == selected_site & strategy == arm], covdata,
               by = c(".uid", "interval"))
    threshold <- q(z$raw_weight, .99)
    high <- z$raw_weight >= threshold
    rbindlist(lapply(trace_vars, function(v) {
      value <- z[[v]]; s <- sd(value)
      data.table(site = site, strategy = arm, variable = v, tail_threshold = threshold,
        tail_rows = sum(high), other_rows = sum(!high),
        mean_in_tail = mean(value[high]), mean_other = mean(value[!high]),
        standardized_difference = if(is.finite(s) && s > 0 && any(!high))
          (mean(value[high])-mean(value[!high]))/s else NA_real_)
    }))
  }))
}))
write_agg(profiles, "extreme_weight_covariate_profiles.csv")

# Reuse the same fitted probabilities to vary only truncation; this is a
# sensitivity analysis, not a search for a preferred treatment-effect sign.
settings <- list(none = c(0,1), p005_p995 = c(.005,.995),
                 original_p01_p99 = c(.01,.99), p025_p975 = c(.025,.975))
sensitivity <- rbindlist(lapply(names(settings), function(label) {
  stats <- lapply(components, function(x) {
    z <- x
    for (nm in c("W1", "W0")) {
      lim <- q(z[[nm]][z$obs], settings[[label]])
      z[[nm]] <- pmin(pmax(z[[nm]],lim[1]),lim[2])
    }
    g6_local_stats(z)
  })
  fit <- g6_central_ccw(stats)
  data.table(truncation = label, death_risk_treatment = fit$psi1,
    death_risk_no_treatment = fit$psi0, risk_difference = fit$RD,
    risk_ratio = fit$RR, rd_ci_lower_fixed_weights = fit$CI_RD[1],
    rd_ci_upper_fixed_weights = fit$CI_RD[2])
}))
write_agg(sensitivity, "ccw_truncation_sensitivity.csv")
theme_set(theme_minimal(base_size=12))
probabilities[, initiation_group := fifelse(initiated == 1L,"Initiated in interval","Did not initiate in interval")]
p <- ggplot(probabilities[within_grace == TRUE],aes(probability,fill=initiation_group))+
  geom_histogram(aes(y=after_stat(density)),binwidth=.025,position="identity",alpha=.45)+
  facet_wrap(~site,scales="free_y")+labs(title="Fitted initiation probabilities during the grace period",
    subtitle="Original model, in-sample; untreated risk-set rows only; patients may contribute multiple intervals",
    x="Probability of starting steroids in this interval",y="Density",fill=NULL)
ggsave(file.path(aggregate_dir,"plot_initiation_probability_overlap.png"),p,width=11,height=7,dpi=180)
weights_long[, log10_weight := log10(weight)]
p <- ggplot(weights_long,aes(log10_weight,color=weight_version))+
  geom_freqpoly(bins=70)+facet_grid(site~strategy,scales="free_y")+
  labs(title="Raw and truncated weights in strategy-specific risk sets",
       x="log10(weight): 0 = 1; 1 = 10; 2 = 100",y="Person-interval count",color=NULL)
ggsave(file.path(aggregate_dir,"plot_raw_and_truncated_weights.png"),p,width=11,height=11,dpi=180)
p <- ggplot(ess,aes(interval_end_hours,ess_fraction,color=weight_version,linetype=strategy))+
  geom_line()+facet_wrap(~site)+geom_vline(xintercept=tau*step,linetype=3)+
  labs(title="Patient risk-set ESS as a fraction of patients at each interval",
       subtitle="ESS = sum(w)^2 / sum(w^2); strategies contain overlapping clones, not independent samples",
       x="Hours since time zero",y="ESS / risk-set patients",color=NULL,linetype=NULL)
ggsave(file.path(aggregate_dir,"plot_ess_over_time.png"),p,width=12,height=8,dpi=180)
p <- ggplot(ess[interval <= tau],aes(interval_end_hours,ess,color=weight_version,linetype=strategy))+
  geom_line()+geom_point()+facet_wrap(~site,scales="free_y")+
  labs(title="Patient risk-set ESS during the grace period",x="Interval end, hours",y="ESS",color=NULL,linetype=NULL)
ggsave(file.path(aggregate_dir,"plot_ess_during_grace.png"),p,width=12,height=8,dpi=180)
write_agg(unique(warnings_log),"model_warnings.csv")
write_agg(data.table(parameter=c("data_dir","cohort_dir","treatment_dir","grace_hours",
  "horizon_days","data_cutoff","max_patients","probability_epsilon","original_truncation",
  "requires_script_07","requires_script_09","model_comparison"),value=as.character(c(
    normalizePath(data_dir),normalizePath(cohort_dir),normalizePath(treatment_dir),tau*step,
    design$horizon_days,cutoff_text,max_patients,1e-6,"0.01;0.99",FALSE,FALSE,FALSE))),"diagnostic_settings.csv")
writeLines(c(
  "Standalone run from raw CLIF tables and scripts 01/04 outputs. No 07 or 09 execution or outputs are needed.",
  "Start with model_warnings, original_model_fit_diagnostics, raw weight tails, probability overlap and patient ESS by interval.",
  "Use original_initiation_model_coefficients and extreme_weight_covariate_profiles to investigate predictors associated with instability; coefficient magnitude alone is not variable importance.",
  "PS here is P(initiation in this interval | not yet initiated, recorded history), not a baseline ever-treated score.",
  "Original hazards are fitted on ALL eligible intervals through day 28, matching script 07, not only the grace period.",
  "Weight summaries and plots use strategy-specific risk sets; clipping thresholds match script 06's all-observed-cell convention.",
  "ESS by interval is patient risk-set ESS; all-interval ESS is person-interval ESS. Never sum arm ESS or interpret repeated intervals as independent patients.",
  "Extreme-weight covariate profiles pool repeated intervals and describe associations, not causal mechanisms.",
  "Small observed-action probabilities, large weight increments, separated/singular fits and sparse site/interval cells can identify sources of instability.",
  "Only the original logistic weight models are fitted; no alternative-model comparison or cross-validation is performed.",
  "No additional confounders are invented here. New variables require clinically justified, pre-treatment definitions and extraction before fitting.",
  "Assess covariate availability in the actual untreated decision risk sets; global 28-day missingness can be driven by post-discharge intervals.",
  "Contemporaneous covariates after initiation can be treatment affected; poor cross-arm balance on these is not by itself proof of inadequate confounding control.",
  "Truncation sensitivity varies only truncation of the SAME reconstructed models. CIs treat fitted weights/cutoffs as fixed, as in 07.",
  "All aggregate plots/CSVs are unsuppressed; patient-level files are confined to sensitive/. No source files or original results were changed."
),file.path(output_dir,"START_HERE.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir,"SESSION_INFO.txt"))
writeLines(paste("Completed", format(Sys.time(), tz="UTC", usetz=TRUE)),
           file.path(output_dir,"COMPLETE.txt"))
message("Diagnostics complete: ", normalizePath(output_dir))
