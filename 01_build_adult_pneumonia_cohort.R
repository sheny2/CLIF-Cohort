#!/usr/bin/env Rscript

# Build an ED-admitted adult pneumonia respiratory-failure cohort from Hopkins
# EHR data represented in CLIF format.
#
# Analytic unit: one randomly selected eligible six-hour hospital block per
# patient. Component hospitalization IDs remain available through a crosswalk.
#
# Default inclusion criteria:
#   1. valid completed hospital block; age >=18 and <120 years
#   2. admission_type_category == "ed" for the first component encounter
#   3. ICD-10-CM J12.*, J13-J16.*, or J18.* present on admission, in any
#      diagnosis position
#   4. time zero at initiation of IMV or NIPPV, or initiation of HFNO with
#      S/F <315 while FiO2 >0.50; no future support-duration requirement
#   5. no DNR/DNI-equivalent code status at time zero
#   6. no concurrent shock proxy at time zero: any norepinephrine/epinephrine
#      on NIPPV/HFNO, norepinephrine/epinephrine >=0.25 mcg/kg/min on IMV,
#      any vasopressin, or lactate >=4 mmol/L
#   7. no recorded death before time zero
#
# Usage:
#   Rscript 01_build_adult_pneumonia_cohort.R \
#     CLIF_DATA_DIR [OUTPUT_DIR] [CASE_DEFINITION] [BLOCK_GAP_HOURS] [SEED] \
#     [FALLBACK_SITE_ID] [HOSPITALIZATION_SITE_COLUMN]
#
# CASE_DEFINITION is retained only for command compatibility and must be
# "poa_any" (or omitted). Older expanded/principal/any definitions are retired.

rm(list = ls())

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
})

# ---- Paths and parameters ----------------------------------------------------

script_arguments <- commandArgs(trailingOnly = FALSE)
script_file_argument <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file_argument)) {
  script_file <- sub("^--file=", "", script_file_argument[1L])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  script_dir <- dirname(normalizePath(script_file))
} else {
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  script_dir <- if (is.null(source_file)) normalizePath(getwd()) else
    dirname(normalizePath(source_file))
}
source(file.path(script_dir, "pneumonia_icd_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else
  Sys.getenv("CLIF_DATA_DIR", unset = "")
if (!nzchar(data_dir)) {
  stop(
    paste(
      "CLIF_DATA_DIR is required. Supply it as the first command-line",
      "argument or set the CLIF_DATA_DIR environment variable."
    ),
    call. = FALSE
  )
}
output_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "pipeline_output", "sensitive", "cohort")

case_definition <- if (length(args) >= 3L && nzchar(args[3L])) args[3L] else
  Sys.getenv("PNEUMONIA_CASE_DEFINITION", unset = "poa_any")
case_definition <- tolower(trimws(case_definition))
if (!case_definition %in% c("poa_any", "poa")) {
  stop(
    "CASE_DEFINITION must be 'poa_any'. The older expanded/principal/any ",
    "definitions are retired for the CAPE COD emulation.",
    call. = FALSE
  )
}
case_definition <- "poa_any"

block_gap_hours <- if (length(args) >= 4L && nzchar(args[4L])) {
  suppressWarnings(as.numeric(args[4L]))
} else {
  suppressWarnings(as.numeric(Sys.getenv("HOSPITAL_BLOCK_GAP_HOURS", unset = "6")))
}
selection_seed <- if (length(args) >= 5L && nzchar(args[5L])) {
  suppressWarnings(as.integer(args[5L]))
} else {
  suppressWarnings(as.integer(Sys.getenv("COHORT_SELECTION_SEED", unset = "32284")))
}
if (!is.finite(block_gap_hours) || block_gap_hours < 0) {
  stop("BLOCK_GAP_HOURS must be nonnegative.", call. = FALSE)
}
if (!is.finite(selection_seed)) stop("SEED must be a finite integer.", call. = FALSE)

fallback_site_id <- if (length(args) >= 6L && nzchar(args[6L])) {
  args[6L]
} else {
  Sys.getenv("CLIF_SITE_ID", unset = "hopkins_unspecified")
}
fallback_site_id <- trimws(as.character(fallback_site_id))
if (!nzchar(fallback_site_id)) {
  stop("FALLBACK_SITE_ID/CLIF_SITE_ID must be nonempty.", call. = FALSE)
}
requested_site_column <- if (length(args) >= 7L && nzchar(args[7L])) {
  args[7L]
} else Sys.getenv("CLIF_SITE_COLUMN", unset = "")
requested_site_column <- trimws(as.character(requested_site_column))

age_min <- 18
age_max_exclusive <- 120
support_start_days <- 7
support_data_days <- support_start_days + 1
baseline_pairing_hours <- 4
fio2_pairing_hours <- 4
fio2_min <- 0.21
fio2_max <- 1
spo2_min <- 80
spo2_max <- 100
hfno_sf_max_exclusive <- 315
hfno_fio2_min_exclusive <- 0.50
shock_lookback_hours <- 6
lactate_lookback_hours <- 6
weight_lookback_hours <- 24
high_dose_pressor_threshold <- 0.25
lactate_exclusion_threshold <- 4

paths <- c(
  patient = file.path(data_dir, "clif_patient.parquet"),
  hospitalization = file.path(data_dir, "clif_hospitalization.parquet"),
  diagnosis = file.path(data_dir, "clif_hospital_diagnosis.parquet"),
  respiratory_support = file.path(data_dir, "clif_respiratory_support.parquet"),
  vitals = file.path(data_dir, "clif_vitals.parquet"),
  labs = file.path(data_dir, "clif_labs.parquet"),
  code_status = file.path(data_dir, "clif_code_status.parquet"),
  continuous_medication = file.path(
    data_dir, "clif_medication_admin_continuous.parquet"
  )
)
missing_files <- paths[!file.exists(paths)]
if (length(missing_files)) {
  stop("Required CLIF file(s) not found:\n  ", paste(missing_files, collapse = "\n  "),
       call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("CLIF data directory: ", normalizePath(data_dir))
message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))
message("Pneumonia case definition: ", case_definition)
message("Admission type restriction: ed")
message("Hospital-block gap: ", block_gap_hours, " hours")
message("Respiratory eligibility: IMV/NIPPV initiation, or HFNO initiation ",
        "with S/F <", hfno_sf_max_exclusive, " and FiO2 >",
        hfno_fio2_min_exclusive, "; no minimum future duration")
message("Shock lookback: ", shock_lookback_hours,
        " hours; lactate lookback: ", lactate_lookback_hours, " hours")
message("One-block-per-patient selection seed: ", selection_seed)
message("Fallback CLIF site ID: ", fallback_site_id)
if (nzchar(requested_site_column)) {
  message("Requested hospitalization site column: ", requested_site_column)
} else {
  message("Site assignment: first ADT hospital when available; otherwise an ",
          "available hospitalization site field; otherwise fallback site ID")
}

# ---- Utilities ---------------------------------------------------------------

check_columns <- function(x, required, table_name) {
  available <- if (is.character(x) && is.null(names(x))) x else names(x)
  missing <- setdiff(required, available)
  if (length(missing)) {
    stop(table_name, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

first_nonmissing <- function(x) {
  index <- which(!is.na(x) & (!is.character(x) | nzchar(x)))
  if (!length(index)) return(x[NA_integer_][1L])
  x[index[1L]]
}

last_nonmissing <- function(x) {
  index <- which(!is.na(x) & (!is.character(x) | nzchar(x)))
  if (!length(index)) return(x[NA_integer_][1L])
  x[index[length(index)]]
}

normalize_category <- function(x) {
  value <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "_", value)
}

normalize_fio2 <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  percent <- which(!is.na(value) & value > 1 & value <= 100)
  value[percent] <- value[percent] / 100
  outlier <- which(!is.na(value) & (value < fio2_min | value > fio2_max))
  value[outlier] <- NA_real_
  value
}

normalize_truth <- function(x) {
  value <- tolower(trimws(as.character(x)))
  as.integer(value %in% c("1", "true", "t", "yes", "y"))
}

assign_blocks_running_end <- function(admission, discharge, threshold_hours) {
  n <- length(admission)
  block <- integer(n)
  gap <- rep(NA_real_, n)
  running_end <- rep(as.POSIXct(NA, tz = "UTC"), n)
  current_block <- 0L
  block_end <- as.POSIXct(NA, tz = "UTC")
  for (i in seq_len(n)) {
    if (i > 1L) {
      gap[i] <- as.numeric(difftime(admission[i], block_end, units = "hours"))
    }
    begins_new <- i == 1L || is.na(gap[i]) || gap[i] > threshold_hours
    if (begins_new) {
      current_block <- current_block + 1L
      block_end <- discharge[i]
    } else if (!is.na(discharge[i]) && discharge[i] > block_end) {
      block_end <- discharge[i]
    }
    block[i] <- current_block
    running_end[i] <- block_end
  }
  list(block = block, gap = gap, running_end = running_end)
}

pair_to_nearest_fio2 <- function(measurements, fio2_records, value_name) {
  if (!nrow(measurements) || !nrow(fio2_records)) return(data.table())
  x <- copy(fio2_records)
  i <- copy(measurements)
  x[, join_dttm := fio2_dttm]
  setkey(x, hospital_block_id, join_dttm)
  paired <- x[
    i,
    on = .(hospital_block_id, join_dttm = measurement_dttm),
    roll = "nearest",
    .(
      hospital_block_id = i.hospital_block_id,
      measurement_dttm = i.measurement_dttm,
      measurement_value = i.measurement_value,
      fio2_dttm = x.fio2_dttm,
      fio2 = x.fio2
    )
  ]
  paired[, hours_from_fio2 := abs(as.numeric(difftime(
    measurement_dttm, fio2_dttm, units = "hours"
  )))]
  setnames(paired, "measurement_value", value_name)
  paired[is.finite(hours_from_fio2) & hours_from_fio2 <= fio2_pairing_hours]
}

write_parquet_dt <- function(x, path) write_parquet(as.data.frame(x), path)

read_filtered_clif <- function(path, id_column, ids, columns) {
  dataset <- open_dataset(path)
  check_columns(dataset$schema$names, columns, basename(path))
  if (!length(ids)) {
    return(as.data.table(
      as.data.frame(setNames(replicate(length(columns), logical(), simplify = FALSE),
                             columns))
    ))
  }
  dataset |>
    filter(.data[[id_column]] %in% ids) |>
    select(all_of(columns)) |>
    collect() |>
    as.data.table()
}

# ---- Read source tables ------------------------------------------------------

message("Reading required CLIF tables...")
hospitalization_schema <- open_dataset(paths[["hospitalization"]])$schema$names
hospitalization_read_columns <- unique(c(
  "hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm",
  "age_at_admission", "admission_type_category", "discharge_category",
  requested_site_column,
  intersect(c("admission_hospital_id", "admission_hospital_name", "hospital_id",
              "hospital_name", "site_id"), hospitalization_schema)
))
hospitalization_read_columns <- intersect(
  hospitalization_read_columns[nzchar(hospitalization_read_columns)],
  hospitalization_schema
)
hospitalization <- open_dataset(paths[["hospitalization"]]) |>
  select(all_of(hospitalization_read_columns)) |>
  collect() |>
  as.data.table()
diagnosis <- open_dataset(paths[["diagnosis"]]) |>
  select(hospitalization_id, diagnosis_code, diagnosis_code_format,
         diagnosis_primary, poa_present) |>
  collect() |>
  as.data.table()

check_columns(hospitalization, c(
  "hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm",
  "age_at_admission", "admission_type_category"
), "clif_hospitalization.parquet")
check_columns(diagnosis, c(
  "hospitalization_id", "diagnosis_code", "diagnosis_code_format",
  "diagnosis_primary", "poa_present"
), "clif_hospital_diagnosis.parquet")
check_columns(open_dataset(paths[["respiratory_support"]])$schema$names, c(
  "hospitalization_id", "recorded_dttm", "device_category", "fio2_set",
  "lpm_set", "tracheostomy"
), "clif_respiratory_support.parquet")
check_columns(open_dataset(paths[["vitals"]])$schema$names, c(
  "hospitalization_id", "recorded_dttm", "vital_category", "vital_value"
), "clif_vitals.parquet")
check_columns(open_dataset(paths[["labs"]])$schema$names, c(
  "hospitalization_id", "lab_collect_dttm", "lab_category", "lab_value_numeric"
), "clif_labs.parquet")
check_columns(open_dataset(paths[["code_status"]])$schema$names, c(
  "patient_id", "start_dttm", "code_status_category"
), "clif_code_status.parquet")
check_columns(open_dataset(paths[["continuous_medication"]])$schema$names, c(
  "hospitalization_id", "admin_dttm", "med_category", "med_dose",
  "med_dose_unit", "mar_action_category"
), "clif_medication_admin_continuous.parquet")

for (x in list(hospitalization, diagnosis)) {
  x[, hospitalization_id := as.character(hospitalization_id)]
}
hospitalization[, `:=`(
  patient_id = as.character(patient_id),
  admission_dttm = as.POSIXct(admission_dttm, tz = "UTC"),
  discharge_dttm = as.POSIXct(discharge_dttm, tz = "UTC")
)]
if (!"discharge_category" %in% names(hospitalization)) {
  hospitalization[, discharge_category := NA_character_]
}

raw_hospitalization_ids <- uniqueN(hospitalization$hospitalization_id, na.rm = TRUE)
raw_patients <- uniqueN(hospitalization$patient_id, na.rm = TRUE)
setorder(hospitalization, patient_id, admission_dttm, hospitalization_id)
duplicate_hospitalization_rows <- hospitalization[
  !is.na(hospitalization_id) & hospitalization_id != "", .N,
  by = hospitalization_id
][N > 1L, sum(N - 1L)]
if (!length(duplicate_hospitalization_rows)) duplicate_hospitalization_rows <- 0L
hospitalization <- hospitalization[
  !is.na(hospitalization_id) & hospitalization_id != "", .SD[1L],
  by = hospitalization_id
]

# Assign the real hospital/site associated with each component encounter. For
# a combined multisite Hopkins extract, this follows the collaborator example:
# the first nonmissing ADT hospital ordered by in_dttm. An explicitly requested
# hospitalization column takes precedence. The constant fallback is intended
# only for site-separated deliveries or records lacking hospital provenance.
site_assignment_source <- "fallback_site_id"
site_source_column <- ""
deferred_adt_site_assignment <- FALSE
hospitalization_site_candidates <- c(
  "admission_hospital_id", "admission_hospital_name", "hospital_id",
  "hospital_name", "site_id"
)

if (nzchar(requested_site_column)) {
  if (!requested_site_column %in% names(hospitalization)) {
    stop("Requested CLIF_SITE_COLUMN is absent from clif_hospitalization: ",
         requested_site_column, call. = FALSE)
  }
  hospitalization[, source_site_id := as.character(get(requested_site_column))]
  site_assignment_source <- "hospitalization_column"
  site_source_column <- requested_site_column
} else {
  adt_path <- file.path(data_dir, "clif_adt.parquet")
  adt_columns <- if (file.exists(adt_path)) open_dataset(adt_path)$schema$names else
    character()
  if (file.exists(adt_path) &&
      all(c("hospitalization_id", "hospital_id", "in_dttm") %in% adt_columns)) {
    # clif_adt can be much larger than all other cohort tables combined. Do not
    # collect it here. Hospital attribution is not an eligibility criterion, so
    # defer the ADT query until eligible blocks are known and then read ADT rows
    # only for those component hospitalization IDs.
    hospitalization[, source_site_id := NA_character_]
    deferred_adt_site_assignment <- TRUE
    site_assignment_source <- "first_adt_hospital_deferred_to_eligible_blocks"
    site_source_column <- "clif_adt.hospital_id"
  } else {
    detected_site_column <- intersect(
      hospitalization_site_candidates, names(hospitalization)
    )
    if (length(detected_site_column)) {
      detected_site_column <- detected_site_column[1L]
      hospitalization[, source_site_id := as.character(get(detected_site_column))]
      site_assignment_source <- "auto_detected_hospitalization_column"
      site_source_column <- detected_site_column
    } else {
      hospitalization[, source_site_id := NA_character_]
    }
  }
}
if (!deferred_adt_site_assignment) {
  hospitalization[
    is.na(source_site_id) | !nzchar(trimws(source_site_id)),
    source_site_id := fallback_site_id
  ]
}
hospitalization[, source_site_id := trimws(as.character(source_site_id))]
message("Site assignment source: ", site_assignment_source,
        if (nzchar(site_source_column)) paste0(" (", site_source_column, ")") else "")
if (deferred_adt_site_assignment) {
  message("ADT site query deferred until eligible hospital blocks are known.")
} else {
  message("Observed component-encounter sites: ",
          uniqueN(hospitalization$source_site_id, na.rm = TRUE))
}

# ---- Construct hospital blocks ----------------------------------------------

valid_hospitalizations <- hospitalization[
  !is.na(patient_id) & patient_id != "" &
    !is.na(admission_dttm) & !is.na(discharge_dttm) &
    admission_dttm < discharge_dttm
]
setorder(valid_hospitalizations, patient_id, admission_dttm, discharge_dttm,
         hospitalization_id)
valid_hospitalizations[, c(
  "hospital_block_number", "gap_hours_from_running_block_end",
  "running_block_end_dttm"
) := {
  assigned <- assign_blocks_running_end(admission_dttm, discharge_dttm,
                                        block_gap_hours)
  list(assigned$block, assigned$gap, assigned$running_end)
}, by = patient_id]
valid_hospitalizations[, hospital_block_id := paste0(
  patient_id, "_block_", hospital_block_number
)]
valid_hospitalizations[, component_order := seq_len(.N),
                       by = .(patient_id, hospital_block_id)]

hospital_block_crosswalk_all <- valid_hospitalizations[, .(
  patient_id, hospital_block_id, hospital_block_number, hospitalization_id,
  component_order,
  source_site_id = as.character(source_site_id),
  component_site_id = as.character(source_site_id),
  hospitalization_admission_dttm = admission_dttm,
  hospitalization_discharge_dttm = discharge_dttm,
  gap_hours_from_running_block_end
)]

hospital_blocks <- valid_hospitalizations[, .(
  index_hospitalization_id = hospitalization_id[1L],
  source_site_id = as.character(first_nonmissing(source_site_id)),
  component_site_ids = paste(unique(source_site_id[!is.na(source_site_id) &
                                                      nzchar(source_site_id)]),
                             collapse = ";"),
  crosses_site_boundary = as.integer(uniqueN(source_site_id[!is.na(source_site_id) &
                                                              nzchar(source_site_id)]) > 1L),
  admission_dttm = min(admission_dttm),
  discharge_dttm = max(discharge_dttm),
  age_at_admission = suppressWarnings(as.numeric(first_nonmissing(age_at_admission))),
  admission_type_category = as.character(first_nonmissing(admission_type_category)),
  discharge_category = as.character(last_nonmissing(discharge_category)),
  component_hospitalizations = .N,
  component_hospitalization_ids = paste(hospitalization_id, collapse = ";"),
  contains_overlapping_hospitalizations = as.integer(any(
    is.finite(gap_hours_from_running_block_end) &
      gap_hours_from_running_block_end < 0
  )),
  contains_link_within_gap_threshold = as.integer(any(
    is.finite(gap_hours_from_running_block_end) &
      gap_hours_from_running_block_end <= block_gap_hours
  ))
), by = .(patient_id, hospital_block_id, hospital_block_number)]
hospital_blocks[, admission_type_category_clean := normalize_category(
  admission_type_category
)]
setorder(hospital_blocks, patient_id, admission_dttm, hospital_block_id)

# ---- Diagnosis eligibility ---------------------------------------------------

message("Applying POA ICD-10-CM pneumonia definition...")
icd <- derive_pneumonia_encounter_flags(diagnosis)
diagnoses_classified <- icd$diagnoses
diagnosis_encounters <- icd$encounters

diagnosis_encounters <- merge(
  diagnosis_encounters,
  hospital_block_crosswalk_all[, .(
    patient_id, hospital_block_id, hospitalization_id
  )],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
diagnoses_classified <- merge(
  diagnoses_classified,
  hospital_block_crosswalk_all[, .(
    patient_id, hospital_block_id, hospitalization_id
  )],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
diagnoses_classified <- diagnoses_classified[
  is_qualifying_pneumonia_code == 1L |
    is_influenza_pneumonia_code == 1L |
    is_explicitly_excluded_pneumonia_code == 1L
]

diagnosis_summary_by_block <- diagnosis_encounters[, .(
  qualifying_poa_pneumonia = as.integer(any(poa_pneumonia_case == 1L)),
  pneumonia_primary = as.integer(any(poa_pneumonia_primary == 1L)),
  pneumonia_case_hospitalizations = uniqueN(
    hospitalization_id[poa_pneumonia_case == 1L]
  ),
  influenza_pneumonia_exclusion = as.integer(any(
    influenza_pneumonia_present == 1L
  )),
  other_pneumonia_code_exclusion = as.integer(any(
    explicitly_excluded_pneumonia_present == 1L
  )),
  pneumonia_case_pathways = "qualifying pneumonia code present on admission"
), by = hospital_block_id]
diagnosis_summary_by_block <- diagnosis_summary_by_block[
  qualifying_poa_pneumonia == 1L & influenza_pneumonia_exclusion == 0L &
    other_pneumonia_code_exclusion == 0L
]
pneumonia_code_summary <- diagnoses_classified[qualifies_poa_pneumonia == 1L, .(
  pneumonia_code_count = uniqueN(paste(
    diagnosis_code_format_clean, diagnosis_code_clean, sep = ":"
  )),
  pneumonia_codes = paste(sort(unique(paste(
    diagnosis_code_format_clean, diagnosis_code_original, sep = ":"
  ))), collapse = ";")
), by = hospital_block_id]
diagnosis_summary_by_block <- merge(
  diagnosis_summary_by_block, pneumonia_code_summary,
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)

# ---- Demographic, ED, and diagnosis candidate blocks ------------------------

blocks_with_age <- hospital_blocks[!is.na(age_at_admission)]
adult_blocks <- blocks_with_age[
  age_at_admission >= age_min & age_at_admission < age_max_exclusive
]
ed_blocks <- adult_blocks[admission_type_category_clean == "ed"]
diagnosis_candidate_blocks <- merge(
  ed_blocks, diagnosis_summary_by_block,
  by = "hospital_block_id", all = FALSE, sort = FALSE
)
if (!nrow(diagnosis_candidate_blocks)) {
  stop("No blocks remain after age, ED-admission, and pneumonia-code criteria.")
}
candidate_crosswalk <- hospital_block_crosswalk_all[
  hospital_block_id %in% diagnosis_candidate_blocks$hospital_block_id
]

# ---- Respiratory-support initiation and time zero ----------------------------

message("Constructing respiratory-support initiation events...")
respiratory <- read_filtered_clif(
  paths[["respiratory_support"]], "hospitalization_id",
  unique(candidate_crosswalk$hospitalization_id),
  c("hospitalization_id", "recorded_dttm", "device_category", "fio2_set",
    "lpm_set", "tracheostomy")
)
respiratory[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  recorded_dttm = as.POSIXct(recorded_dttm, tz = "UTC")
)]
resp <- merge(
  respiratory, candidate_crosswalk[, .(hospitalization_id, hospital_block_id)],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
resp <- merge(
  resp,
  diagnosis_candidate_blocks[, .(
    hospital_block_id, block_admission_dttm = admission_dttm,
    block_discharge_dttm = discharge_dttm
  )],
  by = "hospital_block_id", all = FALSE, sort = FALSE
)
resp[, `:=`(
  device_category_clean = normalize_category(device_category),
  fio2 = normalize_fio2(fio2_set),
  lpm = suppressWarnings(as.numeric(lpm_set)),
  tracheostomy_binary = normalize_truth(tracheostomy),
  support_window_end = pmin(
    block_discharge_dttm,
    block_admission_dttm + support_data_days * 24 * 60 * 60
  )
)]
resp <- resp[
  !is.na(recorded_dttm) & recorded_dttm >= block_admission_dttm &
    recorded_dttm <= support_window_end
]
resp[, support_rank := fcase(
  device_category_clean == "imv", 3L,
  device_category_clean == "nippv", 2L,
  device_category_clean == "high_flow_nc", 1L,
  default = 0L
)]
resp[, qualifying_support := as.integer(support_rank > 0L)]

# At a duplicated timestamp, retain the most intensive documented device.
support_state <- resp[order(-support_rank), .SD[1L],
                      by = .(hospital_block_id, recorded_dttm)]
setorder(support_state, hospital_block_id, recorded_dttm)
support_state[, `:=`(
  previous_device = shift(device_category_clean),
  previous_support = shift(qualifying_support, fill = 0L),
  next_dttm = shift(recorded_dttm, type = "lead")
), by = hospital_block_id]
support_state[is.na(next_dttm), next_dttm := support_window_end]
support_state[, interval_hours := pmax(
  0, as.numeric(difftime(next_dttm, recorded_dttm, units = "hours"))
)]
support_state[, support_episode_number := cumsum(
  qualifying_support == 1L & previous_support == 0L
), by = hospital_block_id]
support_state[, initiation_event := as.integer(
  qualifying_support == 1L &
    (is.na(previous_device) | device_category_clean != previous_device)
)]

support_episode_summary <- support_state[qualifying_support == 1L, .(
  respiratory_support_duration_hours = sum(interval_hours),
  respiratory_failure_stop = max(next_dttm),
  episode_support_devices = paste(sort(unique(device_category_clean)),
                                  collapse = ";")
), by = .(hospital_block_id, support_episode_number)]

candidate_events <- support_state[
  initiation_event == 1L &
    recorded_dttm <= block_admission_dttm + support_start_days * 24 * 60 * 60,
  .(
    hospital_block_id, support_episode_number,
    eligibility_dttm = recorded_dttm,
    respiratory_support_category = device_category_clean,
    tracheostomy_at_time_zero = tracheostomy_binary
  )
]
setorder(candidate_events, hospital_block_id, eligibility_dttm)
candidate_events[, eligibility_event_id := .I]
if (!nrow(candidate_events)) {
  stop("No IMV, NIPPV, or HFNO initiation events were observed.")
}

# Read baseline vitals and labs for candidate blocks. Measurements after time
# zero are not used for shock exclusion; the symmetric window below is used
# only to pair an HFNO initiation with contemporaneous SpO2 and FiO2.
event_crosswalk <- candidate_crosswalk[
  hospital_block_id %in% candidate_events$hospital_block_id
]
vitals <- read_filtered_clif(
  paths[["vitals"]], "hospitalization_id",
  unique(event_crosswalk$hospitalization_id),
  c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value")
)
vitals[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  recorded_dttm = as.POSIXct(recorded_dttm, tz = "UTC")
)]
vitals <- merge(
  vitals, event_crosswalk[, .(hospitalization_id, hospital_block_id)],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
labs <- read_filtered_clif(
  paths[["labs"]], "hospitalization_id",
  unique(event_crosswalk$hospitalization_id),
  c("hospitalization_id", "lab_collect_dttm", "lab_category",
    "lab_value_numeric")
)
labs[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  lab_collect_dttm = as.POSIXct(lab_collect_dttm, tz = "UTC")
)]
labs <- merge(
  labs, event_crosswalk[, .(hospitalization_id, hospital_block_id)],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)

hfno_events <- candidate_events[respiratory_support_category == "high_flow_nc"]
hfno_fio2 <- data.table(
  eligibility_event_id = integer(),
  fio2_dttm = as.POSIXct(character()), fio2 = numeric()
)
hfno_spo2 <- data.table(
  eligibility_event_id = integer(),
  spo2_dttm = as.POSIXct(character()), spo2 = numeric()
)
if (nrow(hfno_events)) {
  hfno_fio2 <- merge(
    hfno_events[, .(eligibility_event_id, hospital_block_id, eligibility_dttm)],
    resp[device_category_clean == "high_flow_nc" & !is.na(fio2), .(
      hospital_block_id, fio2_dttm = recorded_dttm, fio2
    )],
    by = "hospital_block_id", allow.cartesian = TRUE, sort = FALSE
  )
  hfno_fio2[, hours_from_event := abs(as.numeric(difftime(
    fio2_dttm, eligibility_dttm, units = "hours"
  )))]
  hfno_fio2 <- hfno_fio2[
    hours_from_event <= baseline_pairing_hours
  ][order(eligibility_event_id, hours_from_event)][, .SD[1L],
                                                by = eligibility_event_id]

  hfno_spo2 <- merge(
    hfno_events[, .(eligibility_event_id, hospital_block_id, eligibility_dttm)],
    vitals[normalize_category(vital_category) == "spo2", .(
      hospital_block_id, spo2_dttm = recorded_dttm,
      spo2 = suppressWarnings(as.numeric(vital_value))
    )],
    by = "hospital_block_id", allow.cartesian = TRUE, sort = FALSE
  )
  hfno_spo2 <- hfno_spo2[spo2 >= spo2_min & spo2 <= spo2_max]
  hfno_spo2[, hours_from_event := abs(as.numeric(difftime(
    spo2_dttm, eligibility_dttm, units = "hours"
  )))]
  hfno_spo2 <- hfno_spo2[
    hours_from_event <= baseline_pairing_hours
  ][order(eligibility_event_id, hours_from_event)][, .SD[1L],
                                                by = eligibility_event_id]
}

hfno_baseline <- merge(
  hfno_events,
  hfno_fio2[, .(eligibility_event_id, fio2_dttm, fio2)],
  by = "eligibility_event_id", all.x = TRUE, sort = FALSE
)
hfno_baseline <- merge(
  hfno_baseline,
  hfno_spo2[, .(eligibility_event_id, spo2_dttm, spo2)],
  by = "eligibility_event_id", all.x = TRUE, sort = FALSE
)
hfno_baseline[, `:=`(
  sf_ratio = spo2 / fio2,
  sf_fio2_pair_hours = abs(as.numeric(difftime(
    spo2_dttm, fio2_dttm, units = "hours"
  )))
)]
hfno_baseline[, hfno_eligible := as.integer(
  !is.na(sf_ratio) & sf_ratio < hfno_sf_max_exclusive &
    !is.na(fio2) & fio2 > hfno_fio2_min_exclusive &
    !is.na(sf_fio2_pair_hours) & sf_fio2_pair_hours <= fio2_pairing_hours
)]

candidate_events[, respiratory_initiation_eligible := as.integer(
  respiratory_support_category %in% c("imv", "nippv")
)]
if (nrow(hfno_baseline)) {
  candidate_events[hfno_baseline, on = "eligibility_event_id", `:=`(
    baseline_spo2 = i.spo2,
    baseline_fio2 = i.fio2,
    baseline_sf_ratio = i.sf_ratio,
    respiratory_initiation_eligible = i.hfno_eligible
  )]
}
qualifying_events <- candidate_events[respiratory_initiation_eligible == 1L]
setorder(qualifying_events, hospital_block_id, eligibility_dttm)
qualifying_episode <- qualifying_events[, .SD[1L], by = hospital_block_id]
qualifying_episode <- merge(
  qualifying_episode, support_episode_summary,
  by = c("hospital_block_id", "support_episode_number"),
  all.x = TRUE, sort = FALSE
)
qualifying_episode[, `:=`(
  respiratory_failure_start = eligibility_dttm,
  sf_24_min = baseline_sf_ratio,
  pf_24_min = NA_real_
)]
if (!nrow(qualifying_episode)) {
  stop("No blocks meet the initiation-based respiratory eligibility rule.")
}

sf_pairs <- qualifying_episode[
  respiratory_support_category == "high_flow_nc",
  .(
    hospital_block_id, eligibility_event_id,
    measurement_dttm = eligibility_dttm,
    spo2 = baseline_spo2, fio2 = baseline_fio2,
    sf_ratio = baseline_sf_ratio
  )
]
pf_pairs <- data.table(
  hospital_block_id = character(), measurement_dttm = as.POSIXct(character()),
  pao2 = numeric(), fio2 = numeric(), pf_ratio = numeric()
)

# ---- Code status and concurrent shock exclusions at time zero ---------------

episode_windows <- merge(
  qualifying_episode[, .(hospital_block_id, respiratory_failure_start,
                          respiratory_support_category)],
  diagnosis_candidate_blocks[, .(hospital_block_id, patient_id)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)

# Exclude patients whose recorded death precedes time zero. This check belongs
# in cohort construction so impossible temporal records cannot enter treatment
# assignment or downstream target-trial analyses.
patient_death <- read_filtered_clif(
  paths[["patient"]], "patient_id", unique(episode_windows$patient_id),
  c("patient_id", "death_dttm")
)
patient_death[, `:=`(
  patient_id = as.character(patient_id),
  death_dttm = as.POSIXct(death_dttm, tz = "UTC")
)]
patient_death <- patient_death[, .(
  death_dttm = if (any(!is.na(death_dttm))) {
    min(death_dttm, na.rm = TRUE)
  } else as.POSIXct(NA, tz = "UTC")
), by = patient_id]
death_timing <- merge(
  episode_windows[, .(hospital_block_id, patient_id,
                       respiratory_failure_start)],
  patient_death,
  by = "patient_id", all.x = TRUE, sort = FALSE
)
death_timing[, death_before_time_zero := as.integer(
  !is.na(death_dttm) & death_dttm < respiratory_failure_start
)]

# Last code status documented during the prior seven days, through time zero.
code_status <- read_filtered_clif(
  paths[["code_status"]], "patient_id", unique(episode_windows$patient_id),
  c("patient_id", "start_dttm", "code_status_category")
)
code_status[, `:=`(
  patient_id = as.character(patient_id),
  start_dttm = as.POSIXct(start_dttm, tz = "UTC")
)]
status_candidates <- merge(
  code_status, episode_windows,
  by = "patient_id", all = FALSE, allow.cartesian = TRUE, sort = FALSE
)
status_candidates <- status_candidates[
  start_dttm >= respiratory_failure_start - 7 * 24 * 60 * 60 &
    start_dttm <= respiratory_failure_start
]
setorder(status_candidates, hospital_block_id, start_dttm)
status_summary <- status_candidates[, .(
  code_status_category_eligibility = as.character(last_nonmissing(
    code_status_category
  ))
), by = hospital_block_id]
status_summary[, limited_code_status := as.integer(
  normalize_category(code_status_category_eligibility) %in%
    c("and", "dnr", "dni", "dnar", "dnr_dni", "dnar_dni")
)]

# Most recent lactate at or before time zero within the agreed six-hour window.
lactate_candidates <- merge(
  labs[normalize_category(lab_category) == "lactate", .(
    hospital_block_id, lactate_dttm = lab_collect_dttm,
    lactate_value = suppressWarnings(as.numeric(lab_value_numeric))
  )],
  episode_windows[, .(hospital_block_id, respiratory_failure_start)],
  by = "hospital_block_id", allow.cartesian = TRUE, sort = FALSE
)
lactate_candidates <- lactate_candidates[
  lactate_dttm <= respiratory_failure_start &
    lactate_dttm >= respiratory_failure_start - lactate_lookback_hours * 60 * 60 &
    is.finite(lactate_value)
]
setorder(lactate_candidates, hospital_block_id, -lactate_dttm)
lactate_summary <- lactate_candidates[, .SD[1L], by = hospital_block_id]

# Most recent pre-baseline weight permits conversion of non-weight-normalized
# catecholamine rates. Unsupported or missing units are retained as an explicit
# indeterminate-dose diagnostic rather than treated as evidence of shock.
weight_candidates <- merge(
  vitals[normalize_category(vital_category) == "weight_kg", .(
    hospital_block_id, weight_dttm = recorded_dttm,
    weight_kg = suppressWarnings(as.numeric(vital_value))
  )],
  episode_windows[, .(hospital_block_id, respiratory_failure_start)],
  by = "hospital_block_id", allow.cartesian = TRUE, sort = FALSE
)
weight_candidates <- weight_candidates[
  weight_dttm <= respiratory_failure_start &
    weight_dttm >= respiratory_failure_start - weight_lookback_hours * 60 * 60 &
    weight_kg > 20 & weight_kg < 400
]
setorder(weight_candidates, hospital_block_id, -weight_dttm)
weight_summary <- weight_candidates[, .SD[1L], by = hospital_block_id]

continuous_schema <- open_dataset(paths[["continuous_medication"]])$schema$names
pressor_columns <- intersect(c(
  "hospitalization_id", "admin_dttm", "med_category", "med_dose",
  "med_dose_unit", "mar_action_category", "mar_action_name"
), continuous_schema)
pressor <- read_filtered_clif(
  paths[["continuous_medication"]], "hospitalization_id",
  unique(event_crosswalk$hospitalization_id), pressor_columns
)
if (!"mar_action_name" %in% names(pressor)) pressor[, mar_action_name := NA_character_]
pressor[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  admin_dttm = as.POSIXct(admin_dttm, tz = "UTC"),
  med_category_clean = normalize_category(med_category),
  med_dose_numeric = suppressWarnings(as.numeric(med_dose)),
  med_dose_unit_clean = gsub(
    "[^a-z]", "", tolower(gsub("μ", "u", as.character(med_dose_unit), fixed = TRUE))
  ),
  action_clean = normalize_category(fifelse(
    !is.na(mar_action_category), mar_action_category, mar_action_name
  ))
)]
pressor <- pressor[med_category_clean %in%
                     c("norepinephrine", "epinephrine", "vasopressin")]
pressor <- merge(
  pressor, event_crosswalk[, .(hospitalization_id, hospital_block_id)],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
pressor <- merge(
  pressor, episode_windows[, .(hospital_block_id, respiratory_failure_start)],
  by = "hospital_block_id", allow.cartesian = TRUE, sort = FALSE
)
pressor <- pressor[
  admin_dttm <= respiratory_failure_start &
    admin_dttm >= respiratory_failure_start - shock_lookback_hours * 60 * 60
]
setorder(pressor, hospital_block_id, med_category_clean, -admin_dttm)
pressor <- pressor[, .SD[1L], by = .(hospital_block_id, med_category_clean)]
pressor[, pressor_active := as.integer(
  !grepl("stop|discontinu|hold|held|pause|cancel|not_given", action_clean)
)]
pressor <- merge(
  pressor, weight_summary[, .(hospital_block_id, weight_kg)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
pressor[, dose_mcg_kg_min := fcase(
  med_dose_unit_clean %in% c("mcgkgmin", "ugkgmin", "microgramkgmin"),
    med_dose_numeric,
  med_dose_unit_clean %in% c("mgkgmin", "milligramkgmin"),
    med_dose_numeric * 1000,
  med_dose_unit_clean %in% c("mcgkghr", "ugkghr", "microgramkghr"),
    med_dose_numeric / 60,
  med_dose_unit_clean %in% c("mcgmin", "ugmin", "microgrammin") &
    is.finite(weight_kg), med_dose_numeric / weight_kg,
  med_dose_unit_clean %in% c("mghr", "milligramhr") & is.finite(weight_kg),
    med_dose_numeric * 1000 / 60 / weight_kg,
  med_dose_unit_clean %in% c("mcghr", "ughr", "microgramhr") &
    is.finite(weight_kg), med_dose_numeric / 60 / weight_kg,
  default = NA_real_
)]

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}
pressor_summary <- pressor[, .(
  active_norepinephrine_or_epinephrine = as.integer(any(
    pressor_active == 1L & med_category_clean %in%
      c("norepinephrine", "epinephrine")
  )),
  active_vasopressin = as.integer(any(
    pressor_active == 1L & med_category_clean == "vasopressin"
  )),
  max_catecholamine_dose_mcg_kg_min = safe_max(
    dose_mcg_kg_min[pressor_active == 1L & med_category_clean %in%
                      c("norepinephrine", "epinephrine")]
  ),
  catecholamine_dose_indeterminate = as.integer(any(
    pressor_active == 1L & med_category_clean %in%
      c("norepinephrine", "epinephrine") & !is.finite(dose_mcg_kg_min)
  ))
), by = hospital_block_id]

respiratory_eligibility <- Reduce(
  function(x, y) merge(x, y, by = "hospital_block_id", all.x = TRUE,
                       sort = FALSE),
  list(qualifying_episode, status_summary,
       death_timing[, .(hospital_block_id, death_before_time_zero)],
       lactate_summary[, .(hospital_block_id, lactate_dttm, lactate_value)],
       weight_summary[, .(hospital_block_id, weight_dttm, weight_kg)],
       pressor_summary)
)
for (column in c("death_before_time_zero", "limited_code_status",
                 "active_norepinephrine_or_epinephrine",
                 "active_vasopressin", "catecholamine_dose_indeterminate")) {
  if (!column %in% names(respiratory_eligibility)) {
    respiratory_eligibility[, (column) := 0L]
  }
  respiratory_eligibility[is.na(get(column)), (column) := 0L]
}
respiratory_eligibility[, `:=`(
  exclude_pressor_niv_hfno = as.integer(
    respiratory_support_category %in% c("nippv", "high_flow_nc") &
      active_norepinephrine_or_epinephrine == 1L
  ),
  exclude_high_dose_pressor_imv = as.integer(
    respiratory_support_category == "imv" &
      is.finite(max_catecholamine_dose_mcg_kg_min) &
      max_catecholamine_dose_mcg_kg_min >= high_dose_pressor_threshold
  ),
  exclude_vasopressin = as.integer(active_vasopressin == 1L),
  exclude_high_lactate = as.integer(
    is.finite(lactate_value) & lactate_value >= lactate_exclusion_threshold
  ),
  lactate_missing_at_baseline = as.integer(!is.finite(lactate_value))
)]
respiratory_eligibility[, concurrent_shock_exclusion := as.integer(
  exclude_pressor_niv_hfno == 1L |
    exclude_high_dose_pressor_imv == 1L |
    exclude_vasopressin == 1L |
    exclude_high_lactate == 1L
)]
respiratory_eligibility[, respiratory_eligible := as.integer(
  death_before_time_zero == 0L & limited_code_status == 0L &
    concurrent_shock_exclusion == 0L
)]

# ---- Final eligibility and one block per patient ----------------------------

eligible_blocks <- merge(
  diagnosis_candidate_blocks,
  respiratory_eligibility[respiratory_eligible == 1L],
  by = "hospital_block_id", all = FALSE, sort = FALSE
)

# Add real hospital attribution only after eligibility is known. This avoids
# materializing the full five-hospital ADT table, which can exhaust memory or
# abort inside Arrow before R can raise a catchable error.
if (deferred_adt_site_assignment && nrow(eligible_blocks)) {
  eligible_block_ids <- unique(eligible_blocks$hospital_block_id)
  eligible_hospitalization_ids <- unique(hospital_block_crosswalk_all[
    hospital_block_id %in% eligible_block_ids, hospitalization_id
  ])
  message("Reading ADT hospital_id for ",
          length(eligible_hospitalization_ids),
          " hospitalization(s) in eligible blocks...")
  adt_site <- open_dataset(file.path(data_dir, "clif_adt.parquet")) |>
    filter(.data$hospitalization_id %in% eligible_hospitalization_ids) |>
    select(hospitalization_id, hospital_id, in_dttm) |>
    filter(!is.na(.data$hospital_id)) |>
    collect() |>
    as.data.table()
  adt_site[, `:=`(
    hospitalization_id = as.character(hospitalization_id),
    hospital_id = trimws(as.character(hospital_id)),
    in_dttm = as.POSIXct(in_dttm, tz = "UTC")
  )]
  setorder(adt_site, hospitalization_id, in_dttm, na.last = TRUE)
  adt_site <- adt_site[
    !is.na(hospital_id) & nzchar(hospital_id),
    .(adt_source_site_id = first_nonmissing(hospital_id)),
    by = hospitalization_id
  ]

  eligible_crosswalk_index <- which(
    hospital_block_crosswalk_all$hospital_block_id %in% eligible_block_ids
  )
  site_match <- match(
    hospital_block_crosswalk_all$hospitalization_id[eligible_crosswalk_index],
    adt_site$hospitalization_id
  )
  resolved_site <- adt_site$adt_source_site_id[site_match]
  resolved_site[
    is.na(resolved_site) | !nzchar(trimws(resolved_site))
  ] <- fallback_site_id
  hospital_block_crosswalk_all[eligible_crosswalk_index, `:=`(
    source_site_id = resolved_site,
    component_site_id = resolved_site
  )]

  eligible_site_by_block <- hospital_block_crosswalk_all[
    hospital_block_id %in% eligible_block_ids
  ][order(component_order), .(
    source_site_id = as.character(first_nonmissing(component_site_id)),
    component_site_ids = paste(unique(component_site_id[
      !is.na(component_site_id) & nzchar(component_site_id)
    ]), collapse = ";"),
    crosses_site_boundary = as.integer(uniqueN(component_site_id[
      !is.na(component_site_id) & nzchar(component_site_id)
    ]) > 1L)
  ), by = hospital_block_id]

  hospital_blocks[eligible_site_by_block, on = "hospital_block_id", `:=`(
    source_site_id = i.source_site_id,
    component_site_ids = i.component_site_ids,
    crosses_site_boundary = i.crosses_site_boundary
  )]
  eligible_blocks[eligible_site_by_block, on = "hospital_block_id", `:=`(
    source_site_id = i.source_site_id,
    component_site_ids = i.component_site_ids,
    crosses_site_boundary = i.crosses_site_boundary
  )]
  message("Observed eligible-block sites: ",
          uniqueN(eligible_site_by_block$source_site_id, na.rm = TRUE))
  rm(adt_site)
  invisible(gc(verbose = FALSE))
}

eligible_blocks[, eligible_blocks_for_patient := .N, by = patient_id]
setorder(eligible_blocks, patient_id, admission_dttm, hospital_block_id)
set.seed(selection_seed)
cohort <- eligible_blocks[, .SD[sample.int(.N, 1L)], by = patient_id]
setorder(cohort, patient_id, admission_dttm, hospital_block_id)

selected_crosswalk <- hospital_block_crosswalk_all[
  hospital_block_id %in% cohort$hospital_block_id
]
selected_diagnoses <- diagnoses_classified[
  hospital_block_id %in% cohort$hospital_block_id
]
if (anyDuplicated(cohort$patient_id)) {
  stop("Internal error: final cohort has more than one block per patient.")
}

# ---- Attrition and design audit ---------------------------------------------

count_blocks <- function(block_ids) {
  block_ids <- unique(block_ids)
  data.table(
    n_source_hospitalizations = hospital_block_crosswalk_all[
      hospital_block_id %in% block_ids, uniqueN(hospitalization_id)
    ],
    n_hospital_blocks = length(block_ids),
    n_patients = hospital_blocks[
      hospital_block_id %in% block_ids, uniqueN(patient_id)
    ]
  )
}

episode_blocks <- unique(qualifying_episode$hospital_block_id)
alive_at_time_zero_blocks <- respiratory_eligibility[
  death_before_time_zero == 0L,
  hospital_block_id
]
code_status_blocks <- respiratory_eligibility[
  death_before_time_zero == 0L & limited_code_status == 0L,
  hospital_block_id
]
no_pressor_niv_blocks <- respiratory_eligibility[
  death_before_time_zero == 0L & limited_code_status == 0L &
    exclude_pressor_niv_hfno == 0L,
  hospital_block_id
]
no_high_dose_imv_blocks <- respiratory_eligibility[
  death_before_time_zero == 0L & limited_code_status == 0L &
    exclude_pressor_niv_hfno == 0L &
    exclude_high_dose_pressor_imv == 0L,
  hospital_block_id
]
no_vasopressin_blocks <- respiratory_eligibility[
  death_before_time_zero == 0L & limited_code_status == 0L &
    exclude_pressor_niv_hfno == 0L &
    exclude_high_dose_pressor_imv == 0L & exclude_vasopressin == 0L,
  hospital_block_id
]
no_high_lactate_blocks <- respiratory_eligibility[
  death_before_time_zero == 0L & limited_code_status == 0L &
    exclude_pressor_niv_hfno == 0L &
    exclude_high_dose_pressor_imv == 0L & exclude_vasopressin == 0L &
    exclude_high_lactate == 0L,
  hospital_block_id
]

attrition <- rbindlist(list(
  data.table(step = "All source hospitalizations",
             n_source_hospitalizations = raw_hospitalization_ids,
             n_hospital_blocks = NA_integer_, n_patients = raw_patients),
  data.table(step = "Valid completed hospitalizations with patient_id",
             n_source_hospitalizations = uniqueN(valid_hospitalizations$hospitalization_id),
             n_hospital_blocks = NA_integer_,
             n_patients = uniqueN(valid_hospitalizations$patient_id)),
  cbind(data.table(step = paste0("Hospital blocks after combining gaps <= ",
                                 block_gap_hours, " hours")),
        count_blocks(hospital_blocks$hospital_block_id)),
  cbind(data.table(step = "Adult blocks (age >=18 and <120)"),
        count_blocks(adult_blocks$hospital_block_id)),
  cbind(data.table(step = "Adult blocks with admission_type_category == ed"),
        count_blocks(ed_blocks$hospital_block_id)),
  cbind(data.table(step = paste0(
    "ED adult blocks with qualifying ICD-10-CM pneumonia present on admission"
  )),
        count_blocks(diagnosis_candidate_blocks$hospital_block_id)),
  cbind(data.table(step = paste0(
    "Blocks meeting initiation criterion: IMV/NIPPV or HFNO with S/F <",
    hfno_sf_max_exclusive, " and FiO2 >", hfno_fio2_min_exclusive
  )),
        count_blocks(episode_blocks)),
  cbind(data.table(step = "Blocks without recorded death before time zero"),
        count_blocks(alive_at_time_zero_blocks)),
  cbind(data.table(step = "Blocks without DNR/DNI-equivalent status at time zero"),
        count_blocks(code_status_blocks)),
  cbind(data.table(step = "Blocks without norepinephrine/epinephrine on NIPPV/HFNO"),
        count_blocks(no_pressor_niv_blocks)),
  cbind(data.table(step = paste0(
    "Blocks without norepinephrine/epinephrine >=",
    high_dose_pressor_threshold, " mcg/kg/min on IMV"
  )), count_blocks(no_high_dose_imv_blocks)),
  cbind(data.table(step = "Blocks without active vasopressin at time zero"),
        count_blocks(no_vasopressin_blocks)),
  cbind(data.table(step = paste0("Blocks without lactate >=",
                                 lactate_exclusion_threshold, " mmol/L")),
        count_blocks(no_high_lactate_blocks)),
  cbind(data.table(step = "Final cohort: one randomly selected eligible block per patient"),
        count_blocks(cohort$hospital_block_id))
), use.names = TRUE, fill = TRUE)
attrition[, excluded_hospital_blocks_from_previous := NA_integer_]
for (i in 2:nrow(attrition)) {
  if (!is.na(attrition$n_hospital_blocks[i - 1L]) &&
      !is.na(attrition$n_hospital_blocks[i])) {
    attrition$excluded_hospital_blocks_from_previous[i] <-
      attrition$n_hospital_blocks[i - 1L] - attrition$n_hospital_blocks[i]
  }
}

design_summary <- data.table(
  parameter = c(
    "site_assignment_source", "site_source_column", "fallback_site_id",
    "observed_component_site_count", "observed_component_site_ids",
    "cohort_unit", "hospital_block_gap_hours",
    "running_block_end_algorithm",
    "age_minimum", "age_maximum_exclusive", "admission_type_category",
    "pneumonia_case_definition", "qualifying_icd10_families",
    "pneumonia_must_be_present_on_admission", "pneumonia_diagnosis_position",
    "support_devices", "support_start_days", "baseline_pairing_hours",
    "hfno_sf_max_exclusive", "hfno_fio2_min_exclusive",
    "future_support_duration_required", "exclude_death_before_time_zero",
    "exclude_dni_equivalent_at_time_zero",
    "shock_lookback_hours", "lactate_lookback_hours",
    "imv_catecholamine_threshold_mcg_kg_min", "lactate_threshold_mmol_l",
    "active_pressor_definition", "influenza_exclusion_source",
    "missing_shock_data_policy", "one_eligible_block_per_patient",
    "selection_seed", "duplicate_hospitalization_rows_beyond_first",
    "eligible_blocks_before_patient_selection", "final_blocks", "final_patients"
  ),
  value = as.character(c(
    site_assignment_source, site_source_column, fallback_site_id,
    uniqueN(selected_crosswalk$component_site_id, na.rm = TRUE),
    paste(sort(unique(selected_crosswalk$component_site_id[
      !is.na(selected_crosswalk$component_site_id)
    ])), collapse = ";"),
    "hospital_block", block_gap_hours, TRUE,
    age_min, age_max_exclusive, "ed",
    case_definition, "J12.*;J13-J16.*;J18.*", TRUE, "any",
    "imv;nippv;high_flow_nc", support_start_days, baseline_pairing_hours,
    hfno_sf_max_exclusive, hfno_fio2_min_exclusive,
    FALSE, TRUE, TRUE, shock_lookback_hours, lactate_lookback_hours,
    high_dose_pressor_threshold, lactate_exclusion_threshold,
    "most recent continuous-medication record in prior lookback window",
    "present-on-admission ICD-10-CM J09-J11",
    "retain and flag when no exclusion evidence is observed",
    TRUE, selection_seed, duplicate_hospitalization_rows,
    nrow(eligible_blocks), nrow(cohort), uniqueN(cohort$patient_id)
  ))
)

icd_definition <- data.table(
  code_system = "ICD-10-CM",
  qualifying_codes = "J12.*;J13-J16.*;J18.*",
  required_poa_present = TRUE,
  diagnosis_position = "any",
  nonqualifying_or_excluded_codes = paste(
    "J09-J11 influenza;J17 pneumonia in diseases classified elsewhere;",
    "J69 aspiration;J95.851 ventilator-associated pneumonia"
  ),
  note = paste(
    "U07.1 and organism/sepsis codes do not qualify without a qualifying",
    "pneumonia diagnosis; J12.82 is included through J12.*"
  )
)

# ---- Save --------------------------------------------------------------------

output_paths <- c(
  cohort_parquet = file.path(output_dir, "adult_pneumonia_cohort.parquet"),
  cohort_csv = file.path(output_dir, "adult_pneumonia_cohort.csv"),
  eligible_blocks = file.path(output_dir,
    "adult_pneumonia_eligible_blocks_before_patient_selection.parquet"),
  selected_crosswalk = file.path(output_dir,
    "adult_pneumonia_hospitalization_crosswalk.parquet"),
  all_crosswalk = file.path(output_dir, "hospital_block_crosswalk_all.parquet"),
  diagnoses_parquet = file.path(output_dir, "adult_pneumonia_diagnoses.parquet"),
  all_diagnoses_parquet = file.path(output_dir,
    "all_classified_pneumonia_diagnoses.parquet"),
  respiratory_eligibility = file.path(output_dir,
    "adult_pneumonia_respiratory_eligibility.parquet"),
  sf_pairs = file.path(output_dir, "adult_pneumonia_eligible_sf_pairs.parquet"),
  pf_pairs = file.path(output_dir, "adult_pneumonia_eligible_pf_pairs.parquet"),
  attrition_csv = file.path(output_dir, "adult_pneumonia_attrition.csv"),
  design_csv = file.path(output_dir, "adult_pneumonia_cohort_design.csv"),
  icd_definition_csv = file.path(output_dir, "pneumonia_icd_definition.csv")
)

write_parquet_dt(cohort, output_paths[["cohort_parquet"]])
fwrite(cohort, output_paths[["cohort_csv"]], na = "")
write_parquet_dt(eligible_blocks, output_paths[["eligible_blocks"]])
write_parquet_dt(selected_crosswalk, output_paths[["selected_crosswalk"]])
write_parquet_dt(hospital_block_crosswalk_all, output_paths[["all_crosswalk"]])
write_parquet_dt(selected_diagnoses, output_paths[["diagnoses_parquet"]])
write_parquet_dt(diagnoses_classified, output_paths[["all_diagnoses_parquet"]])
write_parquet_dt(respiratory_eligibility,
                 output_paths[["respiratory_eligibility"]])
write_parquet_dt(sf_pairs, output_paths[["sf_pairs"]])
write_parquet_dt(pf_pairs, output_paths[["pf_pairs"]])
fwrite(attrition, output_paths[["attrition_csv"]], na = "")
fwrite(design_summary, output_paths[["design_csv"]], na = "")
fwrite(icd_definition, output_paths[["icd_definition_csv"]], na = "")

message("\nCohort build complete.")
print(attrition)
message("\nFinal cohort contains ", nrow(cohort), " blocks and ",
        uniqueN(cohort$patient_id), " patients.")
message("\nWrote:\n", paste("  ", output_paths, collapse = "\n"))
