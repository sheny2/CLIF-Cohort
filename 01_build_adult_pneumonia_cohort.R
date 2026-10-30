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
#   3. expanded pneumonia diagnosis definition (see pneumonia_icd_helpers.R)
#   4. >=24 continuous hours of IMV, NIPPV, or HFNC >=30 L/min, with the
#      qualifying episode beginning within seven days after block admission
#   5. FiO2 documented and no tracheostomy in the first 24 hours of the episode
#   6. no DNI-equivalent code status in the example Rmd's eligibility window
#   7. eligible P/F or S/F ratio in the first 24 hours: oxygen measurement and
#      valid FiO2 documentation within four hours; S/F additionally requires
#      SpO2 80-97%. P/F uses PaO2 40-700 and excludes ratios >=600.
#
# Usage:
#   Rscript 01_build_adult_pneumonia_cohort.R \
#     CLIF_DATA_DIR [OUTPUT_DIR] [CASE_DEFINITION] [BLOCK_GAP_HOURS] [SEED] \
#     [FALLBACK_SITE_ID] [HOSPITALIZATION_SITE_COLUMN]
#
# CASE_DEFINITION: expanded (default), principal, or any.

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
  file.path(script_dir, "output")

case_definition <- if (length(args) >= 3L && nzchar(args[3L])) args[3L] else
  Sys.getenv("PNEUMONIA_CASE_DEFINITION", unset = "expanded")
case_definition <- tolower(trimws(case_definition))
# Backward compatibility for the former PRIMARY_ONLY Boolean argument.
if (case_definition %in% c("true", "t", "1", "yes", "y")) {
  case_definition <- "principal"
}
if (case_definition %in% c("false", "f", "0", "no", "n")) {
  case_definition <- "any"
}
if (!case_definition %in% c("expanded", "principal", "any")) {
  stop("CASE_DEFINITION must be expanded, principal, or any.", call. = FALSE)
}

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
support_duration_hours <- 24
support_data_days <- support_start_days + support_duration_hours / 24
oxygenation_window_hours <- 24
fio2_pairing_hours <- 4
fio2_min <- 0.21
fio2_max <- 1
spo2_min <- 80
spo2_max_for_sf <- 97
pao2_min <- 40
pao2_max <- 700
pf_max_exclusive <- 600
hfnc_min_lpm <- 30

paths <- c(
  hospitalization = file.path(data_dir, "clif_hospitalization.parquet"),
  diagnosis = file.path(data_dir, "clif_hospital_diagnosis.parquet"),
  respiratory_support = file.path(data_dir, "clif_respiratory_support.parquet"),
  vitals = file.path(data_dir, "clif_vitals.parquet"),
  labs = file.path(data_dir, "clif_labs.parquet"),
  code_status = file.path(data_dir, "clif_code_status.parquet")
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
message("Respiratory support: >=", support_duration_hours,
        " continuous hours starting within ", support_start_days, " days")
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
hospitalization <- as.data.table(read_parquet(paths[["hospitalization"]]))
diagnosis <- as.data.table(read_parquet(paths[["diagnosis"]]))

check_columns(hospitalization, c(
  "hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm",
  "age_at_admission", "admission_type_category"
), "clif_hospitalization.parquet")
check_columns(diagnosis, c(
  "hospitalization_id", "diagnosis_code", "diagnosis_code_format",
  "diagnosis_primary"
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
    message("Deriving encounter site from the first ADT hospital_id...")
    adt_site <- open_dataset(adt_path) |>
      select(hospitalization_id, hospital_id, in_dttm) |>
      filter(!is.na(hospital_id)) |>
      collect() |>
      as.data.table()
    adt_site[, `:=`(
      hospitalization_id = as.character(hospitalization_id),
      hospital_id = as.character(hospital_id),
      in_dttm = as.POSIXct(in_dttm, tz = "UTC")
    )]
    setorder(adt_site, hospitalization_id, in_dttm, na.last = TRUE)
    adt_site <- adt_site[
      !is.na(hospital_id) & nzchar(trimws(hospital_id)),
      .(source_site_id = first_nonmissing(hospital_id)),
      by = hospitalization_id
    ]
    hospitalization <- merge(
      hospitalization, adt_site, by = "hospitalization_id", all.x = TRUE,
      sort = FALSE
    )
    site_assignment_source <- "first_adt_hospital"
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
hospitalization[
  is.na(source_site_id) | !nzchar(trimws(source_site_id)),
  source_site_id := fallback_site_id
]
hospitalization[, source_site_id := trimws(as.character(source_site_id))]
message("Site assignment source: ", site_assignment_source,
        if (nzchar(site_source_column)) paste0(" (", site_source_column, ")") else "")
message("Observed component-encounter sites: ",
        uniqueN(hospitalization$source_site_id, na.rm = TRUE))

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

message("Applying ICD-9-CM/ICD-10-CM pneumonia case definition...")
icd <- derive_pneumonia_encounter_flags(diagnosis)
diagnoses_classified <- icd$diagnoses
diagnosis_encounters <- icd$encounters
diagnosis_encounters[, diagnosis_case := switch(
  case_definition,
  expanded = expanded_pneumonia_case,
  principal = principal_pneumonia,
  any = any_pneumonia
)]

diagnosis_encounters <- merge(
  diagnosis_encounters,
  hospital_block_crosswalk_all[, .(
    patient_id, hospital_block_id, hospitalization_id
  )],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
qualifying_encounters <- diagnosis_encounters[diagnosis_case == 1L]

diagnoses_classified <- merge(
  diagnoses_classified,
  hospital_block_crosswalk_all[, .(
    patient_id, hospital_block_id, hospitalization_id
  )],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
diagnoses_classified <- diagnoses_classified[
  is_pneumonia_code == 1L | is_sepsis_code == 1L |
    is_respiratory_failure_code == 1L
]

diagnosis_summary_by_block <- qualifying_encounters[, .(
  pneumonia_primary = as.integer(any(principal_pneumonia == 1L)),
  pneumonia_case_hospitalizations = uniqueN(hospitalization_id),
  pneumonia_case_pathways = paste(sort(unique(pneumonia_case_pathway)),
                                   collapse = ";")
), by = hospital_block_id]
pneumonia_code_summary <- diagnoses_classified[is_pneumonia_code == 1L, .(
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

# ---- Respiratory-support episodes -------------------------------------------

message("Constructing qualifying respiratory-support episodes...")
respiratory <- read_filtered_clif(
  paths[["respiratory_support"]],
  "hospitalization_id",
  unique(candidate_crosswalk$hospitalization_id),
  c(
    "hospitalization_id", "recorded_dttm", "device_category", "fio2_set",
    "lpm_set", "tracheostomy"
  )
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
resp[
  device_category_clean == "nasal_cannula" & lpm >= 20 & fio2 >= 0.30,
  device_category_clean := "high_flow_nc"
]
resp <- resp[
  !is.na(recorded_dttm) & recorded_dttm >= block_admission_dttm &
    recorded_dttm <= support_window_end
]
setorder(resp, hospital_block_id, device_category_clean, recorded_dttm)
resp[, `:=`(
  lpm_carried = nafill(nafill(lpm, type = "locf"), type = "nocb"),
  fio2_carried = nafill(nafill(fio2, type = "locf"), type = "nocb")
), by = .(hospital_block_id, device_category_clean)]
resp[, qualifying_support := as.integer(
  device_category_clean %in% c("imv", "nippv") |
    (device_category_clean == "high_flow_nc" & lpm_carried >= hfnc_min_lpm)
)]
resp[, support_rank := fcase(
  device_category_clean == "imv", 3L,
  device_category_clean == "nippv", 2L,
  device_category_clean == "high_flow_nc", 1L,
  default = 0L
)]

# Resolve simultaneous device rows to the highest support state.
support_state <- resp[
  order(-qualifying_support, -support_rank), .SD[1L],
  by = .(hospital_block_id, recorded_dttm)
]
setorder(support_state, hospital_block_id, recorded_dttm)
support_state[, previous_support := shift(qualifying_support, fill = 0L),
              by = hospital_block_id]
support_state[, support_episode_number := cumsum(
  qualifying_support == 1L & previous_support == 0L
), by = hospital_block_id]
support_state[, next_dttm := shift(recorded_dttm, type = "lead"),
              by = hospital_block_id]
support_state[is.na(next_dttm), next_dttm := support_window_end]
support_state[, interval_hours := pmax(
  0, as.numeric(difftime(next_dttm, recorded_dttm, units = "hours"))
)]

support_episodes <- support_state[
  qualifying_support == 1L,
  .(
    respiratory_failure_start = min(recorded_dttm),
    respiratory_support_duration_hours = sum(interval_hours),
    respiratory_failure_stop = max(next_dttm),
    episode_support_devices = paste(sort(unique(device_category_clean)),
                                    collapse = ";")
  ),
  by = .(hospital_block_id, support_episode_number)
]
support_episodes <- merge(
  support_episodes,
  diagnosis_candidate_blocks[, .(hospital_block_id, admission_dttm)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
support_episodes <- support_episodes[
  respiratory_support_duration_hours >= support_duration_hours &
    respiratory_failure_start <= admission_dttm + support_start_days * 24 * 60 * 60
]
setorder(support_episodes, hospital_block_id, respiratory_failure_start)
qualifying_episode <- support_episodes[, .SD[1L], by = hospital_block_id]
qualifying_episode[, admission_dttm := NULL]
if (!nrow(qualifying_episode)) {
  stop("No blocks have a qualifying >=24-hour respiratory-support episode.")
}

# ---- First-24-hour respiratory exclusions and oxygenation -------------------

message("Applying first-24-hour FiO2, tracheostomy, code-status, and P/F-S/F criteria...")
resp_first24 <- merge(
  resp, qualifying_episode[, .(hospital_block_id, respiratory_failure_start)],
  by = "hospital_block_id", all = FALSE, sort = FALSE
)
resp_first24 <- resp_first24[
  recorded_dttm >= respiratory_failure_start &
    recorded_dttm <= respiratory_failure_start + oxygenation_window_hours * 60 * 60
]
fio2_records <- unique(resp_first24[
  !is.na(fio2), .(hospital_block_id, fio2_dttm = recorded_dttm, fio2)
])
first24_support_summary <- resp_first24[, .(
  fio2_documented_first24 = as.integer(any(!is.na(fio2))),
  tracheostomy_first24 = as.integer(any(tracheostomy_binary == 1L)),
  maximum_support_rank_first24 = max(support_rank, na.rm = TRUE)
), by = hospital_block_id]
first24_support_summary[, respiratory_support_category := fcase(
  maximum_support_rank_first24 >= 3L, "imv",
  maximum_support_rank_first24 == 2L, "nippv",
  maximum_support_rank_first24 == 1L, "high_flow_nc",
  default = NA_character_
)]

episode_windows <- merge(
  qualifying_episode[, .(hospital_block_id, respiratory_failure_start)],
  diagnosis_candidate_blocks[, .(hospital_block_id, patient_id)],
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)

episode_crosswalk <- candidate_crosswalk[
  hospital_block_id %in% qualifying_episode$hospital_block_id
]
vitals <- read_filtered_clif(
  paths[["vitals"]],
  "hospitalization_id",
  unique(episode_crosswalk$hospitalization_id),
  c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value")
)
vitals[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  recorded_dttm = as.POSIXct(recorded_dttm, tz = "UTC")
)]
labs <- read_filtered_clif(
  paths[["labs"]],
  "hospitalization_id",
  unique(episode_crosswalk$hospitalization_id),
  c(
    "hospitalization_id", "lab_collect_dttm", "lab_category",
    "lab_value_numeric"
  )
)
labs[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  lab_collect_dttm = as.POSIXct(lab_collect_dttm, tz = "UTC")
)]

spo2 <- merge(
  vitals[normalize_category(vital_category) == "spo2", .(
    hospitalization_id,
    measurement_dttm = recorded_dttm,
    measurement_value = suppressWarnings(as.numeric(vital_value))
  )],
  candidate_crosswalk[, .(hospitalization_id, hospital_block_id)],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
spo2 <- merge(spo2, qualifying_episode[, .(
  hospital_block_id, respiratory_failure_start
)], by = "hospital_block_id", all = FALSE, sort = FALSE)
spo2 <- spo2[
  measurement_dttm >= respiratory_failure_start &
    measurement_dttm <= respiratory_failure_start + oxygenation_window_hours * 60 * 60 &
    measurement_value >= spo2_min & measurement_value <= spo2_max_for_sf
]
sf_pairs <- pair_to_nearest_fio2(
  spo2[, .(hospital_block_id, measurement_dttm, measurement_value)],
  fio2_records, "spo2"
)
if (nrow(sf_pairs)) sf_pairs[, sf_ratio := spo2 / fio2]

pao2 <- merge(
  labs[normalize_category(lab_category) == "po2_arterial", .(
    hospitalization_id,
    measurement_dttm = lab_collect_dttm,
    measurement_value = suppressWarnings(as.numeric(lab_value_numeric))
  )],
  candidate_crosswalk[, .(hospitalization_id, hospital_block_id)],
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
pao2 <- merge(pao2, qualifying_episode[, .(
  hospital_block_id, respiratory_failure_start
)], by = "hospital_block_id", all = FALSE, sort = FALSE)
pao2 <- pao2[
  measurement_dttm >= respiratory_failure_start &
    measurement_dttm <= respiratory_failure_start + oxygenation_window_hours * 60 * 60 &
    measurement_value >= pao2_min & measurement_value <= pao2_max
]
pf_pairs <- pair_to_nearest_fio2(
  pao2[, .(hospital_block_id, measurement_dttm, measurement_value)],
  fio2_records, "pao2"
)
if (nrow(pf_pairs)) {
  pf_pairs[, pf_ratio := pao2 / fio2]
  pf_pairs <- pf_pairs[pf_ratio < pf_max_exclusive]
}

sf_summary <- if (nrow(sf_pairs)) sf_pairs[, .(
  eligible_sf_measurements = .N,
  sf_24_min = min(sf_ratio),
  min_sf_spo2 = spo2[which.min(sf_ratio)],
  min_sf_fio2 = fio2[which.min(sf_ratio)],
  min_sf_dttm = measurement_dttm[which.min(sf_ratio)]
), by = hospital_block_id] else data.table(hospital_block_id = character())
pf_summary <- if (nrow(pf_pairs)) pf_pairs[, .(
  eligible_pf_measurements = .N,
  pf_24_min = min(pf_ratio),
  min_pf_pao2 = pao2[which.min(pf_ratio)],
  min_pf_fio2 = fio2[which.min(pf_ratio)],
  min_pf_dttm = measurement_dttm[which.min(pf_ratio)]
), by = hospital_block_id] else data.table(hospital_block_id = character())

respiratory_eligibility <- Reduce(
  function(x, y) merge(x, y, by = "hospital_block_id", all.x = TRUE,
                       sort = FALSE),
  list(qualifying_episode, first24_support_summary, sf_summary, pf_summary)
)
for (column in c("eligible_sf_measurements", "eligible_pf_measurements")) {
  if (!column %in% names(respiratory_eligibility)) {
    respiratory_eligibility[, (column) := 0L]
  }
  set(respiratory_eligibility, which(is.na(respiratory_eligibility[[column]])),
      column, 0L)
}
respiratory_eligibility[, eligible_oxygenation := as.integer(
  eligible_sf_measurements > 0L | eligible_pf_measurements > 0L
)]

# Last DNI-equivalent status from seven days before through 24 hours after start.
code_status <- read_filtered_clif(
  paths[["code_status"]],
  "patient_id",
  unique(episode_windows$patient_id),
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
    start_dttm <= respiratory_failure_start + 24 * 60 * 60
]
setorder(status_candidates, hospital_block_id, start_dttm)
status_summary <- status_candidates[, .(
  code_status_category_eligibility = as.character(last_nonmissing(
    code_status_category
  ))
), by = hospital_block_id]
status_summary[, limited_code_status := as.integer(
  normalize_category(code_status_category_eligibility) %in%
    c("and", "dnr_dni", "dnar_dni")
)]
respiratory_eligibility <- merge(
  respiratory_eligibility, status_summary,
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
respiratory_eligibility[is.na(limited_code_status), limited_code_status := 0L]
respiratory_eligibility[, respiratory_eligible := as.integer(
  fio2_documented_first24 == 1L &
    tracheostomy_first24 == 0L &
    limited_code_status == 0L &
    eligible_oxygenation == 1L
)]

# ---- Final eligibility and one block per patient ----------------------------

eligible_blocks <- merge(
  diagnosis_candidate_blocks,
  respiratory_eligibility[respiratory_eligible == 1L],
  by = "hospital_block_id", all = FALSE, sort = FALSE
)
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
fio2_blocks <- respiratory_eligibility[fio2_documented_first24 == 1L,
                                       hospital_block_id]
no_trach_blocks <- respiratory_eligibility[
  fio2_documented_first24 == 1L & tracheostomy_first24 == 0L,
  hospital_block_id
]
code_status_blocks <- respiratory_eligibility[
  fio2_documented_first24 == 1L & tracheostomy_first24 == 0L &
    limited_code_status == 0L,
  hospital_block_id
]
oxygenation_blocks <- respiratory_eligibility[
  fio2_documented_first24 == 1L & tracheostomy_first24 == 0L &
    limited_code_status == 0L & eligible_oxygenation == 1L,
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
  cbind(data.table(step = paste0("ED adult blocks meeting ", case_definition,
                                 " pneumonia coding definition")),
        count_blocks(diagnosis_candidate_blocks$hospital_block_id)),
  cbind(data.table(step = "Blocks with qualifying >=24-hour respiratory-support episode"),
        count_blocks(episode_blocks)),
  cbind(data.table(step = "Blocks with FiO2 documented in first episode day"),
        count_blocks(fio2_blocks)),
  cbind(data.table(step = "Blocks without tracheostomy in first episode day"),
        count_blocks(no_trach_blocks)),
  cbind(data.table(step = "Blocks without DNI-equivalent status"),
        count_blocks(code_status_blocks)),
  cbind(data.table(step = "Blocks with eligible P/F or S/F measurement"),
        count_blocks(oxygenation_blocks)),
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
    "pneumonia_case_definition", "support_devices", "hfnc_min_lpm",
    "support_duration_hours", "support_start_days", "oxygenation_window_hours",
    "fio2_pairing_hours", "fio2_range", "spo2_range_for_sf",
    "pao2_range_for_pf", "pf_max_exclusive", "exclude_tracheostomy_first24",
    "exclude_dni_equivalent", "one_eligible_block_per_patient",
    "selection_seed", "duplicate_hospitalization_rows_beyond_first",
    "eligible_blocks_before_patient_selection", "final_blocks", "final_patients"
  ),
  value = as.character(c(
    site_assignment_source, site_source_column, fallback_site_id,
    uniqueN(hospitalization$source_site_id, na.rm = TRUE),
    paste(sort(unique(hospitalization$source_site_id)), collapse = ";"),
    "hospital_block", block_gap_hours, TRUE,
    age_min, age_max_exclusive, "ed",
    case_definition, "imv;nippv;high_flow_nc", hfnc_min_lpm,
    support_duration_hours, support_start_days, oxygenation_window_hours,
    fio2_pairing_hours, "0.21-1.00", "80-97", "40-700", pf_max_exclusive,
    TRUE, TRUE, TRUE, selection_seed, duplicate_hospitalization_rows,
    nrow(eligible_blocks), nrow(cohort), uniqueN(cohort$patient_id)
  ))
)

icd_definition <- data.table(
  code_system = c("ICD-9-CM", "ICD-9-CM", "ICD-9-CM",
                  "ICD-10-CM", "ICD-10-CM", "ICD-10-CM"),
  concept = rep(c("pneumonia", "sepsis", "respiratory failure"), 2),
  included_codes = c(
    "481;482;483;485;486",
    "038;995.91;995.92;785.52",
    "518.81;518.82;518.84;799.1",
    "J13-J16;J18",
    "A40;A41;R65.20;R65.21",
    "J80;J96.0;J96.2;J96.9;R09.2"
  ),
  role_in_default_definition = c(
    "principal, or secondary with principal sepsis/respiratory failure",
    "principal with secondary pneumonia",
    "principal with secondary pneumonia",
    "principal, or secondary with principal sepsis/respiratory failure",
    "principal with secondary pneumonia",
    "principal with secondary pneumonia"
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
    "all_classified_pneumonia_sepsis_respiratory_diagnoses.parquet"),
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
