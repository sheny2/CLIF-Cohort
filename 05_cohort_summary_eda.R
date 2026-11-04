#!/usr/bin/env Rscript

# Create disclosure-conscious aggregate summaries and figures for the Hopkins
# CLIF pneumonia/steroid cohort. All patient-level inputs remain on Discovery
# HPC. This script writes CSV summaries and PNG figures only; it never exports
# patient, hospitalization, or hospital-block identifiers.
#
# Usage:
#   Rscript 05_cohort_summary_eda.R \
#     CLIF_DATA_DIR [COHORT_DIR] [TREATMENT_DIR] [OUTPUT_DIR] [MIN_CELL_N] \
#     [DATA_CUTOFF_DATE]

rm(list = ls())
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else normalizePath(getwd())

data_dir <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else
  Sys.getenv("CLIF_DATA_DIR", unset = "")
cohort_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "output")
treatment_dir <- if (length(args) >= 3L && nzchar(args[3L])) args[3L] else
  file.path(script_dir, "steroid_treatment")
output_dir <- if (length(args) >= 4L && nzchar(args[4L])) args[4L] else
  file.path(script_dir, "summary_eda")
min_cell_n <- if (length(args) >= 5L && nzchar(args[5L])) {
  suppressWarnings(as.integer(args[5L]))
} else suppressWarnings(as.integer(
  Sys.getenv("SUMMARY_SMALL_CELL_THRESHOLD", unset = "10")
))
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Supply a valid CLIF_DATA_DIR as argument 1 or an environment variable.",
       call. = FALSE)
}
if (!is.finite(min_cell_n) || min_cell_n < 1L) {
  stop("MIN_CELL_N must be a positive integer.", call. = FALSE)
}
cutoff_text <- if (length(args) >= 6L && nzchar(args[6L])) args[6L] else
  Sys.getenv("CLIF_DATA_CUTOFF_DATE", unset = "")
specified_cutoff_dttm <- if (nzchar(cutoff_text)) {
  as.POSIXct(paste0(as.Date(cutoff_text), " 23:59:59"), tz = "UTC")
} else as.POSIXct(NA, tz = "UTC")
if (nzchar(cutoff_text) && is.na(specified_cutoff_dttm)) {
  stop("DATA_CUTOFF_DATE must be parseable as YYYY-MM-DD.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohort_path <- file.path(cohort_dir, "adult_pneumonia_cohort.parquet")
treatment_path <- file.path(treatment_dir, "steroid_treatment_assignment.parquet")
attrition_path <- file.path(cohort_dir, "adult_pneumonia_attrition.csv")
if (!file.exists(cohort_path)) stop("Missing: ", cohort_path, call. = FALSE)
if (!file.exists(treatment_path)) stop("Missing: ", treatment_path, call. = FALSE)

normalize_text <- function(x) tolower(trimws(as.character(x)))
safe_time <- function(x) as.POSIXct(x, tz = "UTC")
safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probability, na.rm = TRUE, names = FALSE))
}
report_count <- function(x) {
  x <- as.integer(x)
  fifelse(x == 0L, 0L, fifelse(x < min_cell_n, NA_integer_, x))
}
report_percent <- function(numerator, denominator) {
  ifelse(
    denominator < min_cell_n | (numerator > 0 & numerator < min_cell_n),
    NA_real_, 100 * numerator / denominator
  )
}
write_summary <- function(x, filename) {
  forbidden <- intersect(
    c("patient_id", "hospitalization_id", "hospital_block_id"), names(x)
  )
  if (length(forbidden)) {
    stop("Refusing to export identifier column(s) in ", filename, ": ",
         paste(forbidden, collapse = ", "), call. = FALSE)
  }
  fwrite(x, file.path(output_dir, filename), na = "")
}

message("Reading cohort and treatment assignment on HPC...")
cohort <- as.data.table(read_parquet(cohort_path))
treatment <- as.data.table(read_parquet(treatment_path))
required_cohort <- c("patient_id", "hospital_block_id", "respiratory_failure_start")
missing_cohort <- setdiff(required_cohort, names(cohort))
if (length(missing_cohort)) {
  stop("Cohort is missing: ", paste(missing_cohort, collapse = ", "),
       call. = FALSE)
}
if (anyDuplicated(cohort$patient_id) || anyDuplicated(cohort$hospital_block_id)) {
  stop("Final cohort is not one row per patient and hospital block.",
       call. = FALSE)
}
cohort[, `:=`(
  patient_id = as.character(patient_id),
  hospital_block_id = as.character(hospital_block_id),
  respiratory_failure_start = safe_time(respiratory_failure_start)
)]
for (column in intersect(c("admission_dttm", "discharge_dttm"), names(cohort))) {
  set(cohort, j = column, value = safe_time(cohort[[column]]))
}

treatment_keys <- c("patient_id", "hospital_block_id")
if (!all(treatment_keys %in% names(treatment))) {
  stop("Treatment assignment lacks patient_id/hospital_block_id.",
       call. = FALSE)
}
treatment_keep <- intersect(
  c(treatment_keys, "treated_within_grace", "treatment_grace_hours",
    "treatment_admin_dttm", "treatment_hours_from_time_zero",
    "treatment_med_category_clean", "treatment_dose_mg",
    "treatment_med_route_clean", "ever_qualifying_steroid"),
  names(treatment)
)
treatment <- unique(treatment[, ..treatment_keep])
cohort <- merge(cohort, treatment, by = treatment_keys, all.x = TRUE, sort = FALSE)
if (!"treated_within_grace" %in% names(cohort)) {
  cohort[, treated_within_grace := NA_integer_]
}
if (!"ever_qualifying_steroid" %in% names(cohort)) {
  cohort[, ever_qualifying_steroid := NA_integer_]
}
if ("treatment_admin_dttm" %in% names(cohort)) {
  cohort[, treatment_admin_dttm := safe_time(treatment_admin_dttm)]
}
grace_values <- unique(na.omit(cohort$treatment_grace_hours))
grace_hours <- if (length(grace_values)) as.numeric(grace_values[1L]) else NA_real_

# Add patient demographics and death timestamp, but only for cohort patients.
patient_path <- file.path(data_dir, "clif_patient.parquet")
death_source_available <- FALSE
if (file.exists(patient_path)) {
  patient_dataset <- open_dataset(patient_path)
  patient_columns <- patient_dataset$schema$names
  patient_keep <- intersect(
    c("patient_id", "death_dttm", "sex_category", "race_category",
      "ethnicity_category", "language_category"),
    patient_columns
  )
  if ("patient_id" %in% patient_keep) {
    patient <- patient_dataset |>
      filter(.data$patient_id %in% cohort$patient_id) |>
      select(all_of(patient_keep)) |>
      collect() |>
      as.data.table()
    patient[, patient_id := as.character(patient_id)]
    patient <- unique(patient, by = "patient_id")
    if ("death_dttm" %in% names(patient)) {
      patient[, death_dttm := safe_time(death_dttm)]
      death_source_available <- TRUE
    }
    cohort <- merge(cohort, patient, by = "patient_id", all.x = TRUE, sort = FALSE)
  }
}

# Derive follow-up and outcomes. Patient death_dttm is preferred; an Expired
# discharge timestamp is used only when patient death_dttm is unavailable.
cohort[, event_dttm := as.POSIXct(NA, tz = "UTC")]
cohort[, event_source := NA_character_]
if ("death_dttm" %in% names(cohort)) {
  cohort[!is.na(death_dttm), `:=`(
    event_dttm = death_dttm,
    event_source = "patient death_dttm"
  )]
}
if (all(c("discharge_category", "discharge_dttm") %in% names(cohort))) {
  cohort[
    is.na(event_dttm) & normalize_text(discharge_category) == "expired" &
      !is.na(discharge_dttm),
    `:=`(event_dttm = discharge_dttm,
         event_source = "expired discharge fallback")
  ]
}

# A fixed data cutoff is preferred. If it is not supplied, use the latest
# observed admission, discharge, death, or treatment timestamp as a transparent
# empirical proxy and flag that choice in the exported definition.
if (!is.na(specified_cutoff_dttm)) {
  data_cutoff_dttm <- specified_cutoff_dttm
  data_cutoff_source <- "user-specified CLIF_DATA_CUTOFF_DATE"
} else {
  cutoff_candidates <- c(
    if ("admission_dttm" %in% names(cohort)) cohort$admission_dttm else NULL,
    if ("discharge_dttm" %in% names(cohort)) cohort$discharge_dttm else NULL,
    cohort$event_dttm,
    if ("treatment_admin_dttm" %in% names(cohort))
      cohort$treatment_admin_dttm else NULL
  )
  cutoff_candidates <- cutoff_candidates[!is.na(cutoff_candidates)]
  data_cutoff_dttm <- if (length(cutoff_candidates)) {
    max(cutoff_candidates)
  } else as.POSIXct(NA, tz = "UTC")
  data_cutoff_source <- "inferred latest observed timestamp; confirm before analysis"
}
cohort[, event_hours := as.numeric(difftime(
  event_dttm, respiratory_failure_start, units = "hours"
))]
cohort[, complete_28d_followup := as.integer(
  !is.na(data_cutoff_dttm) &
    respiratory_failure_start + 28 * 24 * 60 * 60 <= data_cutoff_dttm
)]
cohort[, complete_grace_followup := if (is.finite(grace_hours)) {
  as.integer(!is.na(data_cutoff_dttm) &
    respiratory_failure_start + grace_hours * 60 * 60 <= data_cutoff_dttm)
} else NA_integer_]
cohort[, death_within_28d := as.integer(
  !is.na(event_hours) & event_hours >= 0 & event_hours <= 28 * 24
)]
cohort[, death_within_grace := if (is.finite(grace_hours)) {
  as.integer(!is.na(event_hours) & event_hours >= 0 & event_hours <= grace_hours)
} else NA_integer_]
cohort[, in_hospital_death := as.integer(
  !is.na(discharge_category) & normalize_text(discharge_category) == "expired"
)]
cohort[, hospice_discharge := as.integer(
  !is.na(discharge_category) & grepl("hospice", normalize_text(discharge_category))
)]
if (all(c("admission_dttm", "discharge_dttm") %in% names(cohort))) {
  cohort[, length_of_stay_days := as.numeric(difftime(
    discharge_dttm, admission_dttm, units = "days"
  ))]
  cohort[, hours_admission_to_time_zero := as.numeric(difftime(
    respiratory_failure_start, admission_dttm, units = "hours"
  ))]
}

# ---- Overall and site-level summaries --------------------------------------

site_column <- if ("source_site_id" %in% names(cohort)) "source_site_id" else NULL
if (is.null(site_column)) cohort[, source_site_id := "site unavailable"]
cohort[, source_site_id := fifelse(
  is.na(source_site_id) | !nzchar(source_site_id), "Missing", source_site_id
)]

overall_summary <- data.table(
  parameter = c(
    "final_cohort_patients", "observed_sites", "admission_date_min",
    "admission_date_max", "data_cutoff", "data_cutoff_source",
    "complete_grace_followup_n", "complete_28d_followup_n",
    "age_median", "age_q1", "age_q3",
    "length_of_stay_median_days", "length_of_stay_q1_days",
    "length_of_stay_q3_days", "treated_within_grace_n",
    "treated_within_grace_percent", "treatment_grace_hours",
    "death_within_grace_n", "death_within_grace_percent",
    "death_within_28d_n", "death_within_28d_percent",
    "in_hospital_death_n", "in_hospital_death_percent",
    "hospice_discharge_n", "hospice_discharge_percent",
    "patient_death_dttm_available", "minimum_reportable_cell"
  ),
  value = as.character(c(
    nrow(cohort), uniqueN(cohort$source_site_id),
    if ("admission_dttm" %in% names(cohort))
      format(min(cohort$admission_dttm, na.rm = TRUE), "%Y-%m-%d") else NA,
    if ("admission_dttm" %in% names(cohort))
      format(max(cohort$admission_dttm, na.rm = TRUE), "%Y-%m-%d") else NA,
    if (!is.na(data_cutoff_dttm)) format(data_cutoff_dttm, "%Y-%m-%d") else NA,
    data_cutoff_source,
    if (is.finite(grace_hours)) sum(cohort$complete_grace_followup) else NA,
    sum(cohort$complete_28d_followup),
    if ("age_at_admission" %in% names(cohort))
      median(cohort$age_at_admission, na.rm = TRUE) else NA,
    if ("age_at_admission" %in% names(cohort))
      safe_quantile(cohort$age_at_admission, 0.25) else NA,
    if ("age_at_admission" %in% names(cohort))
      safe_quantile(cohort$age_at_admission, 0.75) else NA,
    if ("length_of_stay_days" %in% names(cohort))
      median(cohort$length_of_stay_days, na.rm = TRUE) else NA,
    if ("length_of_stay_days" %in% names(cohort))
      safe_quantile(cohort$length_of_stay_days, 0.25) else NA,
    if ("length_of_stay_days" %in% names(cohort))
      safe_quantile(cohort$length_of_stay_days, 0.75) else NA,
    if (all(!is.na(cohort$treated_within_grace)) && is.finite(grace_hours))
      sum(cohort$treated_within_grace[cohort$complete_grace_followup == 1L]) else NA,
    if (all(!is.na(cohort$treated_within_grace)) && is.finite(grace_hours) &&
        sum(cohort$complete_grace_followup) > 0)
      100 * mean(cohort$treated_within_grace[
        cohort$complete_grace_followup == 1L
      ]) else NA,
    grace_hours,
    if (is.finite(grace_hours)) sum(cohort$death_within_grace[
      cohort$complete_grace_followup == 1L
    ]) else NA,
    if (is.finite(grace_hours) && sum(cohort$complete_grace_followup) > 0)
      100 * mean(cohort$death_within_grace[
        cohort$complete_grace_followup == 1L
      ]) else NA,
    sum(cohort$death_within_28d[cohort$complete_28d_followup == 1L]),
    if (sum(cohort$complete_28d_followup) > 0) 100 * mean(
      cohort$death_within_28d[cohort$complete_28d_followup == 1L]
    ) else NA,
    sum(cohort$in_hospital_death), 100 * mean(cohort$in_hospital_death),
    sum(cohort$hospice_discharge), 100 * mean(cohort$hospice_discharge),
    death_source_available, min_cell_n
  ))
)
write_summary(overall_summary, "cohort_overall_summary.csv")

site_summary <- cohort[, .(
  n = .N,
  n_grace_complete = sum(complete_grace_followup, na.rm = TRUE),
  n_28d_complete = sum(complete_28d_followup),
  treated = sum(treated_within_grace[complete_grace_followup == 1L], na.rm = TRUE),
  death_grace = sum(death_within_grace[complete_grace_followup == 1L],
                    na.rm = TRUE),
  death_28d = sum(death_within_28d[complete_28d_followup == 1L]),
  in_hospital_death = sum(in_hospital_death),
  hospice = sum(hospice_discharge),
  age_median = median(age_at_admission, na.rm = TRUE),
  los_median_days = if ("length_of_stay_days" %in% names(cohort))
    median(length_of_stay_days, na.rm = TRUE) else NA_real_
), by = source_site_id]
site_summary[, `:=`(
  n_reportable = report_count(n),
  n_grace_complete_reportable = report_count(n_grace_complete),
  n_28d_complete_reportable = report_count(n_28d_complete),
  treated_reportable = report_count(treated),
  treated_percent = report_percent(treated, n_grace_complete),
  death_within_grace_reportable = report_count(death_grace),
  death_within_grace_percent = report_percent(death_grace, n_grace_complete),
  death_28d_reportable = report_count(death_28d),
  death_28d_percent = report_percent(death_28d, n_28d_complete),
  in_hospital_death_reportable = report_count(in_hospital_death),
  in_hospital_death_percent = report_percent(in_hospital_death, n),
  hospice_reportable = report_count(hospice),
  hospice_percent = report_percent(hospice, n),
  age_median = fifelse(n >= min_cell_n, age_median, NA_real_),
  los_median_days = fifelse(n >= min_cell_n, los_median_days, NA_real_),
  small_cell_threshold = min_cell_n
)]
site_summary[, c("n", "n_grace_complete", "n_28d_complete",
                 "treated", "death_grace", "death_28d",
                 "in_hospital_death", "hospice") := NULL]
setorder(site_summary, source_site_id)
write_summary(site_summary, "cohort_by_site_suppressed.csv")

# ---- Treatment, outcome, baseline, and missingness diagnostics -------------

treatment_group_available <- all(!is.na(cohort$treated_within_grace))
if (treatment_group_available) {
  cohort[, treatment_group := fcase(
    complete_grace_followup != 1L, "Incomplete grace-period follow-up",
    treated_within_grace == 1L, "Treated within grace",
    default = "Not treated within grace"
  )]
} else cohort[, treatment_group := "Grace-period assignment unavailable"]

outcome_summary <- cohort[, .(
  n = .N,
  n_grace_complete = sum(complete_grace_followup, na.rm = TRUE),
  n_28d_complete = sum(complete_28d_followup),
  death_grace = sum(death_within_grace[complete_grace_followup == 1L],
                    na.rm = TRUE),
  death_28d = sum(death_within_28d[complete_28d_followup == 1L]),
  in_hospital_death = sum(in_hospital_death),
  hospice = sum(hospice_discharge)
), by = treatment_group]
for (event in c("death_grace", "death_28d", "in_hospital_death", "hospice")) {
  outcome_summary[, paste0(event, "_reportable") := report_count(get(event))]
  denominator_column <- if (event == "death_grace") {
    "n_grace_complete"
  } else if (event == "death_28d") {
    "n_28d_complete"
  } else "n"
  outcome_summary[, paste0(event, "_percent") := report_percent(
    get(event), get(denominator_column)
  )]
}
outcome_summary[, n_reportable := report_count(n)]
outcome_summary[, n_grace_complete_reportable := report_count(n_grace_complete)]
outcome_summary[, n_28d_complete_reportable := report_count(n_28d_complete)]
outcome_summary[, c("n", "n_grace_complete", "n_28d_complete",
                    "death_grace", "death_28d",
                    "in_hospital_death", "hospice") := NULL]
outcome_summary[, small_cell_threshold := min_cell_n]
write_summary(outcome_summary, "mortality_by_treatment_suppressed.csv")

continuous_candidates <- intersect(c(
  "age_at_admission", "length_of_stay_days", "hours_admission_to_time_zero",
  "respiratory_support_duration_hours", "baseline_sf_ratio",
  "baseline_fio2", "baseline_spo2", "lactate_value", "weight_kg",
  "max_catecholamine_dose_mcg_kg_min"
), names(cohort))
continuous_parts <- lapply(continuous_candidates, function(variable) {
  cohort[, {
    value <- suppressWarnings(as.numeric(get(variable)))
    observed <- value[is.finite(value)]
    n_observed <- length(observed)
    report <- n_observed >= min_cell_n
    .(
      variable = variable,
      level = "Continuous",
      n_nonmissing_reportable = report_count(n_observed),
      missing_percent = 100 * mean(!is.finite(value)),
      mean = if (report) mean(observed) else NA_real_,
      sd = if (report) sd(observed) else NA_real_,
      median = if (report) median(observed) else NA_real_,
      q1 = if (report) safe_quantile(observed, 0.25) else NA_real_,
      q3 = if (report) safe_quantile(observed, 0.75) else NA_real_
    )
  }, by = treatment_group]
})

categorical_candidates <- intersect(c(
  "source_site_id", "sex_category", "race_category", "ethnicity_category",
  "respiratory_support_category", "discharge_category"
), names(cohort))
categorical_parts <- lapply(categorical_candidates, function(variable) {
  x <- cohort[, .(
    level = fifelse(is.na(get(variable)) | !nzchar(as.character(get(variable))),
                    "Missing", as.character(get(variable))),
    treatment_group
  )]
  denominator <- x[, .(denominator = .N), by = treatment_group]
  out <- x[, .(n = .N), by = .(treatment_group, level)]
  out <- merge(out, denominator, by = "treatment_group", all.x = TRUE)
  out[, `:=`(
    variable = variable,
    n_reportable = report_count(n),
    percent = report_percent(n, denominator)
  )]
  out[, c("n", "denominator") := NULL]
  out
})
baseline_continuous <- rbindlist(continuous_parts, fill = TRUE)
baseline_categorical <- rbindlist(categorical_parts, fill = TRUE)
if (nrow(baseline_continuous)) {
  baseline_continuous[, small_cell_threshold := min_cell_n]
  write_summary(baseline_continuous,
                "baseline_continuous_by_treatment_suppressed.csv")
}
if (nrow(baseline_categorical)) {
  baseline_categorical[, small_cell_threshold := min_cell_n]
  write_summary(baseline_categorical,
                "baseline_categorical_by_treatment_suppressed.csv")
}

missingness_variables <- unique(c(
  continuous_candidates, categorical_candidates,
  intersect(c("treated_within_grace", "death_dttm"), names(cohort))
))
missingness_summary <- rbindlist(lapply(missingness_variables, function(variable) {
  value <- cohort[[variable]]
  missing <- is.na(value) | (is.character(value) & !nzchar(value))
  data.table(
    variable = variable,
    n_missing_reportable = report_count(sum(missing)),
    missing_percent = 100 * mean(missing),
    n_nonmissing_reportable = report_count(sum(!missing)),
    small_cell_threshold = min_cell_n
  )
}))
setorder(missingness_summary, -missing_percent)
write_summary(missingness_summary, "missingness_summary_suppressed.csv")

if ("treatment_hours_from_time_zero" %in% names(cohort)) {
  treatment_times <- cohort[
    treated_within_grace == 1L & is.finite(treatment_hours_from_time_zero),
    treatment_hours_from_time_zero
  ]
  treatment_timing <- data.table(
    measure = c("n", "minimum", "q1", "median", "q3", "maximum"),
    value = c(
      if (length(treatment_times) >= min_cell_n) length(treatment_times) else NA,
      if (length(treatment_times) >= min_cell_n) min(treatment_times) else NA,
      if (length(treatment_times) >= min_cell_n)
        safe_quantile(treatment_times, 0.25) else NA,
      if (length(treatment_times) >= min_cell_n) median(treatment_times) else NA,
      if (length(treatment_times) >= min_cell_n)
        safe_quantile(treatment_times, 0.75) else NA,
      if (length(treatment_times) >= min_cell_n) max(treatment_times) else NA
    ),
    units = c("patients", rep("hours from time zero", 5)),
    small_cell_threshold = min_cell_n
  )
  write_summary(treatment_timing, "treatment_timing_summary.csv")
}

discharge_summary <- cohort[, .(
  discharge_category = fifelse(
    is.na(discharge_category) | !nzchar(as.character(discharge_category)),
    "Missing", as.character(discharge_category)
  )
), by = treatment_group][, .N, by = .(treatment_group, discharge_category)]
discharge_summary[, n_reportable := report_count(N)]
discharge_summary[, N := NULL]
discharge_summary[, small_cell_threshold := min_cell_n]
write_summary(discharge_summary, "discharge_outcomes_suppressed.csv")

quality_flags <- data.table(
  check = c(
    "one row per patient", "one row per hospital block",
    "five observed hospitals", "patient death_dttm available",
    "data cutoff explicitly supplied", "complete 28-day follow-up",
    "treatment grace period available", "complete grace-period follow-up",
    "negative treatment times", "death before time zero",
    "invalid length of stay", "missing baseline lactate",
    "indeterminate active catecholamine dose"
  ),
  value = c(
    !anyDuplicated(cohort$patient_id), !anyDuplicated(cohort$hospital_block_id),
    uniqueN(cohort$source_site_id) == 5L, death_source_available,
    !is.na(specified_cutoff_dttm), sum(cohort$complete_28d_followup),
    is.finite(grace_hours), sum(cohort$complete_grace_followup, na.rm = TRUE),
    if ("treatment_hours_from_time_zero" %in% names(cohort))
      sum(cohort$treatment_hours_from_time_zero < 0, na.rm = TRUE) else NA,
    sum(cohort$event_hours < 0, na.rm = TRUE),
    if ("length_of_stay_days" %in% names(cohort))
      sum(cohort$length_of_stay_days < 0, na.rm = TRUE) else NA,
    if ("lactate_missing_at_baseline" %in% names(cohort))
      sum(cohort$lactate_missing_at_baseline == 1L, na.rm = TRUE) else NA,
    if ("catecholamine_dose_indeterminate" %in% names(cohort))
      sum(cohort$catecholamine_dose_indeterminate == 1L, na.rm = TRUE) else NA
  ),
  note = c(
    "Expected TRUE", "Expected TRUE", "Expected TRUE for Hopkins extract",
    "If FALSE, 28-day mortality uses only Expired-discharge fallback",
    "Preferred; otherwise the latest observed timestamp is used",
    "Denominator for crude 28-day mortality",
    "Required for treated-within-grace assignment",
    "Denominator for treatment assignment and grace-period mortality",
    "Expected 0", "Expected 0", "Expected 0",
    "Retained and flagged; report as an eligibility-data diagnostic",
    "Retained and flagged; report as an eligibility-data diagnostic"
  )
)
write_summary(quality_flags, "data_quality_checks.csv")

# ---- Aggregate-only figures -------------------------------------------------

theme_set(theme_minimal(base_size = 13))
if (file.exists(attrition_path)) {
  attrition <- fread(attrition_path)
  attrition_plot <- attrition[!is.na(n_hospital_blocks)]
  attrition_plot[, step := factor(step, levels = rev(step))]
  p <- ggplot(attrition_plot, aes(n_hospital_blocks, step)) +
    geom_col(fill = "#2C7FB8") +
    geom_text(aes(label = format(n_hospital_blocks, big.mark = ",")),
              hjust = -0.08, size = 3.6) +
    expand_limits(x = max(attrition_plot$n_hospital_blocks) * 1.17) +
    labs(title = "Pneumonia cohort eligibility", x = "Hospital blocks", y = NULL)
  ggsave(file.path(output_dir, "plot_cohort_attrition.png"), p,
         width = 11, height = 6.8, dpi = 300)
}

site_plot <- copy(site_summary)[!is.na(n_reportable)]
if (nrow(site_plot)) {
  p <- ggplot(site_plot, aes(source_site_id, n_reportable)) +
    geom_col(fill = "#41AB5D") +
    geom_text(aes(label = format(n_reportable, big.mark = ",")),
              vjust = -0.35) +
    expand_limits(y = max(site_plot$n_reportable) * 1.12) +
    labs(title = "Final cohort by hospital", x = "Hospital", y = "Patients")
  ggsave(file.path(output_dir, "plot_cohort_by_site.png"), p,
         width = 8, height = 5.5, dpi = 300)
}

if ("length_of_stay_days" %in% names(cohort)) {
  los <- cohort[is.finite(length_of_stay_days) & length_of_stay_days >= 0,
                .(bin = cut(pmin(length_of_stay_days, 60),
                            breaks = seq(0, 60, by = 5), include.lowest = TRUE))][,
                  .N, by = bin]
  los[, n_reportable := report_count(N)]
  los <- los[!is.na(n_reportable)]
  if (nrow(los)) {
    p <- ggplot(los, aes(bin, n_reportable)) +
      geom_col(fill = "#756BB1") +
      labs(title = "Hospital length of stay",
           subtitle = "Values above 60 days truncated to the final bin",
           x = "Length of stay, days", y = "Patients") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(output_dir, "plot_length_of_stay.png"), p,
           width = 9, height = 5.8, dpi = 300)
  }
}

if ("treatment_hours_from_time_zero" %in% names(cohort) &&
    is.finite(grace_hours)) {
  timing <- cohort[
    treated_within_grace == 1L & is.finite(treatment_hours_from_time_zero),
    .(bin_start = floor(treatment_hours_from_time_zero / 4) * 4)
  ][, .N, by = bin_start]
  timing[, n_reportable := report_count(N)]
  timing <- timing[!is.na(n_reportable)]
  if (nrow(timing)) {
    p <- ggplot(timing, aes(bin_start, n_reportable)) +
      geom_col(width = 3.7, fill = "#E6550D") +
      labs(title = "Time to qualifying steroid treatment",
           subtitle = paste0("Four-hour bins within the ", grace_hours,
                             "-hour grace period"),
           x = "Hours from respiratory-failure time zero", y = "Patients")
    ggsave(file.path(output_dir, "plot_treatment_timing.png"), p,
           width = 8.5, height = 5.5, dpi = 300)
  }
}

mortality_plot <- copy(outcome_summary)[!is.na(death_28d_percent)]
if (nrow(mortality_plot)) {
  p <- ggplot(mortality_plot,
              aes(treatment_group, death_28d_percent, fill = treatment_group)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%", death_28d_percent)), vjust = -0.35) +
    expand_limits(y = max(mortality_plot$death_28d_percent) * 1.14) +
    labs(title = "Observed 28-day mortality by treatment assignment",
         x = NULL, y = "Mortality (%)") +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  ggsave(file.path(output_dir, "plot_mortality_by_treatment.png"), p,
         width = 8.5, height = 5.5, dpi = 300)
}

missing_plot <- missingness_summary[order(-missing_percent)][1:min(.N, 15L)]
if (nrow(missing_plot)) {
  missing_plot[, variable := factor(variable, levels = rev(variable))]
  p <- ggplot(missing_plot, aes(missing_percent, variable)) +
    geom_col(fill = "#636363") +
    labs(title = "Missingness in selected cohort variables",
         x = "Missing (%)", y = NULL)
  ggsave(file.path(output_dir, "plot_missingness.png"), p,
         width = 8.5, height = 6, dpi = 300)
}

cat("\nAggregate cohort EDA complete.\n")
cat("Final cohort patients: ", nrow(cohort), "\n", sep = "")
cat("Sites: ", paste(sort(unique(cohort$source_site_id)), collapse = ", "),
    "\n", sep = "")
if (treatment_group_available) {
  cat("Treated within ", grace_hours, " hours: ",
      sum(cohort$treated_within_grace[cohort$complete_grace_followup == 1L]),
      " of ", sum(cohort$complete_grace_followup),
      " with complete grace-period follow-up\n", sep = "")
}
cat("Observed deaths within 28 days: ",
    sum(cohort$death_within_28d[cohort$complete_28d_followup == 1L]),
    " of ", sum(cohort$complete_28d_followup),
    " with complete 28-day follow-up\n", sep = "")
cat("Patient-level inputs were not exported.\n")
cat("Summary results: ", normalizePath(output_dir), "\n", sep = "")
