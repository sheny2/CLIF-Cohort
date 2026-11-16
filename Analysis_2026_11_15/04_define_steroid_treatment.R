#!/usr/bin/env Rscript

# Apply the clinician-specified systemic glucocorticoid dose/route definition
# to candidate administrations extracted by 03. Patient-level parquet outputs
# must remain on the approved Discovery HPC. CSV outputs use small-cell
# suppression and are intended for review as potentially portable summaries.
#
# Usage:
#   Rscript 04_define_steroid_treatment.R \
#     [CANDIDATE_DIR] [COHORT_OUTPUT_DIR] [OUTPUT_DIR] [GRACE_HOURS]
#
# GRACE_HOURS defaults to 24 hours and can also be set with
# TREATMENT_GRACE_HOURS. All qualifying administrations remain available for
# the independent 24/96-hour EDA, regardless of this assignment window.

rm(list = ls())
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_file[1L]),
                             fixed = TRUE)))
} else normalizePath(getwd())

candidate_dir <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else
  file.path(script_dir, "pipeline_output", "sensitive", "steroid_candidates")
cohort_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "pipeline_output", "sensitive", "cohort")
output_dir <- if (length(args) >= 3L && nzchar(args[3L])) args[3L] else
  file.path(script_dir, "pipeline_output", "sensitive", "steroid_treatment")
grace_text <- if (length(args) >= 4L && nzchar(args[4L])) args[4L] else
  Sys.getenv("TREATMENT_GRACE_HOURS", unset = "24")
grace_hours <- if (nzchar(grace_text)) suppressWarnings(as.numeric(grace_text)) else
  NA_real_
if (nzchar(grace_text) && (!is.finite(grace_hours) || grace_hours < 0)) {
  stop("GRACE_HOURS/TREATMENT_GRACE_HOURS must be nonnegative.", call. = FALSE)
}
small_cell_threshold <- suppressWarnings(as.integer(
  Sys.getenv("SUMMARY_SMALL_CELL_THRESHOLD", unset = "10")
))
if (!is.finite(small_cell_threshold) || small_cell_threshold < 1L) {
  stop("SUMMARY_SMALL_CELL_THRESHOLD must be a positive integer.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

candidate_path <- file.path(
  candidate_dir, "candidate_steroid_administrations_in_final_cohort.parquet"
)
cohort_path <- file.path(cohort_dir, "adult_pneumonia_cohort.parquet")
if (!file.exists(candidate_path)) stop("Missing: ", candidate_path, call. = FALSE)
if (!file.exists(cohort_path)) stop("Missing: ", cohort_path, call. = FALSE)

administrations <- as.data.table(read_parquet(candidate_path))
cohort <- as.data.table(read_parquet(cohort_path))
required_administration_columns <- c(
  "patient_id", "hospital_block_id", "admin_dttm", "med_group",
  "med_category", "med_dose", "med_dose_unit", "med_route_category"
)
missing_administration_columns <- setdiff(
  required_administration_columns, names(administrations)
)
if (length(missing_administration_columns)) {
  stop("Candidate administrations are missing: ",
       paste(missing_administration_columns, collapse = ", "), call. = FALSE)
}
if (!all(c("patient_id", "hospital_block_id") %in% names(cohort))) {
  stop("Final cohort must contain patient_id and hospital_block_id.",
       call. = FALSE)
}

normalize_text <- function(x) tolower(trimws(as.character(x)))
normalize_drug <- function(x) gsub("[^a-z]", "", normalize_text(x))
normalize_unit <- function(x) {
  value <- normalize_text(x)
  value <- gsub("μ", "u", value, fixed = TRUE)
  gsub("[^a-z]", "", value)
}

# Thresholds are per administration, as specified by the clinical team. The
# hydrocortisone 50-mg threshold reflects expected q6h/q8h dosing and the goal
# of identifying regimens exceeding approximately 150 mg/day hydrocortisone.
dose_rules <- data.table(
  med_category_clean = c(
    "dexamethasone", "hydrocortisone", "prednisone",
    "methylprednisolone", "prednisolone"
  ),
  minimum_dose_mg = c(6, 50, 40, 32, 40),
  eligible_routes = c("enteral;iv", "iv", "enteral", "enteral;iv", "enteral;iv"),
  rationale = c(
    "Clinician-specified minimum dose",
    "50 mg per administration; typically q6h/q8h",
    "Clinician-specified minimum dose",
    "Clinician-specified minimum dose",
    "Clinician-specified minimum dose"
  )
)
fwrite(dose_rules, file.path(output_dir, "steroid_dose_route_rules.csv"), na = "")

administrations[, `:=`(
  med_group_clean = normalize_text(med_group),
  med_category_clean = normalize_drug(med_category),
  med_route_clean = normalize_text(med_route_category),
  med_dose_numeric = suppressWarnings(as.numeric(med_dose)),
  med_dose_unit_clean = normalize_unit(med_dose_unit),
  admin_dttm = as.POSIXct(admin_dttm, tz = "UTC")
)]
administrations[, dose_mg := fcase(
  med_dose_unit_clean %in% c("mg", "milligram", "milligrams"), med_dose_numeric,
  med_dose_unit_clean %in% c("mcg", "ug", "microgram", "micrograms"),
    med_dose_numeric / 1000,
  med_dose_unit_clean %in% c("g", "gram", "grams"), med_dose_numeric * 1000,
  default = NA_real_
)]
administrations[dose_rules, on = "med_category_clean", `:=`(
  minimum_dose_mg = i.minimum_dose_mg,
  eligible_routes = i.eligible_routes
)]
administrations[, route_eligible_for_drug := mapply(
  function(route, allowed) {
    if (is.na(route) || is.na(allowed)) return(FALSE)
    route %in% strsplit(allowed, ";", fixed = TRUE)[[1L]]
  }, med_route_clean, eligible_routes
)]
administrations[, qualification_reason := fcase(
  !med_group_clean %in% c("steroid", "steroids"), "not_steroid_med_group",
  is.na(minimum_dose_mg), "med_category_not_in_definition",
  is.na(med_dose_numeric) | med_dose_numeric <= 0, "missing_or_nonpositive_dose",
  is.na(dose_mg), "unsupported_or_missing_dose_unit",
  !route_eligible_for_drug, "route_not_eligible_for_drug",
  dose_mg < minimum_dose_mg, "dose_below_drug_threshold",
  default = "qualifying_treatment_administration"
)]
administrations[, qualifies_as_treatment := as.integer(
  qualification_reason == "qualifying_treatment_administration"
)]

if (!"hours_from_time_zero" %in% names(administrations) &&
    "respiratory_failure_start" %in% names(administrations)) {
  administrations[, respiratory_failure_start := as.POSIXct(
    respiratory_failure_start, tz = "UTC"
  )]
  administrations[, hours_from_time_zero := as.numeric(difftime(
    admin_dttm, respiratory_failure_start, units = "hours"
  ))]
}

qualifying <- administrations[qualifies_as_treatment == 1L]
setorder(qualifying, patient_id, hospital_block_id, admin_dttm)
first_qualifying <- qualifying[, .SD[1L], by = .(patient_id, hospital_block_id)]

write_parquet(
  as.data.frame(administrations),
  file.path(output_dir, "steroid_administrations_with_treatment_eligibility.parquet")
)
write_parquet(
  as.data.frame(qualifying),
  file.path(output_dir, "qualifying_steroid_administrations.parquet")
)
write_parquet(
  as.data.frame(first_qualifying),
  file.path(output_dir, "first_qualifying_steroid_by_block.parquet")
)

# Aggregate diagnostic tables use suppression so they can be reviewed for
# export without exposing rare drug/route combinations.
qualification_summary <- administrations[, .(
  n_administration_rows = .N,
  n_hospital_blocks = uniqueN(hospital_block_id),
  n_patients = uniqueN(patient_id)
), by = .(med_category_clean, med_route_clean, qualification_reason)]
for (column in c("n_administration_rows", "n_hospital_blocks", "n_patients")) {
  qualification_summary[, paste0(column, "_reportable") := fifelse(
    get(column) < small_cell_threshold, NA_integer_, as.integer(get(column))
  )]
  qualification_summary[, (column) := NULL]
}
qualification_summary[, small_cell_threshold := small_cell_threshold]
fwrite(
  qualification_summary,
  file.path(output_dir, "steroid_qualification_summary_suppressed.csv"), na = ""
)

cohort_columns <- intersect(
  c("patient_id", "hospital_block_id", "source_site_id",
    "respiratory_failure_start", "admission_dttm", "discharge_dttm"),
  names(cohort)
)
treatment_assignment <- unique(cohort[, ..cohort_columns])
treatment_assignment[, ever_qualifying_steroid := as.integer(
  hospital_block_id %in% qualifying$hospital_block_id
)]

if (is.finite(grace_hours)) {
  if (!"hours_from_time_zero" %in% names(qualifying)) {
    stop("A grace period was requested, but hours_from_time_zero cannot be derived.",
         call. = FALSE)
  }
  within_grace <- qualifying[
    !is.na(hours_from_time_zero) & hours_from_time_zero >= 0 &
      hours_from_time_zero <= grace_hours
  ]
  setorder(within_grace, patient_id, hospital_block_id, admin_dttm)
  first_within_grace <- within_grace[, .SD[1L],
                                    by = .(patient_id, hospital_block_id)]
  treatment_fields <- intersect(
    c("patient_id", "hospital_block_id", "admin_dttm", "hours_from_time_zero",
      "med_category_clean", "dose_mg", "med_route_clean"),
    names(first_within_grace)
  )
  first_within_grace <- first_within_grace[, ..treatment_fields]
  setnames(
    first_within_grace,
    setdiff(names(first_within_grace), c("patient_id", "hospital_block_id")),
    paste0("treatment_", setdiff(names(first_within_grace),
                                 c("patient_id", "hospital_block_id")))
  )
  treatment_assignment <- merge(
    treatment_assignment, first_within_grace,
    by = c("patient_id", "hospital_block_id"), all.x = TRUE, sort = FALSE
  )
  treatment_assignment[, treated_within_grace := as.integer(
    !is.na(treatment_admin_dttm)
  )]
  treatment_assignment[, treatment_grace_hours := grace_hours]
}
write_parquet(
  as.data.frame(treatment_assignment),
  file.path(output_dir, "steroid_treatment_assignment.parquet")
)

definition <- rbindlist(list(
  data.table(parameter = c(
    "medication_table", "med_group", "administration_time", "dose_rule",
    "summary_small_cell_threshold", "grace_hours"
  ), value = c(
    "clif_medication_admin_intermittent", "steroid;steroids", "admin_dttm",
    "drug-specific per-administration threshold; see steroid_dose_route_rules.csv",
    as.character(small_cell_threshold),
    if (is.finite(grace_hours)) as.character(grace_hours) else "NOT SPECIFIED"
  ))
), use.names = TRUE)
fwrite(definition, file.path(output_dir, "steroid_treatment_definition.csv"),
       na = "")

cat("\nSteroid treatment definition applied.\n")
cat("Candidate administration rows: ", nrow(administrations), "\n", sep = "")
cat("Qualifying administration rows: ", nrow(qualifying), "\n", sep = "")
cat("Blocks with qualifying steroid: ",
    uniqueN(qualifying$hospital_block_id), "\n", sep = "")
if (is.finite(grace_hours)) {
  cat("Blocks treated within ", grace_hours, " hours: ",
      sum(treatment_assignment$treated_within_grace), "\n", sep = "")
} else {
  cat("No grace period supplied; treated_within_grace was not assigned.\n")
}
cat("Patient-level parquet files must remain on the Discovery HPC.\n")
cat("Results: ", normalizePath(output_dir), "\n", sep = "")
