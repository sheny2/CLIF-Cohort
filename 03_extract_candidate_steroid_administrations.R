#!/usr/bin/env Rscript

# Extract candidate systemic glucocorticoid administrations for the final CLIF
# pneumonia cohort using the source fields confirmed by the Hopkins team:
#   clif_medication_admin_intermittent
#   med_group in c("steroid", "steroids")
#   med_route_category in c("enteral", "iv")
#   admin_dttm is the administration time
#
# This script intentionally does NOT classify treatment. Drug-specific dose and
# route rules are applied by 04_define_steroid_treatment.R.
#
# Usage:
#   Rscript 03_extract_candidate_steroid_administrations.R \
#     CLIF_DATA_DIR [COHORT_OUTPUT_DIR] [OUTPUT_DIR]

rm(list = ls())
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else
  Sys.getenv("CLIF_DATA_DIR", unset = "")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Supply a valid CLIF_DATA_DIR as argument 1 or an environment variable.",
       call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else normalizePath(getwd())
cohort_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "output")
output_dir <- if (length(args) >= 3L && nzchar(args[3L])) args[3L] else
  file.path(script_dir, "steroid_candidates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

medication_path <- file.path(
  data_dir, "clif_medication_admin_intermittent.parquet"
)
crosswalk_path <- file.path(
  cohort_dir, "adult_pneumonia_hospitalization_crosswalk.parquet"
)
cohort_path <- file.path(cohort_dir, "adult_pneumonia_cohort.parquet")
required_paths <- c(medication_path, crosswalk_path, cohort_path)
if (any(!file.exists(required_paths))) {
  stop("Required input file(s) missing:\n  ",
       paste(required_paths[!file.exists(required_paths)], collapse = "\n  "),
       call. = FALSE)
}

dataset <- open_dataset(medication_path)
available <- dataset$schema$names
required_columns <- c(
  "hospitalization_id", "med_group", "admin_dttm", "med_route_category"
)
missing_columns <- setdiff(required_columns, available)
if (length(missing_columns)) {
  stop("Intermittent medication table is missing: ",
       paste(missing_columns, collapse = ", "), call. = FALSE)
}

# Discover exact stored capitalization before applying Arrow filters.
group_values <- dataset |>
  select(med_group) |>
  distinct() |>
  collect() |>
  as.data.table()
steroid_group_values <- as.character(group_values$med_group)[
  tolower(trimws(as.character(group_values$med_group))) %in%
    c("steroid", "steroids")
]
route_values <- dataset |>
  select(med_route_category) |>
  distinct() |>
  collect() |>
  as.data.table()
eligible_route_values <- as.character(route_values$med_route_category)[
  tolower(trimws(as.character(route_values$med_route_category))) %in%
    c("enteral", "iv")
]
if (!length(steroid_group_values)) {
  stop("No med_group value equal to 'steroid' or 'steroids' was found.",
       call. = FALSE)
}
if (!length(eligible_route_values)) {
  stop("No enteral or iv med_route_category value was found.", call. = FALSE)
}

optional_columns <- c(
  "med_category", "med_name", "medication_name", "med_dose",
  "med_dose_unit", "med_amount", "med_amount_unit", "mar_action_category",
  "mar_action_name", "order_dttm", "recorded_dttm", "start_dttm", "end_dttm"
)
selected_columns <- unique(c(required_columns,
                             intersect(optional_columns, available)))

crosswalk <- as.data.table(read_parquet(crosswalk_path))
cohort <- as.data.table(read_parquet(cohort_path))
if (!all(c("hospitalization_id", "hospital_block_id", "patient_id") %in%
         names(crosswalk))) {
  stop("Selected cohort crosswalk lacks hospitalization/block/patient IDs.",
       call. = FALSE)
}
if (!"hospital_block_id" %in% names(cohort)) {
  stop("Final cohort lacks hospital_block_id.", call. = FALSE)
}
crosswalk[, hospitalization_id := as.character(hospitalization_id)]
cohort_hospitalization_ids <- unique(crosswalk$hospitalization_id[
  !is.na(crosswalk$hospitalization_id) & nzchar(crosswalk$hospitalization_id)
])

# Keep the large medication table in Arrow. Only rows belonging to the final
# cohort are collected into memory; source-wide attrition counts are aggregated
# by Arrow first.
steroid_query <- dataset |>
  filter(.data$med_group %in% steroid_group_values)
route_query <- steroid_query |>
  filter(.data$med_route_category %in% eligible_route_values)
valid_time_query <- route_query |>
  filter(!is.na(.data$admin_dttm))

query_counts <- function(query) {
  n_rows <- query |>
    summarise(n_rows = n()) |>
    collect() |>
    as.data.table()
  n_hospitalizations <- query |>
    select(hospitalization_id) |>
    filter(!is.na(.data$hospitalization_id)) |>
    distinct() |>
    summarise(n_hospitalizations = n()) |>
    collect() |>
    as.data.table()
  c(
    n_rows = as.integer(n_rows$n_rows[1L]),
    n_hospitalizations = as.integer(n_hospitalizations$n_hospitalizations[1L])
  )
}
all_counts <- query_counts(steroid_query)
route_counts <- query_counts(route_query)
valid_time_counts <- query_counts(valid_time_query)

candidate <- valid_time_query |>
  filter(.data$hospitalization_id %in% cohort_hospitalization_ids) |>
  select(all_of(selected_columns)) |>
  collect() |>
  as.data.table()
candidate[, `:=`(
  hospitalization_id = as.character(hospitalization_id),
  admin_dttm = as.POSIXct(admin_dttm, tz = "UTC")
)]

crosswalk_columns <- intersect(
  c("hospitalization_id", "hospital_block_id", "patient_id", "source_site_id",
    "component_site_id"),
  names(crosswalk)
)
cohort_columns <- intersect(
  c("hospital_block_id", "respiratory_failure_start", "admission_dttm",
    "discharge_dttm"),
  names(cohort)
)
candidate <- merge(
  candidate,
  unique(crosswalk[, ..crosswalk_columns]),
  by = "hospitalization_id", all = FALSE, sort = FALSE
)
candidate <- merge(
  candidate,
  unique(cohort[, ..cohort_columns]),
  by = "hospital_block_id", all.x = TRUE, sort = FALSE
)
if ("respiratory_failure_start" %in% names(candidate)) {
  candidate[, respiratory_failure_start := as.POSIXct(
    respiratory_failure_start, tz = "UTC"
  )]
  candidate[, hours_from_time_zero := as.numeric(
    difftime(admin_dttm, respiratory_failure_start, units = "hours")
  )]
}
setorder(candidate, patient_id, hospital_block_id, admin_dttm)

audit <- data.table(
  step = c(
    "All intermittent rows with med_group in {steroid, steroids}",
    "Steroid rows with enteral or iv route",
    "Route-eligible steroid rows with nonmissing admin_dttm",
    "Candidate administrations linked to final cohort"
  ),
  n_administration_rows = c(
    all_counts[["n_rows"]], route_counts[["n_rows"]],
    valid_time_counts[["n_rows"]], nrow(candidate)
  ),
  n_hospitalizations = c(
    all_counts[["n_hospitalizations"]],
    route_counts[["n_hospitalizations"]],
    valid_time_counts[["n_hospitalizations"]],
    uniqueN(candidate$hospitalization_id, na.rm = TRUE)
  ),
  n_patients = c(NA_integer_, NA_integer_, NA_integer_,
                 uniqueN(candidate$patient_id, na.rm = TRUE))
)

definition <- data.table(
  parameter = c(
    "source_table", "med_group", "eligible_routes", "administration_time",
    "dose_route_definition", "treatment_classification"
  ),
  value = c(
    "clif_medication_admin_intermittent", "steroid;steroids", "enteral;iv",
    "admin_dttm", "See 04_define_steroid_treatment.R", "NOT YET PERFORMED"
  )
)

write_parquet(as.data.frame(candidate), file.path(
  output_dir, "candidate_steroid_administrations_in_final_cohort.parquet"
))
fwrite(audit, file.path(output_dir, "candidate_steroid_attrition.csv"), na = "")
fwrite(definition, file.path(output_dir, "candidate_steroid_definition.csv"),
       na = "")

cat("\nCandidate steroid extraction complete.\n")
print(audit)
cat("\nNo treatment indicator was assigned here; run ",
    "04_define_steroid_treatment.R next.\n", sep = "")
cat("Results: ", normalizePath(output_dir), "\n", sep = "")
