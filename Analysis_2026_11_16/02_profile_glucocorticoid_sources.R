#!/usr/bin/env Rscript

# Targeted source check: intermittent med_group steroid.
# Confirmed source label: steroid; report its observed labels and counts.
# Dose/route eligibility is applied by 04, not by this source check.
# Usage:
#   Rscript 02_profile_glucocorticoid_sources.R CLIF_DATA_DIR [OUTPUT_DIR]

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
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_file[1L]), fixed = TRUE)))
} else normalizePath(getwd())
output_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "results", "default", "sensitive",
            "glucocorticoid_source_profile")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

path <- file.path(data_dir, "clif_medication_admin_intermittent.parquet")
if (!file.exists(path)) stop("Missing: ", path, call. = FALSE)
dataset <- open_dataset(path)
required <- c("hospitalization_id", "med_group", "admin_dttm", "med_category",
              "med_route_category", "med_dose", "med_dose_unit", "mar_action_category")
missing <- setdiff(required, dataset$schema$names)
if (length(missing)) stop("Intermittent table lacks: ", paste(missing, collapse = ", "))

# Confirmed source label; all counts describe administration rows.
steroid_rows <- dataset |> filter(.data$med_group == "steroid")
groups <- steroid_rows |> group_by(med_group) |>
  summarise(n_rows = n(), .groups = "drop") |> collect() |> as.data.table()
fwrite(groups,
       file.path(output_dir, "steroid_med_group_counts.csv"), na = "")
if (!nrow(groups)) stop("No med_group 'steroid' found; review source mapping.")

routes <- steroid_rows |> group_by(med_route_category) |>
  summarise(n_rows = n(), .groups = "drop") |> collect() |> as.data.table()
routes[, systemic_candidate := tolower(trimws(as.character(med_route_category))) %in% c("enteral", "iv")]
fwrite(routes, file.path(output_dir, "steroid_route_counts.csv"), na = "")
# Counts across the source table, not patients or qualifying treatment starts.
profile <- steroid_rows |> group_by(med_category, med_route_category, med_dose_unit) |>
  summarise(n_rows = n(), .groups = "drop") |> collect() |> as.data.table()
fwrite(profile, file.path(output_dir, "steroid_drug_route_unit_counts.csv"), na = "")
actions <- steroid_rows |> group_by(mar_action_category) |>
  summarise(n_rows = n(), .groups = "drop") |> collect() |> as.data.table()
fwrite(actions, file.path(output_dir, "steroid_administration_action_counts.csv"), na = "")
checks <- data.table(
  rule = c("source_table", "medication_group", "administration_time", "candidate_routes", "dose_threshold"),
  specified_value = c("clif_medication_admin_intermittent.parquet", "steroid",
                      "admin_dttm", "enteral;iv", "See 04: drug-specific dose and route rules"),
  available_or_confirmed = c(TRUE, TRUE, TRUE, any(routes$systemic_candidate), NA))
fwrite(checks, file.path(output_dir, "confirmed_steroid_source_check.csv"), na = "")
cat("Observed steroid med_group labels and administration-row counts:\n")
print(groups)
cat("\nSource check complete. Results: ", normalizePath(output_dir), "\n", sep = "")
cat("Only med_group steroid is accepted. Scripts 03/04 apply cohort, timing, drug, dose and route rules.\n")
