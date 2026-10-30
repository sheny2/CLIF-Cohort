#!/usr/bin/env Rscript

# Validate that a Hopkins CLIF extract contains the files and columns required
# by 01_build_adult_pneumonia_cohort.R. This script reads schemas and small
# category samples only; it does not construct or export a patient cohort.
#
# Usage:
#   Rscript 00_check_clif_cohort_inputs.R CLIF_DATA_DIR [OUTPUT_DIR]

rm(list = ls())
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else
  Sys.getenv("CLIF_DATA_DIR", unset = "")
if (!nzchar(data_dir)) {
  stop(
    "Supply CLIF_DATA_DIR as the first argument or set the environment variable.",
    call. = FALSE
  )
}
if (!dir.exists(data_dir)) stop("CLIF data directory not found: ", data_dir)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else normalizePath(getwd())
output_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "input_check")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

requirements <- list(
  clif_hospitalization.parquet = c(
    "hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm",
    "age_at_admission", "admission_type_category"
  ),
  clif_hospital_diagnosis.parquet = c(
    "hospitalization_id", "diagnosis_code", "diagnosis_code_format",
    "diagnosis_primary"
  ),
  clif_respiratory_support.parquet = c(
    "hospitalization_id", "recorded_dttm", "device_category", "fio2_set",
    "lpm_set", "tracheostomy"
  ),
  clif_vitals.parquet = c(
    "hospitalization_id", "recorded_dttm", "vital_category", "vital_value"
  ),
  clif_labs.parquet = c(
    "hospitalization_id", "lab_collect_dttm", "lab_category",
    "lab_value_numeric"
  ),
  clif_code_status.parquet = c(
    "patient_id", "start_dttm", "code_status_category"
  )
)

schema_rows <- rbindlist(lapply(names(requirements), function(filename) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) {
    return(data.table(
      file = filename, exists = FALSE, required_columns_present = FALSE,
      missing_columns = paste(requirements[[filename]], collapse = ";"),
      available_column_count = NA_integer_
    ))
  }
  available <- open_dataset(path)$schema$names
  missing <- setdiff(requirements[[filename]], available)
  data.table(
    file = filename,
    exists = TRUE,
    required_columns_present = !length(missing),
    missing_columns = paste(missing, collapse = ";"),
    available_column_count = length(available)
  )
}))
fwrite(schema_rows, file.path(output_dir, "required_schema_check.csv"), na = "")

sample_column <- function(filename, columns, n = 100000L) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) return(data.table())
  available <- open_dataset(path)$schema$names
  columns <- intersect(columns, available)
  if (!length(columns)) return(data.table())
  open_dataset(path) |>
    select(all_of(columns)) |>
    head(n) |>
    collect() |>
    as.data.table()
}

category_specs <- list(
  admission_type_category = list(
    file = "clif_hospitalization.parquet", column = "admission_type_category"
  ),
  diagnosis_code_format = list(
    file = "clif_hospital_diagnosis.parquet", column = "diagnosis_code_format"
  ),
  device_category = list(
    file = "clif_respiratory_support.parquet", column = "device_category"
  ),
  vital_category = list(
    file = "clif_vitals.parquet", column = "vital_category"
  ),
  lab_category = list(
    file = "clif_labs.parquet", column = "lab_category"
  ),
  code_status_category = list(
    file = "clif_code_status.parquet", column = "code_status_category"
  )
)

category_rows <- rbindlist(lapply(names(category_specs), function(label) {
  spec <- category_specs[[label]]
  x <- sample_column(spec$file, spec$column)
  if (!nrow(x) || !spec$column %in% names(x)) return(data.table())
  values <- as.character(x[[spec$column]])
  data.table(value = values)[, .N, by = value][order(-N)][, `:=`(
    field = label,
    sampled_rows = length(values)
  )]
}), fill = TRUE)
if (!nrow(category_rows)) {
  category_rows <- data.table(
    field = character(), value = character(), N = integer(),
    sampled_rows = integer()
  )
} else {
  setcolorder(category_rows, c("field", "value", "N", "sampled_rows"))
}
fwrite(category_rows, file.path(output_dir, "sampled_category_values.csv"),
       na = "")

required_category_checks <- data.table(
  check = c(
    "ED admission category observed",
    "ICD-9-CM or ICD-10-CM format observed",
    "qualifying respiratory device category observed",
    "SpO2 vital category observed",
    "arterial PaO2 lab category observed"
  ),
  passed = c(
    category_rows[field == "admission_type_category",
                  any(tolower(trimws(value)) == "ed")],
    category_rows[field == "diagnosis_code_format",
                  any(grepl("ICD.?9|ICD.?10|^9$|^10$", value,
                            ignore.case = TRUE))],
    category_rows[field == "device_category",
                  any(tolower(trimws(value)) %in%
                        c("imv", "nippv", "high_flow_nc", "nasal_cannula"))],
    category_rows[field == "vital_category",
                  any(tolower(trimws(value)) == "spo2")],
    category_rows[field == "lab_category",
                  any(tolower(trimws(value)) == "po2_arterial")]
  ),
  note = c(
    "Required by the current ED-only cohort definition",
    "Required by the pneumonia ICD classifier",
    "Required by the respiratory-support eligibility rule",
    "Required for the S/F branch",
    "Required for the P/F branch; the cohort can still qualify through S/F"
  )
)
fwrite(required_category_checks,
       file.path(output_dir, "category_compatibility_check.csv"), na = "")

cat("\nRequired schema check:\n")
print(schema_rows)
cat("\nCategory compatibility check (based on at most 100,000 rows/table):\n")
print(required_category_checks)

if (!all(schema_rows$exists & schema_rows$required_columns_present)) {
  stop(
    "One or more required CLIF files/columns are unavailable. Review ",
    file.path(output_dir, "required_schema_check.csv"),
    call. = FALSE
  )
}

cat("\nInput check completed. Results saved to:\n",
    normalizePath(output_dir), "\n", sep = "")
