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
    "diagnosis_primary", "poa_present"
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
  ),
  clif_medication_admin_continuous.parquet = c(
    "hospitalization_id", "admin_dttm", "med_category", "med_dose",
    "med_dose_unit", "mar_action_category"
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

# Multisite provenance is not part of the minimum cohort schema because some
# CLIF deliveries are split into one directory per site. In a combined Hopkins
# extract, prefer the admission hospital from the first ADT record. A site
# column carried directly by hospitalization is also supported.
hospitalization_path <- file.path(data_dir, "clif_hospitalization.parquet")
hospitalization_columns <- if (file.exists(hospitalization_path)) {
  open_dataset(hospitalization_path)$schema$names
} else character()
adt_path <- file.path(data_dir, "clif_adt.parquet")
adt_columns <- if (file.exists(adt_path)) open_dataset(adt_path)$schema$names else
  character()
site_column_candidates <- c(
  "admission_hospital_id", "admission_hospital_name", "hospital_id",
  "hospital_name", "site_id"
)
site_check <- data.table(
  source = c("clif_adt.parquet", "clif_hospitalization.parquet", "fallback"),
  available = c(
    file.exists(adt_path) && all(c("hospitalization_id", "hospital_id", "in_dttm") %in%
                                  adt_columns),
    any(site_column_candidates %in% hospitalization_columns),
    nzchar(Sys.getenv("CLIF_SITE_ID", unset = ""))
  ),
  selected_or_available_field = c(
    if ("hospital_id" %in% adt_columns) "hospital_id (first ADT location)" else "",
    paste(intersect(site_column_candidates, hospitalization_columns), collapse = ";"),
    Sys.getenv("CLIF_SITE_ID", unset = "")
  ),
  note = c(
    "Preferred for a combined multisite Hopkins extract",
    "Used when CLIF_SITE_COLUMN names one of these fields",
    "Use only for a directory already separated to one site"
  )
)
fwrite(site_check, file.path(output_dir, "site_provenance_check.csv"), na = "")

optional_medication_files <- c(
  "clif_medication_admin_intermittent.parquet",
  "clif_medication_admin_continuous.parquet"
)
medication_check <- rbindlist(lapply(optional_medication_files, function(filename) {
  path <- file.path(data_dir, filename)
  available <- if (file.exists(path)) open_dataset(path)$schema$names else character()
  data.table(
    file = filename,
    exists = file.exists(path),
    has_hospitalization_id = "hospitalization_id" %in% available,
    medication_text_fields = paste(
      intersect(c("med_group", "med_category", "med_name", "medication_name"),
                available),
      collapse = ";"
    ),
    has_confirmed_steroid_group_field = "med_group" %in% available,
    has_confirmed_administration_time = "admin_dttm" %in% available,
    has_confirmed_route_field = "med_route_category" %in% available,
    administration_time_fields = paste(
      intersect(c("admin_dttm", "recorded_dttm", "start_dttm"), available),
      collapse = ";"
    ),
    route_fields = paste(
      intersect(c("med_route_category", "med_route_name"), available),
      collapse = ";"
    )
  )
}), fill = TRUE)
fwrite(medication_check,
       file.path(output_dir, "optional_medication_table_check.csv"), na = "")

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
  poa_present = list(
    file = "clif_hospital_diagnosis.parquet", column = "poa_present"
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
  ),
  continuous_med_category = list(
    file = "clif_medication_admin_continuous.parquet", column = "med_category"
  ),
  continuous_med_dose_unit = list(
    file = "clif_medication_admin_continuous.parquet", column = "med_dose_unit"
  ),
  continuous_mar_action_category = list(
    file = "clif_medication_admin_continuous.parquet",
    column = "mar_action_category"
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
    "pneumonia present-on-admission field is populated",
    "qualifying respiratory device category observed",
    "SpO2 vital category observed",
    "lactate lab category observed",
    "shock vasopressor medication category observed"
  ),
  passed = c(
    category_rows[field == "admission_type_category",
                  any(tolower(trimws(value)) == "ed")],
    category_rows[field == "diagnosis_code_format",
                  any(grepl("ICD.?9|ICD.?10|^9$|^10$", value,
                            ignore.case = TRUE))],
    category_rows[field == "poa_present", any(!is.na(value) & nzchar(value))],
    category_rows[field == "device_category",
                  any(tolower(trimws(value)) %in%
                        c("imv", "nippv", "high_flow_nc", "nasal_cannula"))],
    category_rows[field == "vital_category",
                  any(tolower(trimws(value)) == "spo2")],
    category_rows[field == "lab_category",
                  any(tolower(trimws(value)) == "lactate")],
    category_rows[field == "continuous_med_category",
                  any(tolower(trimws(value)) %in%
                        c("norepinephrine", "epinephrine", "vasopressin"))]
  ),
  note = c(
    "Required by the current ED-only cohort definition",
    "Required by the pneumonia ICD classifier",
    "Required by the POA pneumonia definition",
    "Required by the initiation-based respiratory-support rule",
    "Required for the HFNO S/F branch",
    "Required by the concurrent-shock exclusion",
    "Required by the concurrent-shock exclusion"
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
