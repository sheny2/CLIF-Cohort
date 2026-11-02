#!/usr/bin/env Rscript

# Confirm and profile glucocorticoid administrations in a CLIF delivery. The
# Hopkins definition supplied by the clinical team is intermittent medication
# administration with med_group equal to "steroid" or "steroids" and an
# enteral or IV route. This
# remains a candidate exposure definition until the minimum dose and any dose-
# equivalence rules are supplied. Aggregate values and routes only are exported.
#
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
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else normalizePath(getwd())
output_dir <- if (length(args) >= 2L && nzchar(args[2L])) args[2L] else
  file.path(script_dir, "glucocorticoid_source_check")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

medication_files <- list.files(
  data_dir, pattern = "^clif_medication.*\\.parquet$", full.names = TRUE
)
if (!length(medication_files)) {
  stop("No clif_medication*.parquet tables were found in ", data_dir,
       call. = FALSE)
}

# Broad discovery vocabulary. The systemic/inhaled distinction is provisional:
# actual exposure status still requires route, MAR action, timing, and dose
# review with the clinical team.
dictionary <- data.table(
  search_term = c(
    "methylprednisolone", "prednisolone", "prednisone", "hydrocortisone",
    "dexamethasone", "betamethasone", "cortisone", "triamcinolone",
    "budesonide", "fluticasone", "beclomethasone", "mometasone",
    "glucocorticoid", "corticosteroid", "steroid"
  ),
  provisional_class = c(
    rep("usually systemic candidate", 7),
    "route-dependent candidate",
    rep("usually inhaled candidate", 4),
    rep("class-level term; manual review required", 3)
  )
)
fwrite(dictionary,
       file.path(output_dir, "glucocorticoid_search_dictionary.csv"), na = "")

text_candidates <- c(
  "med_group", "med_category", "med_name", "medication_name", "generic_name",
  "ingredient", "order_name"
)
time_candidates <- c(
  "admin_dttm", "recorded_dttm", "start_dttm", "order_dttm", "end_dttm"
)
route_candidates <- c("med_route_category", "med_route_name", "route")
action_candidates <- c("mar_action_category", "mar_action_name")
dose_candidates <- c("med_dose", "med_dose_unit", "med_amount", "med_amount_unit")

inventory_parts <- list()
value_parts <- list()
route_parts <- list()

classify_value <- function(value) {
  clean <- tolower(trimws(as.character(value)))
  hit <- vapply(dictionary$search_term, function(term) {
    grepl(term, clean, fixed = TRUE)
  }, logical(length(clean)))
  # vapply transposes when clean has length >1; normalize to rows=values.
  if (length(clean) == 1L) hit <- matrix(hit, nrow = 1L)
  if (nrow(hit) != length(clean)) hit <- t(hit)
  first_hit <- apply(hit, 1L, function(x) {
    index <- which(x)
    if (length(index)) index[1L] else NA_integer_
  })
  data.table(
    matched_term = dictionary$search_term[first_hit],
    provisional_class = dictionary$provisional_class[first_hit]
  )
}

for (path in medication_files) {
  table_name <- sub("\\.parquet$", "", basename(path))
  dataset <- open_dataset(path)
  available <- dataset$schema$names
  text_fields <- intersect(text_candidates, available)
  route_fields <- intersect(route_candidates, available)

  inventory_parts[[table_name]] <- data.table(
    source_table = table_name,
    column_count = length(available),
    has_hospitalization_id = "hospitalization_id" %in% available,
    medication_text_fields = paste(text_fields, collapse = ";"),
    administration_time_fields = paste(intersect(time_candidates, available),
                                       collapse = ";"),
    route_fields = paste(route_fields, collapse = ";"),
    mar_action_fields = paste(intersect(action_candidates, available),
                              collapse = ";"),
    dose_fields = paste(intersect(dose_candidates, available), collapse = ";")
  )

  for (field in text_fields) {
    message("Profiling ", table_name, "$", field, "...")
    counts <- tryCatch(
      dataset |>
        filter(!is.na(.data[[field]])) |>
        group_by(.data[[field]]) |>
        summarise(n_rows = n(), .groups = "drop") |>
        collect() |>
        as.data.table(),
      error = function(e) {
        warning("Could not aggregate ", table_name, "$", field,
                "; collecting distinct values without counts: ",
                conditionMessage(e))
        out <- dataset |>
          select(all_of(field)) |>
          distinct() |>
          collect() |>
          as.data.table()
        out[, n_rows := NA_integer_]
        out
      }
    )
    setnames(counts, field, "source_value")
    counts[, source_value := as.character(source_value)]
    classified <- classify_value(counts$source_value)
    counts <- cbind(counts, classified)
    counts <- counts[!is.na(matched_term)]
    if (!nrow(counts)) next
    counts[, `:=`(source_table = table_name, source_field = field)]
    value_parts[[paste(table_name, field, sep = "::")]] <- counts

    matched_values <- unique(counts$source_value)
    for (route_field in route_fields) {
      route_counts <- tryCatch(
        dataset |>
          filter(.data[[field]] %in% matched_values) |>
          group_by(.data[[field]], .data[[route_field]]) |>
          summarise(n_rows = n(), .groups = "drop") |>
          collect() |>
          as.data.table(),
        error = function(e) data.table()
      )
      if (!nrow(route_counts)) next
      setnames(route_counts, c(field, route_field),
               c("source_value", "route_value"))
      route_counts[, `:=`(
        source_table = table_name,
        source_field = field,
        route_field = route_field
      )]
      route_parts[[paste(table_name, field, route_field, sep = "::")]] <-
        route_counts
    }
  }
}

inventory <- rbindlist(inventory_parts, use.names = TRUE, fill = TRUE)
candidate_values <- if (length(value_parts)) {
  rbindlist(value_parts, use.names = TRUE, fill = TRUE)
} else data.table()
candidate_routes <- if (length(route_parts)) {
  rbindlist(route_parts, use.names = TRUE, fill = TRUE)
} else data.table()
if (!nrow(candidate_values)) {
  candidate_values <- data.table(
    source_table = character(), source_field = character(),
    source_value = character(), n_rows = integer(), matched_term = character(),
    provisional_class = character()
  )
}
if (!nrow(candidate_routes)) {
  candidate_routes <- data.table(
    source_table = character(), source_field = character(),
    source_value = character(), route_field = character(),
    route_value = character(), n_rows = integer()
  )
}
setcolorder(candidate_values, c(
  "source_table", "source_field", "source_value", "matched_term",
  "provisional_class", "n_rows"
))
setorder(candidate_values, source_table, source_field, provisional_class,
         -n_rows)
setorder(candidate_routes, source_table, source_field, source_value,
         route_field, -n_rows)

fwrite(inventory, file.path(output_dir, "medication_table_inventory.csv"), na = "")
fwrite(candidate_values,
       file.path(output_dir, "glucocorticoid_candidate_values.csv"), na = "")
fwrite(candidate_routes,
       file.path(output_dir, "glucocorticoid_candidate_routes.csv"), na = "")

confirmed_source_check <- data.table(
  rule = c(
    "source_table", "medication_group", "administration_time",
    "eligible_routes", "dose_threshold"
  ),
  specified_value = c(
    "clif_medication_admin_intermittent.parquet", "med_group in {steroid, steroids}",
    "admin_dttm", "med_route_category in {enteral, iv}",
    "PENDING CLINICAL SPECIFICATION"
  ),
  available_or_confirmed = c(
    "clif_medication_admin_intermittent" %in% inventory$source_table,
    any(candidate_values$source_table == "clif_medication_admin_intermittent" &
          candidate_values$source_field == "med_group" &
          tolower(trimws(candidate_values$source_value)) %in%
            c("steroid", "steroids")),
    any(inventory$source_table == "clif_medication_admin_intermittent" &
          grepl("(^|;)admin_dttm($|;)", inventory$administration_time_fields)),
    any(inventory$source_table == "clif_medication_admin_intermittent" &
          grepl("(^|;)med_route_category($|;)", inventory$route_fields)),
    FALSE
  )
)
fwrite(confirmed_source_check,
       file.path(output_dir, "confirmed_steroid_source_check.csv"), na = "")

cat("\nGlucocorticoid source profiling complete.\n")
cat("Medication tables inspected:", nrow(inventory), "\n")
cat("Candidate vocabulary values:", nrow(candidate_values), "\n")
cat("Results:", normalizePath(output_dir), "\n")
cat(paste(
  "Confirmed candidate rule: intermittent med_group in {'steroid','steroids'}, admin_dttm,",
  "and med_route_category in {'enteral','iv'}. The minimum qualifying dose,",
  "dose-equivalence rules, and treatment grace period remain to be finalized.\n"
))
