#!/usr/bin/env Rscript
# Separate paired sensitivity analysis; never edits the parent scripts or inputs.
suppressPackageStartupMessages({library(arrow); library(data.table); library(dplyr)})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 5L) stop(paste(
  "Usage: Rscript historical_proxy_sensitivity/run_analysis.R",
  "CLIF_DATA_DIR COHORT_DIR TREATMENT_DIR HISTORY_CLIF_DIR NEW_OUTPUT_ROOT"))
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
folder <- dirname(self)
parent_dir <- dirname(folder)
source(file.path(folder, "history_helpers.R"))
history_mode <- Sys.getenv("HISTORY_MODE", "observed_extract")
if (!history_mode %in% c("observed_extract", "complete_history")) stop("HISTORY_MODE must be observed_extract or complete_history.")
# The project documents a January 2021 extraction start. Override if different.
coverage_text <- Sys.getenv("HISTORY_COVERAGE_START", if (history_mode == "observed_extract") "2021-01-01" else "")
coverage_start <- as.POSIXct(paste(coverage_text, "00:00:00"), tz = "UTC")
if (!nzchar(coverage_text) || is.na(coverage_start)) stop("Set HISTORY_COVERAGE_START=YYYY-MM-DD to the extract start.")
if (history_mode == "complete_history" && Sys.getenv("HISTORY_SOURCE_COMPLETE", "") != "yes") stop(
  "Complete-history mode requires HISTORY_SOURCE_COMPLETE=yes and unfiltered historical tables.")
if (history_mode == "observed_extract") message(
  "Using observed history in the supplied extract. Zero means no recorded qualifying prior event; it does not establish absence of prior hospitalization.")
grace <- suppressWarnings(as.numeric(Sys.getenv("TREATMENT_GRACE_HOURS", "96")))
if (length(grace) != 1L || !is.finite(grace) || grace < 6 || grace >= 672 || grace %% 6 != 0) stop("Grace must be a multiple of six hours in [6, 672).")
# A new root prevents accidental overwrite of the original analysis or a prior run.
output <- argv[5]
if (dir.exists(output) || file.exists(output)) stop("NEW_OUTPUT_ROOT already exists; choose a fresh directory.")
cohort <- as.data.table(read_parquet(file.path(argv[2], "adult_pneumonia_cohort.parquet")))
history_dir <- argv[4]
patient_ids <- as.character(cohort$patient_id)
read_history <- function(file, columns, ids = NULL) {
  ds <- open_dataset(file.path(history_dir, file))
  if (!all(columns %in% ds$schema$names)) stop("Missing history columns in ", file)
  q <- ds |> select(all_of(columns))
  if (!is.null(ids)) q <- q |> filter(.data$hospitalization_id %in% ids)
  else q <- q |> filter(.data$patient_id %in% patient_ids)
  as.data.table(collect(q))
}
h <- read_history("clif_hospitalization.parquet", c("hospitalization_id", "patient_id", "admission_dttm", "discharge_dttm", "admission_type_category"))
ids <- as.character(h$hospitalization_id)
a <- read_history("clif_adt.parquet", c("hospitalization_id", "in_dttm", "location_category"), ids)
r <- read_history("clif_respiratory_support.parquet", c("hospitalization_id", "recorded_dttm", "device_category"), ids)
d <- read_history("clif_hospital_diagnosis.parquet", c("hospitalization_id", "diagnosis_code", "diagnosis_code_format"), ids)
proxies <- build_history_proxies(cohort, h, a, r, d, coverage_start, history_mode = history_mode)
if (history_mode == "complete_history" && !any(proxies$full_calendar_lookback == 1L)) stop("No patients have a full 365-day historical coverage window.")
dir.create(file.path(output, "sensitive"), recursive = TRUE)
write_parquet(proxies, file.path(output, "sensitive", "historical_proxies.parquet"))
fwrite(data.table(file = c("07_prepare_and_run_five_methods.R", "06_generalized_five_methods.R"),
                  md5 = unname(tools::md5sum(file.path(parent_dir, c("07_prepare_and_run_five_methods.R", "06_generalized_five_methods.R"))))),
       file.path(output, "sensitive", "parent_source_versions.csv"))
# Fail closed if the parent's insertion points change. Transform parsed R calls,
# not text, and evaluate isolated copies entirely in memory.
program <- parse(file.path(parent_dir, "07_prepare_and_run_five_methods.R"))
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
base_index <- which(vapply(program, function(x) is_assignment(x, "baseline_covariates"), logical(1)))[1]
stopifnot(is.finite(base_index), sum(vapply(program, function(x) is_assignment(x, "grace_hours"), logical(1))) == 1L)
for (mode in c("reference", "augmented")) {
  e <- new.env(parent = globalenv())
  e$script_dir <- parent_dir
  parent_args <- c(argv[1:3], file.path(output, "sensitive", mode),
                   file.path(output, "shareable", mode), Sys.getenv("CCW_MAX_PATIENTS", "0"),
                   Sys.getenv("CLIF_DATA_CUTOFF_DATE", "2026-10-31"))
  for (j in seq_along(program)) {
    expr <- program[[j]]
    # Parent cleanup must not remove wrapper state or injected configuration.
    if (is.call(expr) && identical(expr[[1]], as.name("rm"))) next
    e$script_dir <- parent_dir
    if (is_assignment(expr, "args")) { e$args <- parent_args; next }
    if (is_assignment(expr, "grace_hours")) { e$grace_hours <- as.integer(grace); next }
    if (j == base_index && mode == "augmented") {
      e$long <- merge(e$long, proxies[, c("hospital_block_id", historical_proxy_names,
                                          "prior_history_short_lookback", "history_lookback_days"), with = FALSE],
                      by = "hospital_block_id", all.x = TRUE, sort = FALSE)
      setorder(e$long, site, hospital_block_id, m)
      for (nm in historical_proxy_names) setnafill(e$long, "const", 0, cols = nm)
    }
    eval(expr, envir = e)
    if (j == base_index && mode == "augmented") {
      e$baseline_covariates <- c(e$baseline_covariates, historical_proxy_names,
                                 "prior_history_short_lookback", "history_lookback_days")
    }
  }
  writeLines(c(paste("mode:", mode), paste("grace_hours:", grace),
               paste("history_coverage_start:", coverage_text), paste("history_mode:", history_mode),
               paste("history_source:", normalizePath(history_dir)),
               "Proxies describe events observed in the supplied source, not complete lifetime or outside-system history.",
               "Short windows: observed mode retains observed positives; complete mode zero-fills unknowns with a short-lookback indicator.",
               "Models also adjust for available lookback days; extract filters may make proxies redundant.",
               "Historical proxies augment baseline models; this is not NN/MNN effect calibration."),
             file.path(output, "sensitive", mode, "historical_proxy_design.txt"))
}
# Small-cell-controlled history summaries, with nested/complementary counts handled.
threshold <- as.integer(Sys.getenv("SUMMARY_SMALL_CELL_THRESHOLD", "10"))
if (!is.finite(threshold) || threshold < 1) stop("Invalid suppression threshold.")
summary <- rbindlist(lapply(historical_proxy_names, function(nm) {
  x <- proxies[!is.na(get(nm)), get(nm)]
  report <- sum(x) >= threshold && sum(1L - x) >= threshold
  data.table(variable = nm, history_mode = history_mode, extract_calendar_start = coverage_text, zero_interpretation = "no qualifying event observed in the supplied extract/window", n_evaluable = if (length(x) >= threshold) length(x) else NA_integer_,
             n_positive = if (report) sum(x) else NA_integer_,
             percent_positive = if (report) 100 * mean(x) else NA_real_)
}))
fwrite(summary, file.path(output, "shareable", "historical_proxy_summary.csv"))
ref <- fread(file.path(output, "shareable", "reference", "five_method_tidy_estimates.csv"))
aug <- fread(file.path(output, "shareable", "augmented", "five_method_tidy_estimates.csv"))
comparison <- merge(ref, aug, by = c("method", "estimand"), suffixes = c("_reference", "_augmented"))
fwrite(comparison, file.path(output, "shareable", "reference_vs_augmented.csv"))
cat("Paired historical-proxy sensitivity analysis saved to ", output, "\n", sep = "")
