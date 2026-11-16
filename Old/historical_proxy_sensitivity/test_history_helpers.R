suppressPackageStartupMessages(library(data.table))
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(self)), "history_helpers.R"))
tm <- function(x) as.POSIXct(x, tz = "UTC")
c <- data.table(hospital_block_id = c("b1", "b2", "b3"), patient_id = c("p1", "p2", "p3"),
                admission_dttm = tm(c("2025-06-01", "2025-06-01", "2024-06-01")))
h <- data.table(hospitalization_id = c("old", "transfer", "current", "too_old", "future"),
                patient_id = "p1", admission_type_category = c("ed", "other", "ed", "ed", "ed"),
                admission_dttm = tm(c("2025-01-01", "2025-05-31", "2025-06-01", "2023-01-01", "2025-07-01")),
                discharge_dttm = tm(c("2025-01-03", "2025-06-01", "2025-06-10", "2023-01-03", "2025-07-03")))
a <- data.table(hospitalization_id = c("old", "future"), in_dttm = tm(c("2025-01-02", "2025-07-02")), location_category = "icu")
r <- data.table(hospitalization_id = c("old", "transfer"), recorded_dttm = tm(c("2025-01-02", "2025-05-31")), device_category = "imv")
d <- data.table(hospitalization_id = "old", diagnosis_code = "J96.01", diagnosis_code_format = "ICD10CM")
x <- build_history_proxies(c, h, a, r, d, tm("2024-01-01"), history_mode = "complete_history")
stopifnot(all(unlist(x[hospital_block_id == "b1", ..historical_proxy_names]) == 1),
          all(unlist(x[hospital_block_id == "b2", ..historical_proxy_names]) == 0),
          all(is.na(unlist(x[hospital_block_id == "b3", ..historical_proxy_names]))),
          x[hospital_block_id == "b3", prior_history_short_lookback] == 1L)
# Transfer-connected, future and too-old episodes must not become prior events.
y <- build_history_proxies(c, h[hospitalization_id != "old"], a, r, d, tm("2024-01-01"), history_mode = "complete_history")
stopifnot(all(unlist(y[hospital_block_id == "b1", ..historical_proxy_names]) == 0))
# An observation timestamp outside the recorded prior admission does not qualify.
a[, in_dttm := tm("2025-08-01")]
z <- build_history_proxies(c, h, a, r, d, tm("2024-01-01"), history_mode = "complete_history")
stopifnot(z[hospital_block_id == "b1", observed_prior_icu] == 0L)
# Refuse ambiguous duplicate encounter records.
stopifnot(inherits(try(build_history_proxies(c, rbind(h, h[1]), a, r, d, tm("2024-01-01"), history_mode = "complete_history"), silent = TRUE), "try-error"))
cat("PASS: prior-event flags, absence, unavailable history, transfer exclusion, temporal limits, duplicate rejection.\n")

# Filtered-extract mode preserves positive evidence in partial calendar windows.
partial <- copy(c)
partial[hospital_block_id == "b1", admission_dttm := tm("2025-03-01")]
obs <- build_history_proxies(partial, h, a, r, d, tm("2025-01-01"))
stopifnot(obs[hospital_block_id == "b1", observed_prior_imv] == 1L,
          obs[hospital_block_id == "b1", prior_history_short_lookback] == 1L,
          obs[hospital_block_id == "b1", history_lookback_days] == 59,
          obs[hospital_block_id == "b2", observed_prior_imv] == 0L,
          all(obs$history_mode == "observed_extract"),
          !anyNA(obs[, ..historical_proxy_names]))
stopifnot(inherits(try(build_history_proxies(c, h, a, r, d, tm("2024-01-01"), history_mode = "invalid"), silent = TRUE), "try-error"))
cat("PASS: filtered-extract labels, short-window positives, zero semantics, lookback days and mode validation.\n")
