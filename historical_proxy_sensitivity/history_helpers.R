# Pure historical-feature construction, independently testable without Arrow.
historical_proxy_names <- c("prior_icu", "prior_respiratory_hospitalization",
                            "prior_emergency_hospitalization", "prior_imv")

build_history_proxies <- function(cohort, hospitalization, adt, respiratory,
                                  diagnosis, coverage_start, gap_hours = 6) {
  library(data.table)
  c <- copy(as.data.table(cohort))
  h <- copy(as.data.table(hospitalization))
  required <- c("hospitalization_id", "patient_id", "admission_dttm",
                "discharge_dttm", "admission_type_category")
  stopifnot(all(required %in% names(h)), !anyDuplicated(h$hospitalization_id),
            !anyDuplicated(c$hospital_block_id), is.finite(as.numeric(coverage_start)))
  for (nm in c("admission_dttm", "discharge_dttm")) {
    h[, (nm) := as.POSIXct(get(nm), tz = "UTC")]
  }
  h[, `:=`(patient_id = as.character(patient_id),
           hospitalization_id = as.character(hospitalization_id))]
  c[, `:=`(patient_id = as.character(patient_id),
           admission_dttm = as.POSIXct(admission_dttm, tz = "UTC"))]
  h <- h[patient_id %in% c$patient_id]
  if (h[is.na(admission_dttm) | is.na(discharge_dttm) |
        discharge_dttm < admission_dttm, .N]) {
    stop("Historical hospitalizations have invalid admission/discharge times.")
  }
  h[, `:=`(icu = 0L, imv = 0L, resp = 0L,
            emergency = as.integer(tolower(trimws(admission_type_category)) == "ed"))]
  if (anyNA(h$emergency)) stop("Missing historical admission type; cannot infer absence of emergency admission.")
  mark_events <- function(events, time_name, value_name, values) {
    e <- copy(as.data.table(events))
    e[, hospitalization_id := as.character(hospitalization_id)]
    e <- merge(e, h[, .(hospitalization_id, admission_dttm, discharge_dttm)],
               by = "hospitalization_id")
    e[, event_time := as.POSIXct(get(time_name), tz = "UTC")]
    unique(e[!is.na(event_time) & event_time >= admission_dttm &
               event_time <= discharge_dttm &
               tolower(trimws(get(value_name))) %in% values, hospitalization_id])
  }
  h[hospitalization_id %in% mark_events(adt, "in_dttm", "location_category", "icu"), icu := 1L]
  h[hospitalization_id %in% mark_events(respiratory, "recorded_dttm", "device_category", "imv"), imv := 1L]
  d <- copy(as.data.table(diagnosis))
  d[, code := gsub("[^A-Z0-9]", "", toupper(as.character(diagnosis_code)))]
  d[, fmt := gsub("[^A-Z0-9]", "", toupper(as.character(diagnosis_code_format)))]
  # Any diagnosis in the completed prior episode; not restricted to POA.
  ids <- d[(fmt %in% c("ICD10CM", "ICD10", "10") & grepl("^J(1[2-8]|96)", code)) |
             (fmt %in% c("ICD9CM", "ICD9", "9") &
                grepl("^(48[0-6]|51881|51882|51884)", code)), hospitalization_id]
  h[hospitalization_id %in% as.character(ids), resp := 1L]
  # Merge transfers/overlaps first, so current-episode components cannot be history.
  setorder(h, patient_id, admission_dttm, discharge_dttm)
  h[, previous_end := shift(cummax(as.numeric(discharge_dttm))), by = patient_id]
  h[, episode := cumsum(is.na(previous_end) |
                         as.numeric(admission_dttm) - previous_end > gap_hours * 3600),
    by = patient_id]
  episodes <- h[, .(history_start = min(admission_dttm),
                     history_end = max(discharge_dttm),
                     prior_icu = max(icu), prior_respiratory_hospitalization = max(resp),
                     prior_emergency_hospitalization = max(emergency), prior_imv = max(imv)),
                 by = .(patient_id, episode)]
  out <- c[, .(hospital_block_id, patient_id, index_admission = admission_dttm)]
  out[, history_available := as.integer(index_admission - 365 * 86400 >= coverage_start)]
  joined <- merge(out, episodes, by = "patient_id", allow.cartesian = TRUE)
  prior <- joined[history_start >= index_admission - 365 * 86400 &
                    history_end < index_admission - gap_hours * 3600]
  flags <- prior[, lapply(.SD, function(x) if (length(x)) max(x) else NA_integer_), by = hospital_block_id,
                 .SDcols = historical_proxy_names]
  out <- merge(out, flags, by = "hospital_block_id", all.x = TRUE, sort = FALSE)
  for (nm in historical_proxy_names) {
    out[history_available == 1L & is.na(get(nm)), (nm) := 0L]
    out[history_available == 0L, (nm) := NA_integer_]
  }
  out[, prior_history_unavailable := 1L - history_available]
  out
}
