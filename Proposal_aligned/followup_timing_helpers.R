# Six-hour descriptive counts; each patient contributes at most one event per type.
build_followup_timing <- function(cohort, administrations, cutoff, min_cell_n = 10L) {
  library(data.table)
  c <- copy(as.data.table(cohort))
  a <- copy(as.data.table(administrations))
  stopifnot(!anyDuplicated(c$hospital_block_id), min_cell_n >= 1L)
  c[, `:=`(hospital_block_id = as.character(hospital_block_id),
            respiratory_failure_start = as.POSIXct(respiratory_failure_start, tz = "UTC"),
            event_dttm = as.POSIXct(event_dttm, tz = "UTC"),
            discharge_dttm = as.POSIXct(discharge_dttm, tz = "UTC"))]
  c[, death_hours := as.numeric(difftime(event_dttm, respiratory_failure_start, units = "hours"))]
  death <- c[is.finite(death_hours) & death_hours >= 0 & death_hours <= 672 &
               event_dttm <= cutoff, .(hours = death_hours)]
  a[, `:=`(hospital_block_id = as.character(hospital_block_id),
            admin_dttm = as.POSIXct(admin_dttm, tz = "UTC"))]
  a <- merge(a[, .(hospital_block_id, admin_dttm)], c[, .(
    hospital_block_id, respiratory_failure_start, event_dttm, discharge_dttm)],
    by = "hospital_block_id")
  a[, hours := as.numeric(difftime(admin_dttm, respiratory_failure_start, units = "hours"))]
  a <- a[is.finite(hours) & hours >= 0 & hours <= 672 & admin_dttm <= cutoff &
           !is.na(discharge_dttm) & admin_dttm <= discharge_dttm &
           (is.na(event_dttm) | admin_dttm <= event_dttm)]
  # Do not use the grace-limited treatment assignment or the first-ever dose:
  # a pre-time-zero dose must not hide a later qualifying administration.
  first <- a[, .(hours = if (length(hours)) min(hours) else NA_real_), by = hospital_block_id]
  events <- rbindlist(list(data.table(event_type = rep("Death", nrow(death)), hours = death$hours),
                          data.table(event_type = rep("First qualifying steroid administration", nrow(first)), hours = first$hours)))
  events[, interval := pmax(1L, as.integer(ceiling(hours / 6)))]
  counts <- events[, .(n = .N), by = .(event_type, interval)]
  bins <- CJ(event_type = c("Death", "First qualifying steroid administration"), interval = 1:112)
  bins <- merge(bins, counts, by = c("event_type", "interval"), all.x = TRUE)
  bins[is.na(n), n := 0L]
  bins[, `:=`(interval_start_hours = (interval - 1L) * 6,
               interval_end_hours = interval * 6,
               n_reportable = fifelse(n > 0 & n < min_cell_n, NA_integer_, n))]
  # Complementary suppression prevents the 24/96-hour totals from revealing a
  # single hidden six-hour cell by subtraction within either grace window.
  bins[, segment := fifelse(interval_end_hours <= 24, 1L,
                             fifelse(interval_end_hours <= 96, 2L, 3L))]
  bins[, n_reportable := {
    z <- n_reportable
    if (sum(is.na(z)) == 1L && any(!is.na(z) & z > 0L)) {
      candidates <- which(!is.na(z) & z > 0L)
      z[candidates[which.min(z[candidates])]] <- NA_integer_
    }
    z
  }, by = .(event_type, segment)]
  totals <- rbindlist(lapply(c(24, 96), function(g) {
    data.table(grace_hours = g, event_type = unique(bins$event_type))[, n :=
      vapply(event_type, function(type) events[event_type == type & hours <= g, .N], integer(1))]
  }))
  totals[, `:=`(n_reportable = fifelse(n > 0 & n < min_cell_n, NA_integer_, n),
                 cohort_n = nrow(c))]
  # If a segment has only one positive small cell and no complementary cell,
  # hide affected cumulative totals too (otherwise 96h minus 24h reveals it).
  unresolved <- bins[segment <= 2L, .(hidden = sum(is.na(n_reportable))),
                     by = .(event_type, segment)][hidden == 1L]
  if (nrow(unresolved)) for (j in seq_len(nrow(unresolved))) {
    type <- unresolved$event_type[j]
    lower <- if (unresolved$segment[j] == 1L) 24 else 96
    totals[event_type == type & grace_hours >= lower, n_reportable := NA_integer_]
  }
  totals[, count_suppressed := is.na(n_reportable)]
  totals[, percent_of_cohort := fifelse(is.na(n_reportable) | cohort_n < min_cell_n,
                                       NA_real_, 100 * n_reportable / cohort_n)]
  bins[, `:=`(count_suppressed = is.na(n_reportable), small_cell_threshold = min_cell_n,
               interval_definition = "First bin [0,6]; subsequent bins (start,end]")]
  totals[, small_cell_threshold := min_cell_n]
  list(intervals = bins[, !c("n", "segment"), with = FALSE],
       grace_comparison = totals[, !"n", with = FALSE])
}
