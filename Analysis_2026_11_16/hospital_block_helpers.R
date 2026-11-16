# Hospital-block assignment without native data.table grouping or sorting.
assign_blocks_running_end <- function(admission, discharge, threshold_hours) {
  n <- length(admission)
  block <- integer(n)
  gap <- rep(NA_real_, n)
  running_end <- rep(as.POSIXct(NA, tz = "UTC"), n)
  current_block <- 0L
  block_end <- as.POSIXct(NA, tz = "UTC")
  for (i in seq_len(n)) {
    if (i > 1L) {
      gap[i] <- as.numeric(difftime(admission[i], block_end, units = "hours"))
    }
    begins_new <- i == 1L || is.na(gap[i]) || gap[i] > threshold_hours
    if (begins_new) {
      current_block <- current_block + 1L
      block_end <- discharge[i]
    } else if (!is.na(discharge[i]) && discharge[i] > block_end) {
      block_end <- discharge[i]
    }
    block[i] <- current_block
    running_end[i] <- block_end
  }
  list(block = block, gap = gap, running_end = running_end)
}

assign_patient_blocks_base <- function(x, threshold_hours) {
  x <- as.data.frame(x)
  if(!nrow(x))stop("No valid completed hospitalizations to construct blocks")
  o <- order(as.character(x$patient_id),as.numeric(x$admission_dttm),
             as.numeric(x$discharge_dttm),as.character(x$hospitalization_id),method="radix")
  x <- x[o,,drop=FALSE]
  patient <- as.character(x$patient_id)
  starts <- c(1L,which(patient[-1L]!=patient[-length(patient)])+1L)
  ends <- c(starts[-1L]-1L,nrow(x))
  block <- component <- integer(nrow(x));gap <- running_end <- rep(NA_real_,nrow(x))
  for(j in seq_along(starts)) {
    ii <- seq.int(starts[j],ends[j])
    a <- assign_blocks_running_end(x$admission_dttm[ii],x$discharge_dttm[ii],threshold_hours)
    block[ii] <- a$block;gap[ii] <- a$gap;running_end[ii] <- as.numeric(a$running_end)
    component[ii] <- sequence(rle(a$block)$lengths)
  }
  x$hospital_block_number <- block
  x$gap_hours_from_running_block_end <- gap
  x$running_block_end_dttm <- as.POSIXct(running_end,origin="1970-01-01",tz="UTC")
  x$hospital_block_id <- paste0(patient,"_block_",block)
  x$component_order <- component
  data.table::as.data.table(x)
}
