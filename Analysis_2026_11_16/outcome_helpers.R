# Shared primary endpoint for EDA and every estimator in this dated pipeline.
PRIMARY_OUTCOME <- "In-hospital death or hospice discharge by day 28"
PRIMARY_EVENT_LABEL <- "In-hospital death or hospice discharge"

derive_primary_outcome <- function(cohort, time_zero_col="respiratory_failure_start") {
  x <- data.table::copy(data.table::as.data.table(cohort))
  required <- c("admission_dttm","discharge_dttm","discharge_category",time_zero_col)
  if(!all(required %in% names(x)))stop("Missing outcome inputs: ",paste(setdiff(required,names(x)),collapse=", "))
  utc <- function(z)as.POSIXct(z,tz="UTC")
  if(!"death_dttm" %in% names(x))x[,death_dttm:=as.POSIXct(NA,tz="UTC")]
  x[,`:=`(death_dttm=utc(death_dttm),admission_dttm=utc(admission_dttm),discharge_dttm=utc(discharge_dttm))]
  x[,outcome_discharge_category:=tolower(trimws(as.character(discharge_category))) ]
  in_block <- !is.na(x$death_dttm)&!is.na(x$admission_dttm)&!is.na(x$discharge_dttm)&
    x$death_dttm>=x$admission_dttm & x$death_dttm<=x$discharge_dttm
  death <- rep(NA_real_,nrow(x));death[in_block] <- as.numeric(x$death_dttm[in_block])
  expired <- !is.na(x$outcome_discharge_category)&x$outcome_discharge_category=="expired"
  fallback <- expired&!is.finite(death)
  death[fallback] <- as.numeric(x$discharge_dttm[fallback])
  hospice <- !is.na(x$outcome_discharge_category)&grepl("hospice",x$outcome_discharge_category)
  hospice_time <- rep(NA_real_,nrow(x));hospice_time[hospice] <- as.numeric(x$discharge_dttm[hospice])
  event <- pmin(death,hospice_time,na.rm=TRUE);event[!is.finite(event)] <- NA_real_
  kind <- ifelse(is.finite(event),ifelse(is.finite(death)&event==death,"in_hospital_death","hospice_discharge"),NA_character_)
  source <- ifelse(kind=="hospice_discharge","hospice discharge timestamp",
    ifelse(fallback,"expired discharge timestamp","patient death_dttm within selected hospital block"))
  x[,`:=`(in_hospital_death_dttm=as.POSIXct(death,origin="1970-01-01",tz="UTC"),
    hospice_discharge_dttm=as.POSIXct(hospice_time,origin="1970-01-01",tz="UTC"),
    event_dttm=as.POSIXct(event,origin="1970-01-01",tz="UTC"),event_type=kind,event_source=source,
    post_discharge_death_ignored=!is.na(death_dttm)&!is.na(discharge_dttm)&death_dttm>discharge_dttm)]
  x[,event_hours:=as.numeric(difftime(event_dttm,get(time_zero_col),units="hours"))]
  x[!is.finite(event_hours)|event_hours<0,`:=`(event_dttm=as.POSIXct(NA,tz="UTC"),
    event_hours=NA_real_,event_type=NA_character_,event_source=NA_character_)]
  x
}

write_outcome_definition <- function(x,out,cutoff) {
  definition <- data.table::data.table(parameter=c("outcome","time_origin","in_hospital_death",
   "hospice_discharge","post_discharge_death","other_alive_discharge","administrative_cutoff","scope"),
   value=c(PRIMARY_OUTCOME,"respiratory-failure start", "Death timestamp within selected block; otherwise Expired-discharge time",
    "Discharge time for any discharge category containing hospice",
    "Excluded even when a later death timestamp is recorded",
    "No primary event; not censored at discharge; remains event-free through day 28 unless administratively censored",
    format(cutoff,"%Y-%m-%d"),"Selected hospital block only; subsequent separate readmissions are not added"))
  data.table::fwrite(definition,file.path(out,"outcome_definition.csv"),na="")
  x <- data.table::copy(x)
  x[,event_observed_by28:=is.finite(event_hours)&event_hours<=672&event_dttm<=cutoff]
  data.table::fwrite(x[,.(patients=.N,events_by28=sum(event_observed_by28),
    post_discharge_death_timestamps_ignored=sum(post_discharge_death_ignored)),
    by=.(discharge_category,event_type)],file.path(out,"outcome_ascertainment_audit.csv"),na="")
}
