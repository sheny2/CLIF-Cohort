# Historical NCOs: starts in the preceding year, completed before index
# admission; reuse 09's episode construction and uncertainty/coverage rules.
build_historical_nco_panel <- function(cohort,crosswalk,hosp,coverage_start=as.Date(NA)) {
  built <- build_prior_admission_nco(cohort,crosswalk,hosp,coverage_start)
  c <- built$cohort; history <- built$episodes; h <- copy(as.data.table(hosp))
  for (v in c("patient_id","hospitalization_id")) h[,(v):=as.character(get(v))]
  for (v in c("admission_dttm","discharge_dttm")) h[,(v):=as.POSIXct(get(v),tz="UTC")]
  optional <- intersect(c("admission_type_category","discharge_category"),names(h))
  if(length(optional)) {
    conflicts <- h[,.(versions=uniqueN(.SD)),by=hospitalization_id,.SDcols=optional]
    if(any(conflicts$versions>1L))stop("Conflicting historical encounter category fields")
  }
  h <- unique(h,by="hospitalization_id")
  valid <- h[!is.na(admission_dttm)&!is.na(discharge_dttm)&admission_dttm<discharge_dttm]
  components <- assign_patient_blocks_base(valid,6)
  setnames(components,"hospital_block_id","expanded_episode_id")
  setorder(components,patient_id,expanded_episode_id,component_order)
  clean <- function(x)tolower(trimws(gsub("[_-]+"," ",as.character(x))))
  flag <- function(x,pattern) {
    z<-clean(x);answer<-as.integer(grepl(pattern,z))
    answer[is.na(z)|z%in%c("","unknown","unspecified","other","na","n/a")]<-NA_integer_
    answer
  }
  metadata <- components[,.(
    ed=if("admission_type_category"%in%optional)
      flag(admission_type_category[1L],"^(ed|emergency|emergency department)$") else NA_integer_,
    snf_rehab=if("discharge_category"%in%optional)
      flag(discharge_category[.N],"skilled nursing|nursing facility|rehab|^snf$") else NA_integer_
  ),by=.(patient_id,expanded_episode_id)]
  history<-merge(history,metadata,by=c("patient_id","expanded_episode_id"),all.x=TRUE)
  history[,los_days:=as.numeric(difftime(episode_discharge,episode_admission,units="days"))]
  # Absence of usable history is not evidence that every encounter category
  # was observed; three-valued logic preserves unknowns unless a positive exists.
  any_flag<-function(z)if(any(z==1L,na.rm=TRUE))1L else if(anyNA(z))NA_integer_ else 0L
  counts<-history[,.(prior_episodes=.N,prolonged=as.integer(any(los_days>=7)),
    ed=any_flag(ed),snf_rehab=any_flag(snf_rehab)),by=patient_id]
  at<-match(c$patient_id,counts$patient_id)
  c[,prior_episode_count:=fifelse(is.na(at),0L,counts$prior_episodes[at])]
  c[,nco_any_admission:=prior_year_admission]
  c[,nco_recurrent_admission:=as.integer(prior_episode_count>=2L)]
  for(nm in c("prolonged","ed","snf_rehab")) {
    z<-counts[[nm]][at];z[is.na(at)]<-0L
    c[,(paste0("nco_",nm)):=z]
  }
  cols<-c("nco_any_admission","nco_recurrent_admission","nco_prolonged","nco_ed","nco_snf_rehab")
  for(nm in cols) {
    # Incomplete dates can hide an additional episode or a qualifying feature.
    c[history_date_uncertain & (is.na(get(nm))|get(nm)==0L),(nm):=NA_integer_]
    c[!is.na(coverage_verified)&!coverage_verified,(nm):=NA_integer_]
  }
  if(!"admission_type_category"%in%optional)c[,nco_ed:=NA_integer_]
  if(!"discharge_category"%in%optional)c[,nco_snf_rehab:=NA_integer_]
  definitions<-data.table(nco=cols,label=c("Any prior-year admission","At least two prior-year admissions",
    "Prior-year stay at least 7 days","Prior-year ED-origin admission","Prior-year SNF/rehab discharge"),
    required_field=c("timestamps","timestamps","timestamps","admission_type_category","discharge_category"),
    field_available=c(TRUE,TRUE,TRUE,"admission_type_category"%in%optional,"discharge_category"%in%optional))
  mapping<-rbindlist(lapply(optional,function(v)h[,.(rows=.N),by=.(source_value=as.character(get(v)))][,field:=v]),fill=TRUE)
  list(cohort=c,history=history,definitions=definitions,mapping=mapping,source_summary=built$source_summary)
}

historical_log_pr <- function(y,w1,w0) {
  d<-nco_weighted_contrast(y,w1,w0)
  ok<-is.finite(y);w1<-ifelse(ok&is.finite(w1)&w1>0,w1,0);w0<-ifelse(ok&is.finite(w0)&w0>0,w0,0)
  d[,`:=`(positive_early=sum(y[w1>0]==1),negative_early=sum(y[w1>0]==0),
    positive_control=sum(y[w0>0]==1),negative_control=sum(y[w0>0]==0),
    log_pr=NA_real_,se_log_pr=NA_real_,prevalence_ratio=NA_real_,pr_lower=NA_real_,pr_upper=NA_real_)]
  if(d$status!="ok"||!is.finite(d$p_early)||!is.finite(d$p_control)||d$p_early<=0||d$p_control<=0)return(d)
  y[!ok]<-0
  influence<-w1*(y-d$p_early)/sum(w1)/d$p_early-w0*(y-d$p_control)/sum(w0)/d$p_control
  est<-log(d$p_early/d$p_control);log_pr_se_value<-sqrt(sum(influence^2))
  # d already has a column named se (RD scale). A distinct local name avoids
  # data.table resolving that column instead of the log-ratio SE.
  d[,`:=`(log_pr=est,se_log_pr=log_pr_se_value,prevalence_ratio=exp(est),
    pr_lower=exp(est-1.96*log_pr_se_value),pr_upper=exp(est+1.96*log_pr_se_value))]
  d
}
