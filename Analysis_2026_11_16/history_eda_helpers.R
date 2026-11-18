build_history_eda <- function(cohort,crosswalk,hospitalizations,coverage_start=as.Date(NA)) {
  c <- copy(as.data.table(cohort));cw <- copy(as.data.table(crosswalk))
  h <- copy(as.data.table(hospitalizations))
  for(v in c("patient_id","hospitalization_id"))h[,(v):=as.character(get(v))]
  for(v in c("admission_dttm","discharge_dttm"))h[,(v):=as.POSIXct(get(v),tz="UTC")]
  c[,patient_id:=as.character(patient_id)]
  if(anyDuplicated(c$patient_id))stop("Expected one index block per patient")
  valid_id <- !is.na(h$hospitalization_id)&nzchar(h$hospitalization_id)
  versions <- h[valid_id,.(versions=uniqueN(.SD)),by=hospitalization_id,
    .SDcols=c("patient_id","admission_dttm","discharge_dttm","admission_type_category")]
  if(any(versions$versions>1L,na.rm=TRUE))stop("Conflicting duplicate hospitalization IDs")
  source_rows <- nrow(h)
  h <- rbind(h[!valid_id],unique(h[valid_id],by="hospitalization_id"))
  if(!all(c$patient_id%in%h$patient_id)||!all(cw$hospitalization_id%in%h$hospitalization_id))
    stop("Expanded table does not contain all cohort/index encounters")
  valid <- !is.na(h$hospitalization_id)&nzchar(h$hospitalization_id)&
    !is.na(h$admission_dttm)&!is.na(h$discharge_dttm)&h$admission_dttm<h$discharge_dttm
  if(!all(cw$hospitalization_id%in%h$hospitalization_id[valid]))stop("Index encounters have invalid dates/IDs")
  components <- assign_patient_blocks_base(h[valid],6)
  setnames(components,"hospital_block_id","expanded_episode_id")
  setorder(components,patient_id,expanded_episode_id,component_order)
  index_ids <- unique(components[hospitalization_id%in%cw$hospitalization_id,expanded_episode_id])
  episodes <- components[,.(episode_admission=min(admission_dttm),episode_discharge=max(discharge_dttm),
    component_encounters=.N,first_admission_type=as.character(admission_type_category[1L])),
    by=.(patient_id,expanded_episode_id)]
  episodes[,ed_origin:=fcase(tolower(trimws(first_admission_type))=="ed",TRUE,
    tolower(trimws(first_admission_type))%in%c("direct","osh","facility"),FALSE,default=NA)]
  c[,index_admission:=as.POSIXct(admission_dttm,tz="UTC")]
  episodes <- merge(episodes,c[,.(patient_id,index_admission)],by="patient_id")
  episodes[,index_connected:=expanded_episode_id%in%index_ids]
  episodes[,completed_prior:=!index_connected&episode_admission<index_admission&episode_discharge<index_admission]
  prior <- episodes[completed_prior==TRUE]
  prior[,`:=`(los_days=as.numeric(difftime(episode_discharge,episode_admission,units="days")),
    within_year=episode_admission>=index_admission-365*86400,
    days_before_index=as.numeric(difftime(index_admission,episode_admission,units="days")),
    overlap_year_days=pmax(0,as.numeric(difftime(pmin(episode_discharge,index_admission),
      pmax(episode_admission,index_admission-365*86400),units="days"))))]
  bad <- merge(h[!valid],c[,.(patient_id,index_admission)],by="patient_id")
  result <- list()
  for(window in c("All available prior history","Previous 365 days"))for(type in c("All hospitalizations","ED-origin hospitalizations")) {
    z <- if(window=="Previous 365 days")prior[within_year==TRUE] else copy(prior)
    uncertain_ed <- z[is.na(ed_origin),.N,by=patient_id]
    if(type=="ED-origin hospitalizations")z <- z[ed_origin==TRUE]
    totals <- z[,.(observed_episode_count=.N,total_los_days=sum(los_days),
      mean_los_days=mean(los_days),longest_los_days=max(los_days),
      days_since_last_discharge=min(as.numeric(difftime(index_admission,episode_discharge,units="days")))),by=patient_id]
    x <- merge(c,totals,by="patient_id",all.x=TRUE,sort=FALSE)
    for(v in c("observed_episode_count","total_los_days","longest_los_days"))x[is.na(get(v)),(v):=0]
    potential <- bad[is.na(admission_dttm)|admission_dttm<index_admission]
    if(window=="Previous 365 days")potential <- potential[is.na(admission_dttm)|admission_dttm>=index_admission-365*86400]
    x[,`:=`(window=window,encounter_type=type,invalid_possible_prior_rows=0L,unknown_ed_episodes=0L)]
    at <- match(x$patient_id,uncertain_ed$patient_id)
    if(type=="ED-origin hospitalizations")x[,unknown_ed_episodes:=fifelse(is.na(at),0L,uncertain_ed$N[at])]
    nbad <- potential[,.N,by=patient_id];at <- match(x$patient_id,nbad$patient_id)
    x[,invalid_possible_prior_rows:=fifelse(is.na(at),0L,nbad$N[at])]
    x[,year_coverage_verified:=if(is.na(coverage_start))NA else as.Date(index_admission-365*86400)>=coverage_start]
    result[[length(result)+1L]] <- x
  }
  # Calendar-window duration includes older completed episodes overlapping the year.
  overlap <- rbindlist(lapply(c("All hospitalizations","ED-origin hospitalizations"),function(type) {
    z <- prior;if(type=="ED-origin hospitalizations")z <- z[ed_origin==TRUE]
    z[,.(observed_days_in_previous_year=sum(overlap_year_days)),by=patient_id][,encounter_type:=type]
  }))
  audit <- data.table(measure=c("source_rows","duplicate_rows_removed","invalid_date_or_id_rows",
    "cohort_patients","merged_episodes","index_connected_episodes_excluded","completed_prior_episodes",
    "nonprior_nonindex_episodes_excluded","history_coverage_start"),
    value=as.character(c(source_rows,source_rows-nrow(h),sum(!valid),nrow(c),nrow(episodes),
      sum(episodes$index_connected),nrow(prior),sum(!episodes$index_connected&!episodes$completed_prior),as.character(coverage_start))))
  list(patient=rbindlist(result),prior=prior,overlap=overlap,audit=audit,
    invalid=bad,mapping=h[,.(rows=.N),by=admission_type_category])
}

history_summary <- function(x) {
  finite <- x[is.finite(x)]
  if(!length(finite))return(list(n=0L,mean=NA_real_,sd=NA_real_,min=NA_real_,p25=NA_real_,median=NA_real_,p75=NA_real_,p95=NA_real_,max=NA_real_))
  q <- as.numeric(quantile(finite,c(0,.25,.5,.75,.95,1)))
  list(n=length(finite),mean=mean(finite),sd=sd(finite),min=q[1],p25=q[2],median=q[3],p75=q[4],p95=q[5],max=q[6])
}
