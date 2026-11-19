build_steroid_start_audit <- function(cohort,qualifying,administrations,design) {
  c <- data.table::copy(data.table::as.data.table(cohort))
  q <- data.table::copy(data.table::as.data.table(qualifying))
  a <- data.table::copy(data.table::as.data.table(administrations))
  cutoff <- as.POSIXct(paste(design$data_cutoff,"23:59:59"),tz="UTC")
  if(is.na(cutoff))stop("Invalid saved data cutoff")
  if(anyDuplicated(c$hospital_block_id)||anyDuplicated(c$patient_id))stop("Expected one index block per patient")
  prepare <- function(z) {
    z[,admin_dttm:=as.POSIXct(admin_dttm,tz="UTC")]
    z <- z[hospital_block_id%in%c$hospital_block_id]
    z <- merge(z[,.(hospital_block_id,admin_dttm,med_category_clean,med_route_clean,dose_mg,
      mar_action_category)],c[,.(hospital_block_id,index_dttm,admission_dttm,discharge_dttm,event_dttm)],by="hospital_block_id")
    z <- z[tolower(trimws(as.character(mar_action_category)))=="given" & !is.na(admin_dttm) &
      admin_dttm<=cutoff & admin_dttm<=discharge_dttm & (is.na(event_dttm)|admin_dttm<=event_dttm)]
    z[,hours_from_zero:=as.numeric(difftime(admin_dttm,index_dttm,units="hours"))]
    z
  }
  q <- prepare(q)
  a <- a[route_eligible_for_drug==TRUE & is.finite(dose_mg) & dose_mg>0 & !is.na(minimum_dose_mg)]
  a <- prepare(a)
  pre <- q[hours_from_zero>= -24 & hours_from_zero<0]
  pre_summary <- pre[,.(prior_qualifying_doses_24h=.N,last_prior_hours=max(hours_from_zero),
    prior_qualifying_6h=any(hours_from_zero>= -6),prior_qualifying_12h=any(hours_from_zero>= -12),
    prior_dose_before_block_admission=any(admin_dttm<admission_dttm)),by=hospital_block_id]
  c <- merge(c,pre_summary,by="hospital_block_id",all.x=TRUE,sort=FALSE)
  c[is.na(prior_qualifying_doses_24h),prior_qualifying_doses_24h:=0L]
  for(v in c("prior_qualifying_6h","prior_qualifying_12h","prior_dose_before_block_admission"))c[is.na(get(v)),(v):=FALSE]
  c[,prior_qualifying_24h:=prior_qualifying_doses_24h>0L]
  c[,any_qualifying_before_zero:=hospital_block_id%in%q[hours_from_zero<0,hospital_block_id]]
  c[,any_given_systemic_24h:=hospital_block_id%in%a[hours_from_zero>= -24&hours_from_zero<0,hospital_block_id]]
  c[,nominal_pre_zero_hours:=pmin(24,pmax(0,as.numeric(difftime(index_dttm,admission_dttm,units="hours"))))]
  c[,full_24h_in_selected_block:=nominal_pre_zero_hours>=24]
  post <- q[hours_from_zero>=0 & hours_from_zero<=design$horizon_days*24]
  first <- post[,.(first_post_hours=min(hours_from_zero)),by=hospital_block_id]
  post <- merge(post,first,by="hospital_block_id")
  first_drugs <- post[hours_from_zero==first_post_hours,.(first_post_drugs=paste(sort(unique(med_category_clean)),collapse=";"),
    multiple_first_post_drugs=data.table::uniqueN(med_category_clean)>1L),by=hospital_block_id]
  c <- merge(c,first,by="hospital_block_id",all.x=TRUE,sort=FALSE)
  c <- merge(c,first_drugs,by="hospital_block_id",all.x=TRUE,sort=FALSE)
  c[,rebuilt_S:=ifelse(is.finite(first_post_hours),pmax(1L,as.integer(ceiling(first_post_hours/design$step_hours))),design$horizon_intervals+1L)]
  c[,saved_start_matches_source:=S==rebuilt_S]
  c[,interval1_starter:=S==1L]
  c[,start_group:=data.table::fcase(S==1L,"Interval 1 (0–6 h)",S<=design$tau_intervals,"Later within grace",
    S<=design$horizon_intervals,"After grace",default="No post-zero initiation")]
  c[,prior_status:=data.table::fcase(prior_qualifying_24h,"Observed qualifying prior dose",
    full_24h_in_selected_block,"No observed prior dose; 24 h within block",
    default="No observed prior dose; shorter block lookback")]
  c[,`:=`(same_drug_prior_24h=FALSE,prior_to_post_gap_hours=first_post_hours-last_prior_hours)]
  if(nrow(pre)) {
    drugs <- pre[,.(prior_drugs=list(unique(med_category_clean))),by=hospital_block_id]
    for(i in which(c$prior_qualifying_24h & !is.na(c$first_post_drugs))) {
      j <- match(c$hospital_block_id[i],drugs$hospital_block_id)
      c$same_drug_prior_24h[i] <- length(intersect(strsplit(c$first_post_drugs[i],";",fixed=TRUE)[[1L]],drugs$prior_drugs[[j]]))>0L
    }
  }
  list(cohort=c,qualifying=q,prior=pre)
}
