# Primary NCO: a distinct within-system admission starting in the 365 days
# before the selected index BLOCK admission. Index-connected episodes excluded.
build_prior_admission_nco <- function(cohort,crosswalk,hosp,coverage_start=as.Date(NA)) {
 c<-copy(as.data.table(cohort));h<-copy(as.data.table(hosp));cw<-copy(as.data.table(crosswalk))
 for(v in c("patient_id","hospitalization_id"))h[,(v):=as.character(get(v))]
 for(v in c("admission_dttm","discharge_dttm"))h[,(v):=as.POSIXct(get(v),tz="UTC")]
 if(anyDuplicated(h$hospitalization_id)) {
  conflict<-h[,.(versions=uniqueN(paste(patient_id,admission_dttm,discharge_dttm))),by=hospitalization_id][versions>1]
  if(nrow(conflict))stop("Conflicting hospitalization IDs in expanded table; deduplicate source first")
  h<-unique(h,by="hospitalization_id")
 }
 if(!all(as.character(c$patient_id)%in%h$patient_id))stop("Expanded table does not contain all analyzed patients")
 if(!all(as.character(cw$hospitalization_id)%in%h$hospitalization_id))stop("Expanded table lacks selected index component encounter IDs")
 valid<-h[!is.na(admission_dttm)&!is.na(discharge_dttm)&admission_dttm<discharge_dttm]
 episodes<-assign_patient_blocks_base(valid,6)
 setnames(episodes,"hospital_block_id","expanded_episode_id")
 index_episodes<-unique(episodes[hospitalization_id%in%cw$hospitalization_id,expanded_episode_id])
 episodes<-episodes[,.(episode_admission=min(admission_dttm),episode_discharge=max(discharge_dttm)),
   by=.(patient_id,expanded_episode_id)]
 c[,`:=`(window_end=as.POSIXct(admission_dttm,tz="UTC"),window_start=as.POSIXct(admission_dttm,tz="UTC")-365*86400)]
 history<-merge(episodes,c[,.(patient_id,hospital_block_id,window_start,window_end)],by="patient_id",allow.cartesian=TRUE)
 history<-history[!expanded_episode_id%in%index_episodes & episode_admission>=window_start &
   episode_admission<window_end & episode_discharge<window_end]
 c[,prior_year_admission:=as.integer(patient_id%in%history$patient_id)]
 invalid<-h[is.na(admission_dttm)|is.na(discharge_dttm)|admission_dttm>=discharge_dttm]
 invalid<-merge(invalid,c[,.(patient_id,window_start,window_end)],by="patient_id")
 uncertain<-invalid[is.na(admission_dttm)|(admission_dttm>=window_start&admission_dttm<window_end),unique(patient_id)]
 c[,history_date_uncertain:=patient_id%in%uncertain]
 c[prior_year_admission==0L & history_date_uncertain,prior_year_admission:=NA_integer_]
 c[,coverage_verified:=if(is.na(coverage_start))NA else as.Date(window_start)>=coverage_start]
 c[,nco_eligible:=!is.na(prior_year_admission) & (is.na(coverage_verified)|coverage_verified)]
 list(cohort=c,episodes=history,source_summary=data.table(measure=c("expanded_patient_rows","expanded_patients",
  "invalid_date_rows","verified_coverage_start"),value=c(nrow(h),uniqueN(h$patient_id),nrow(invalid),as.character(coverage_start))))
}

# Patient-level influence contributions combine both clones before variance:
# the two arms can contain the same original patient.
nco_weighted_contrast <- function(y,w1,w0) {
 keep<-is.finite(y);w1<-ifelse(keep&is.finite(w1)&w1>0,w1,0);w0<-ifelse(keep&is.finite(w0)&w0>0,w0,0)
 y[!keep]<-0
 n1<-sum(w1>0);n0<-sum(w0>0)
 if(!sum(w1)||!sum(w0))return(data.table(n_early=n1,n_control=n0,p_early=NA_real_,p_control=NA_real_,
  difference=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_,ess_early=NA_real_,ess_control=NA_real_,status="Empty arm"))
 p1<-sum(w1*y)/sum(w1);p0<-sum(w0*y)/sum(w0)
 influence<-w1*(y-p1)/sum(w1)-w0*(y-p0)/sum(w0)
 se<-sqrt(sum(influence^2));d<-p1-p0
 result<-data.table(n_early=n1,n_control=n0,p_early=p1,p_control=p0,difference=d,se=se,
  lower=d-qnorm(.975)*se,upper=d+qnorm(.975)*se,
  ess_early=sum(w1)^2/sum(w1^2),ess_control=sum(w0)^2/sum(w0^2),
  status=if(length(unique(y[keep&(w1+w0)>0]))<2)"NCO has no variation" else "ok")
 if(result$status!="ok")result[,`:=`(se=NA_real_,lower=NA_real_,upper=NA_real_)]
 result
}
