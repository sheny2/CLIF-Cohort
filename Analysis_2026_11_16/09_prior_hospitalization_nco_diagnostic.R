#!/usr/bin/env Rscript
# Rscript 09_prior_hospitalization_nco_diagnostic.R RUN_ROOT EXPANDED_HOSPITALIZATION_PARQUET [SYSTEM_COVERAGE_START_DATE]
# Diagnostic only: no primary-effect calibration or refitting of nuisance models.
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
setDTthreads(1)
a<-commandArgs(TRUE)
if(length(a)<2)stop("Usage: 09 RUN_ROOT EXPANDED_HOSPITALIZATION_PARQUET [SYSTEM_COVERAGE_START_DATE]")
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);root<-dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE)))
for(h in c("06_generalized_five_methods.R","hospital_block_helpers.R","nco_history_helpers.R"))source(file.path(root,h))
run<-normalizePath(a[1]);expanded<-normalizePath(a[2]);coverage<-if(length(a)>=3)as.Date(a[3]) else as.Date(NA)
if(length(a)>=3&&is.na(coverage))stop("Coverage start must be YYYY-MM-DD")
inp<-file.path(run,"sensitive","analysis")
if(!file.exists(file.path(run,"aggregate","analysis","COMPLETE.txt")))stop("Run 07 to completion first")
out<-file.path(run,"aggregate","nco_prior_hospitalization");sensitive<-file.path(run,"sensitive","nco_prior_hospitalization")
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
unlink(file.path(out,"COMPLETE.txt"))
write<-function(x,f)fwrite(x,file.path(out,f),na="")
c<-as.data.table(read_parquet(file.path(inp,"pneumonia_ccw_cohort.parquet")))
cw<-as.data.table(read_parquet(file.path(inp,"sampled_hospital_block_crosswalk.parquet")))
ids<-as.character(c$patient_id)
ds<-open_dataset(expanded);needed<-c("patient_id","hospitalization_id","admission_dttm","discharge_dttm")
if(!all(needed%in%ds$schema$names))stop("Expanded hospitalization table lacks required columns")
hosp<-ds |> filter(.data$patient_id%in%ids) |> select(all_of(needed)) |> collect() |> as.data.table()
constructed<-build_prior_admission_nco(c,cw,hosp,coverage);c<-constructed$cohort
write(constructed$source_summary,"historical_source_audit.csv")
write(c[,.(patients=.N,observed_prior_admission=sum(prior_year_admission==1L,na.rm=TRUE),
 unknown_nco=sum(is.na(prior_year_admission)),eligible_for_diagnostic=sum(nco_eligible),
 coverage_unverified=sum(is.na(coverage_verified)),short_system_coverage=sum(coverage_verified==FALSE,na.rm=TRUE)),by=site],"nco_availability_by_site.csv")
write_parquet(c,file.path(sensitive,"cohort_with_primary_nco.parquet"))
write_parquet(constructed$episodes,file.path(sensitive,"qualifying_prior_admission_episodes.parquet"))
weights<-readRDS(file.path(inp,"weight_diagnostics_patient_level.rds"));cfg<-weights$config
long<-as.data.frame(read_parquet(file.path(inp,"pneumonia_ccw_long.parquet")))
p<-g6_prepare_data(long,"site","hospital_block_id","m","event_m","censor_m","init_m",NULL,
 cfg$baseline,cfg$tv,cfg$horizon,cfg$lag_initial)
map<-unique(as.data.table(p$data)[,.(uid=as.character(.uid),hospital_block_id=as.character(.id))])
map<-merge(map,c[,.(hospital_block_id,prior_year_admission,nco_eligible)],by="hospital_block_id")
map[nco_eligible==FALSE,prior_year_admission:=NA_integer_]
if(!all(unlist(lapply(weights$site_components,`[[`,"uid"))%in%map$uid))stop("Saved clone IDs do not match prepared cohort")
results<-list()
collect_arms<-function(m,version,sites) {
 rbindlist(lapply(sites,function(s) {
  z<-weights$site_components[[s]];y<-map$prior_year_admission[match(z$uid,map$uid)]
  w1<-if(version=="raw")z$raw_W1[,m] else z$W1[,m]
  w0<-if(version=="raw")z$raw_W0[,m] else z$W0[,m]
  data.table(y=y,w1=w1*z$R1[,m],w0=w0*z$R0[,m])
 }))
}
site_names<-names(weights$site_components)
for(scope in c(site_names,"All sites")) {
 sites<-if(scope=="All sites")site_names else scope
 # Full-cohort observed-group comparison: descriptive, unweighted, not CCW.
 z<-c[site%in%sites];y<-z$prior_year_admission;y[!z$nco_eligible]<-NA
 results[[length(results)+1]]<-cbind(data.table(site=scope,assessment="Observed groups: full baseline cohort",version="unweighted",interval=0L,hours=NA_real_),
  nco_weighted_contrast(y,as.numeric(z$S<=cfg$tau),as.numeric(z$S>cfg$tau)))
 for(m in seq_len(cfg$tau))for(version in c("raw","clipped")) {
  z<-collect_arms(m,version,sites)
  results[[length(results)+1]]<-cbind(data.table(site=scope,assessment="CCW adherence/survivor risk set",version=version,
   interval=m,hours=(m-1)*6),nco_weighted_contrast(z$y,z$w1,z$w0))
 }
 # Saved landmark weights already contain identical confounders and setup.
 for(version in c("raw","clipped")) {
  z<-rbindlist(lapply(sites,function(s) {
   l<-weights$landmark[[s]];if(is.null(l))return(NULL)
   y<-map$prior_year_admission[match(l$uid,map$uid)]
   w<-if(version=="raw")l$raw_weights else l$weights
   data.table(y=y,w1=w*(l$A==1),w0=w*(l$A==0))
  }))
  if(nrow(z))results[[length(results)+1]]<-cbind(data.table(site=scope,assessment="Landmark survivor comparison",version=version,
   interval=cfg$tau,hours=cfg$tau*6),nco_weighted_contrast(z$y,z$w1,z$w0))
 }
}
res<-rbindlist(results,fill=TRUE);write(res,"primary_nco_weighted_contrasts.csv")
write(res[site=="All sites" & (interval==cfg$tau|interval==0)],"primary_nco_summary.csv")
plot<-ggplot(res[assessment=="CCW adherence/survivor risk set"],aes(hours,100*difference,color=version))+
 geom_hline(yintercept=0,linetype=3)+geom_errorbar(aes(ymin=100*lower,ymax=100*upper),width=.3)+
 geom_line()+geom_point()+facet_wrap(~site,ncol=2)+theme_minimal()+
 labs(title="Prior-year hospitalization NCO: weighted CCW risk-set differences",
 subtitle="Historical outcome; differences may reflect residual confounding and selection into these risk sets",
 x="Decision time, hours since time zero",y="NCO prevalence difference (percentage points)",color="Weight")
ggsave(file.path(out,"primary_nco_ccw_diagnostic.png"),plot,width=12,height=9,dpi=180)
writeLines(c("Primary NCO: any distinct within-system admission STARTING within 365 days before index block admission.",
 "Overlapping/<=6h-gap encounters are merged; the full expanded episode containing an index encounter is excluded.",
 "Only completed prior episodes count. Potential prior records with invalid dates make a negative NCO unknown, rather than zero.",
 "Future encounters never count. This is observed within-system utilization, not all hospitalization anywhere.",
 paste("Explicit system history-coverage start:",as.character(coverage)),
 "Without an explicit coverage start, zero means no observed admission; complete historical capture is NOT independently verified.",
 "With a coverage start, patients whose full 365-day window is not covered are excluded from diagnostic contrasts, not from 07.",
 "07 weights/confounders are reused unchanged. The NCO is not added to the weight model or treated as a day-28 survival event.",
 "CCW contrasts condition on evolving adherence/survival risk sets; zero is NOT automatically guaranteed even without confounding.",
 "Initial identical clones have identical baseline NCOs by construction; this equality is not evidence against confounding.",
 "Landmark comparisons also select survivors. Observed full-cohort initiation groups are descriptive and can involve treatment-timing selection.",
 "CIs combine the two clones at the original-patient level, treating fitted weights and coverage/selection as fixed; no bootstrap/refitting.",
 "These are bias/selection diagnostics, not a formal identified causal NCO effect, a definitive bias test, or calibration of primary estimates.",
 "Primary-effect estimates and all 07/08 files remain unchanged."),file.path(out,"METHOD_NOTES.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"));print(res[site=="All sites"])
