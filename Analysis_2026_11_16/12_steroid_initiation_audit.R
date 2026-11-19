#!/usr/bin/env Rscript
# Rscript 12_steroid_initiation_audit.R RUN_ROOT [STEROID_TREATMENT_DIR]
suppressPackageStartupMessages({library(arrow);library(data.table);library(ggplot2)})
setDTthreads(1L)
a <- commandArgs(TRUE)
if(!length(a))stop("Supply the completed 07 RUN_ROOT")
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root <- dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE)))
source(file.path(root,"site_helpers.R"));source(file.path(root,"steroid_start_helpers.R"))
run <- normalizePath(a[1]);input <- file.path(run,"sensitive","analysis")
treatment <- if(length(a)>=2L)a[2L] else file.path(run,"sensitive","steroid_treatment")
if(!file.exists(file.path(run,"aggregate","analysis","COMPLETE.txt")))stop("Complete 07 first")
paths <- file.path(treatment,c("qualifying_steroid_administrations.parquet","steroid_administrations_with_treatment_eligibility.parquet"))
if(any(!file.exists(paths)))stop("Supply the 04 output directory used for this 07 run; missing: ",paste(paths[!file.exists(paths)],collapse=", "))
c <- as.data.table(read_parquet(file.path(input,"pneumonia_ccw_cohort.parquet")))
assert_mapped_sites(c$site)
design <- readRDS(file.path(input,"analysis_design.rds"))
q <- as.data.table(read_parquet(paths[1]));doses <- as.data.table(read_parquet(paths[2]))
required <- c("hospital_block_id","admin_dttm","med_category_clean","med_route_clean","dose_mg","mar_action_category")
if(!all(required%in%names(q))||!all(c(required,"route_eligible_for_drug","minimum_dose_mg")%in%names(doses)))
  stop("Rerun current 03/04; drug/dose/action eligibility fields are required")
audit <- build_steroid_start_audit(c,q,doses,design);c <- audit$cohort
out <- file.path(run,"aggregate","steroid_initiation_audit");sensitive <- file.path(run,"sensitive","steroid_initiation_audit")
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
unlink(file.path(out,"COMPLETE.txt"))
unlink(file.path(out,"interval1_doses_around_zero.png"))
write <- function(z,name) {
  if(any(c("patient_id","hospital_block_id","hospitalization_id")%in%names(z)))stop("Identifier in aggregate output")
  fwrite(z,file.path(out,name),na="")
}
write_parquet(c,file.path(sensitive,"patient_steroid_start_audit.parquet"))
write_parquet(audit$qualifying,file.path(sensitive,"qualifying_dose_timing.parquet"))
scopes <- rbind(copy(c)[,scope:="All sites"],copy(c)[,scope:=site])
summarize <- function(z)z[,.(patients=.N,prior_qualifying_24h=sum(prior_qualifying_24h),
  prior_qualifying_24h_percent=100*mean(prior_qualifying_24h),
  prior_qualifying_6h=sum(prior_qualifying_6h),prior_qualifying_12h=sum(prior_qualifying_12h),
  any_qualifying_before_zero=sum(any_qualifying_before_zero),any_given_systemic_24h=sum(any_given_systemic_24h),
  same_drug_prior_24h=sum(same_drug_prior_24h),full_24h_in_selected_block=sum(full_24h_in_selected_block),
  prior_dose_before_block_admission=sum(prior_dose_before_block_admission),
  exactly_time_zero=sum(first_post_hours==0,na.rm=TRUE),source_mismatches=sum(!saved_start_matches_source)),by=.(scope,start_group)]
summary <- summarize(scopes)
write(summary,"pre_zero_exposure_by_start_group.csv")
write(summary[start_group=="Interval 1 (0–6 h)"],"interval1_prior_exposure_summary.csv")
status <- scopes[,.(patients=.N),by=.(scope,start_group,prior_status)]
status[,percent:=100*patients/sum(patients),by=.(scope,start_group)]
write(status,"prior_exposure_and_lookback.csv")
write(scopes[,.(patients=.N,mismatches=sum(!saved_start_matches_source)),by=scope],"saved_treatment_source_agreement.csv")
if(any(!c$saved_start_matches_source))warning("Saved 07 treatment timing differs from the supplied 04 source. Inspect the sensitive audit and rerun matching inputs before interpretation.")
interval1 <- merge(audit$prior,c[interval1_starter==TRUE,.(hospital_block_id,site)],by="hospital_block_id")
drugs <- rbind(copy(interval1)[,scope:="All sites"],copy(interval1)[,scope:=site])
drug_counts <- drugs[,.(patients=uniqueN(hospital_block_id),administrations=.N),by=.(scope,med_category_clean,med_route_clean)]
write(drug_counts,"interval1_prior_drugs_and_routes.csv")
timing <- merge(audit$qualifying[hours_from_zero>= -24&hours_from_zero<=6],c[,.(hospital_block_id,site,interval1_starter)],by="hospital_block_id")
timing[,hour_bin:=pmin(5L,as.integer(floor(hours_from_zero)))]
timing_scopes <- rbind(copy(timing)[,scope:="All sites"],copy(timing)[,scope:=site])
hist <- timing_scopes[,.(administrations=.N,patients=uniqueN(hospital_block_id)),by=.(scope,interval1_starter,hour_bin)]
write(hist,"qualifying_doses_minus24_to_plus6h.csv")
long <- as.data.table(read_parquet(file.path(input,"pneumonia_ccw_long.parquet")))
long <- merge(long[,intersect(c("hospital_block_id","site","m","S","init_m","adherence_end_m"),names(long)),with=FALSE],
  c[,.(hospital_block_id,prior_qualifying_24h)],by="hospital_block_id")
decision <- long[m<=design$tau_intervals&m<=S]
if("adherence_end_m"%in%names(decision))decision <- decision[m<=adherence_end_m]
risk <- rbind(copy(decision)[,scope:="All sites"],copy(decision)[,scope:=site])
hazard <- risk[,.(eligible_decisions=.N,starts=sum(init_m),starts_with_prior_qualifying_24h=sum(init_m==1L&prior_qualifying_24h)),by=.(scope,m)]
hazard[,`:=`(interval_start_hours=(m-1L)*design$step_hours,interval_end_hours=m*design$step_hours,
  observed_start_hazard=starts/eligible_decisions,
  percent_starts_with_prior24=ifelse(starts>0,100*starts_with_prior_qualifying_24h/starts,NA_real_))]
write(hazard,"observed_start_hazards_and_prior_exposure.csv")
compare <- scopes[,.(patients=.N,early_initiators=sum(S<=design$tau_intervals),
  early_percent=100*mean(S<=design$tau_intervals)),by=.(scope,prior_qualifying_24h,full_24h_in_selected_block)]
write(compare,"descriptive_prior_exposure_strata.csv")
saveplot <- function(p,name)ggsave(file.path(out,name),p,width=12,height=7,dpi=180)
saveplot(ggplot(status[start_group=="Interval 1 (0–6 h)"],aes(scope,percent,fill=prior_status))+
  geom_col()+theme_minimal()+labs(title="Interval-1 starters: pre-time-zero qualifying steroid exposure",x=NULL,y="Starters (%)",fill=NULL)+
  theme(legend.position="bottom"),"interval1_prior_exposure.png")
if(any(hist$interval1_starter))saveplot(ggplot(hist[interval1_starter==TRUE],aes(hour_bin,administrations))+
  geom_col(fill="#2878A5")+geom_vline(xintercept=0,linetype=3)+facet_wrap(~scope,scales="free_y")+theme_minimal()+
  labs(title="Qualifying doses around time zero among interval-1 starters",subtitle="Counts administrations, including repeat doses; zero belongs to post-time-zero",
    x="Hours from time zero (one-hour bins; final bin includes hour 6)",y="Administrations"),"interval1_doses_around_zero.png")
saveplot(ggplot(hazard,aes(interval_start_hours,100*observed_start_hazard,color=scope))+
  geom_line()+geom_point()+theme_minimal()+labs(title="Observed initiation hazard during grace",x="Interval start, hours",y="Starts / eligible decisions (%)",color=NULL),"observed_start_hazards.png")
writeLines(c("Interval 1 is the saved 07 S=1 group: first qualifying post-zero dose in [0,6] hours, not necessarily a new treatment course.",
  "Primary lookback is [-24,0) hours. Additional counts use [-12,0), [-6,0), and all observed pre-zero doses linked to the selected block.",
  "Qualifying doses reuse 04 drug/dose/route thresholds and require given. An additional any-given-systemic measure includes positive subthreshold doses of the five eligible drugs/routes.",
  "No doses from after discharge, primary event, or administrative cutoff count. Records linked to selected encounters but timestamped before block admission are retained and flagged.",
  "Nominal 24h within block means time zero minus admission >=24h, not verified complete MAR coverage. A shorter interval and no observed dose do not establish steroid-naive status.",
  "No observed prior dose cannot exclude outpatient, prehospital, outside-system or unlinked ED treatment. ED-origin admission or initial ED location does not identify where a dose was given.",
  "Prior qualifying exposure and same-drug repetition support possible continuation, not proof of a chronic regimen or proof of the indication.",
  "Drug/route counts can count a patient in multiple categories. Timing histograms count administrations, not first starts; exact hour 6 is included in the last bin.",
  "Start hazards use saved untreated decision rows during grace, stopping at alive discharge when the saved adherence_end_m is available. They are unweighted observed hazards, not fitted probabilities.",
  "Saved 07 treatment/source agreement is audited. Source mismatch requires reconciling the selected run and 04 inputs.",
  "Observed groups and hypothetical exclusion strata are descriptive; no eligibility restriction, treatment reassignment, or effect-model refit is applied.",
  "01--11 outputs remain unchanged. Aggregate cells are unsuppressed; individual records stay in sensitive/."),file.path(out,"METHOD_NOTES.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"))
print(summary[start_group=="Interval 1 (0–6 h)"]);cat("\nResults:",out,"\n")
