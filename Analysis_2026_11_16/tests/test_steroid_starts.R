#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow);library(data.table)})
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root <- dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"steroid_start_helpers.R"))
t0 <- as.POSIXct("2024-01-02",tz="UTC")
c <- data.table(patient_id=LETTERS[1:7],hospital_block_id=paste0("b",1:7),site=rep(c("JHH","HCG"),length.out=7),
  index_dttm=t0,admission_dttm=t0-c(48,2,48,48,48,48,48)*3600,discharge_dttm=t0+100*3600,
  event_dttm=as.POSIXct(NA,tz="UTC"),S=c(1L,1L,2L,1L,1L,1L,113L))
c[patient_id=="D",event_dttm:=t0+5*3600]
design <- list(data_cutoff="2024-01-03",horizon_days=28L,horizon_intervals=112L,tau_intervals=4L,step_hours=6L)
doses <- data.table(hospital_block_id=c("b1","b1","b1","b2","b2","b3","b4","b4","b4","b5","b5","b6","b6","b7"),
  admin_dttm=t0+c(-24,-1,0,-25,3,7,-2,4,10,-2,6,-12,1,80)*3600,
  med_category_clean=c(rep("hydrocortisone",11),"methylprednisolone","hydrocortisone","hydrocortisone"),
  med_route_clean="iv",dose_mg=c(rep(50,9),25,50,32,50,50),
  mar_action_category=c(rep("given",6),"not_given",rep("given",7)),
  route_eligible_for_drug=TRUE,minimum_dose_mg=c(rep(50,11),32,50,50))
q <- doses[dose_mg>=minimum_dose_mg]
x <- build_steroid_start_audit(c,q,doses,design)$cohort
stopifnot(x[patient_id=="A",prior_qualifying_doses_24h]==2L,
  x[patient_id=="A",same_drug_prior_24h],x[patient_id=="A",first_post_hours]==0,
  x[patient_id=="A",prior_to_post_gap_hours]==1,
  !x[patient_id=="B",prior_qualifying_24h],x[patient_id=="B",any_qualifying_before_zero],
  !x[patient_id=="B",full_24h_in_selected_block],
  !x[patient_id=="D",prior_qualifying_24h],
  !x[patient_id=="E",prior_qualifying_24h],x[patient_id=="E",any_given_systemic_24h],
  x[patient_id=="E",rebuilt_S]==1L,x[patient_id=="F",prior_qualifying_24h],
  !x[patient_id=="F",same_drug_prior_24h],all(x$saved_start_matches_source))
bad <- copy(c);bad[patient_id=="G",S:=2L]
stopifnot(!build_steroid_start_audit(bad,q,doses,design)$cohort[patient_id=="G",saved_start_matches_source])
a <- commandArgs(TRUE);work <- if(length(a))a[1L] else tempfile("steroid_start_audit_")
inp <- file.path(work,"sensitive","analysis");dir.create(inp,recursive=TRUE)
tx <- file.path(work,"sensitive","steroid_treatment");dir.create(tx,recursive=TRUE)
agg <- file.path(work,"aggregate","analysis");dir.create(agg,recursive=TRUE)
writeLines("Synthetic",file.path(agg,"COMPLETE.txt"))
write_parquet(c,file.path(inp,"pneumonia_ccw_cohort.parquet"));saveRDS(design,file.path(inp,"analysis_design.rds"))
long <- c[,.(m=seq_len(if(patient_id=="D")1L else 112L)),by=.(hospital_block_id,site,S)]
long[,init_m:=as.integer(S==m)];write_parquet(long,file.path(inp,"pneumonia_ccw_long.parquet"))
write_parquet(q,file.path(tx,"qualifying_steroid_administrations.parquet"))
write_parquet(doses,file.path(tx,"steroid_administrations_with_treatment_eligibility.parquet"))
before <- tools::md5sum(file.path(agg,"COMPLETE.txt"));log <- file.path(work,"12.log")
status <- system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(root,"12_steroid_initiation_audit.R")),shQuote(work)),stdout=log,stderr=log)
if(status!=0){cat(readLines(log),sep="\n");stop("12 integration failed")}
out <- file.path(work,"aggregate","steroid_initiation_audit")
summary <- fread(file.path(out,"interval1_prior_exposure_summary.csv"))[scope=="All sites"]
stopifnot(summary$patients==5L,summary$prior_qualifying_24h==2L,
  summary$prior_qualifying_24h_percent==40,summary$source_mismatches==0L,
  file.exists(file.path(out,"COMPLETE.txt")),length(list.files(out,pattern="[.]png$"))==3L,
  identical(before,tools::md5sum(file.path(agg,"COMPLETE.txt"))))
cat("PASS: -24/0/6h boundaries, not_given, subthreshold exposure, same-drug checks, incomplete lookback, source mismatches, and script 12.\nOutputs:",out,"\n")
