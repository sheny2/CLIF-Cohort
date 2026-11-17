#!/usr/bin/env Rscript
# Optional argument: existing SYNTHETIC test_pipeline.R result directory.
suppressPackageStartupMessages({library(data.table);library(arrow)})
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root<-dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
for(h in c("hospital_block_helpers.R","nco_history_helpers.R","nco_panel_helpers.R","nco_calibration_helpers.R"))source(file.path(root,h))
t0<-as.POSIXct("2024-06-01",tz="UTC")
c<-data.table(patient_id=LETTERS[1:5],hospital_block_id=paste0("i",1:5),admission_dttm=t0)
cw<-data.table(patient_id=c$patient_id,hospital_block_id=c$hospital_block_id,hospitalization_id=c$hospital_block_id)
h<-data.table(patient_id=c$patient_id,hospitalization_id=c$hospital_block_id,admission_dttm=t0,
  discharge_dttm=t0+86400,admission_type_category="ed",discharge_category="home")
h<-rbind(h,data.table(patient_id=c("A","A","B","C","D"),hospitalization_id=paste0("h",1:5),
  admission_dttm=t0-c(100,200,50,30,20)*86400,discharge_dttm=t0-c(90,198,49,29,19)*86400,
  admission_type_category=c("ed","elective",NA,"elective","ed"),discharge_category=c("snf","home","rehabilitation",NA,"home")))
x<-build_historical_nco_panel(c,cw,h)$cohort
stopifnot(x[patient_id=="A",nco_recurrent_admission]==1L,x[patient_id=="A",nco_prolonged]==1L,
  is.na(x[patient_id=="B",nco_ed]),is.na(x[patient_id=="C",nco_snf_rehab]),x[patient_id=="E",nco_any_admission]==0L)
y<-build_historical_nco_panel(c,cw,h[,!c("admission_type_category","discharge_category")])$cohort
stopifnot(all(is.na(y$nco_ed)),all(is.na(y$nco_snf_rehab)))
y<-build_historical_nco_panel(c,cw,h,as.Date("2024-01-01"))$cohort
stopifnot(all(is.na(y$nco_any_admission)))
z<-historical_log_pr(c(0,1,0,1),rep(1,4),rep(1,4))
stopifnot(z$log_pr==0,z$se_log_pr==0)
z<-historical_log_pr(c(0,1,0,1),c(0,1,0,1),c(1,0,1,0))
stopifnot(is.na(z$log_pr)) # no continuity correction when comparator prevalence is zero
set.seed(123)
est<-c(rnorm(25,-.15,.05),rnorm(8,.6,.07))
fit<-fit_mnn_calibration(est,rep(.03,length(est)))
stopifnot(fit$converged,abs(fit$mu_valid+.15)<.1,fit$pi_valid>.5)
cal<-calibrate_primary(log(1.2),.05,fit,"MNN")
stopifnot(abs(cal$estimate_logRR-(log(1.2)-fit$mu_valid))<1e-12,
  abs(cal$se_logRR-sqrt(.05^2+fit$sigma_valid^2))<1e-12)
cat("PASS: panel definitions, unknown fields, coverage, paired clones, zero prevalence, majority mixture and calibration arithmetic.\n")
a<-commandArgs(TRUE)
if(length(a)) {
  run<-normalizePath(a[1]);work<-tempfile("synthetic_nco_panel_");dir.create(work)
  c<-as.data.table(read_parquet(file.path(run,"sensitive","analysis","pneumonia_ccw_cohort.parquet")))
  cw<-as.data.table(read_parquet(file.path(run,"sensitive","analysis","sampled_hospital_block_crosswalk.parquet")))
  raw<-as.data.table(read_parquet(file.path(dirname(run),"data","clif_hospitalization.parquet")))
  h<-raw[patient_id%in%c$patient_id,c("patient_id","hospitalization_id","admission_dttm","discharge_dttm","admission_type_category","discharge_category"),with=FALSE]
  set.seed(923)
  prior<-rbindlist(lapply(seq_len(nrow(c)),function(i) {
    n<-rpois(1,1.4);if(!n)return(NULL)
    start<-as.POSIXct(c$admission_dttm[i],tz="UTC")-(40+seq_len(n)*30)*86400
    data.table(patient_id=c$patient_id[i],hospitalization_id=paste0("prior_",i,"_",seq_len(n)),
      admission_dttm=start,discharge_dttm=start+sample(c(2,9),n,TRUE)*86400,
      admission_type_category=sample(c("ed","elective"),n,TRUE),
      discharge_category=sample(c("home","skilled_nursing_facility","rehabilitation"),n,TRUE))
  }))
  path<-file.path(work,"expanded.parquet");write_parquet(rbind(h,prior),path)
  reference<-file.path(run,"aggregate","analysis","five_method_results.csv");before<-tools::md5sum(reference)
  for(script in c("09_prior_hospitalization_nco_diagnostic.R","10_historical_nco_calibration.R")) {
    log<-file.path(work,paste0(script,".log"))
    status<-system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(root,script)),shQuote(run),shQuote(path),"2020-01-01"),stdout=log,stderr=log)
    if(status!=0){cat(readLines(log),sep="\n");stop(script," integration failed")}
  }
  out<-file.path(run,"aggregate","nco_calibration")
  fits<-fread(file.path(out,"nn_mnn_fit_diagnostics.csv"))
  z<-fread(file.path(out,"exploratory_calibrated_ccw_rr.csv"))
  stopifnot(identical(before,tools::md5sum(reference)),all(c("NN","MNN")%in%fits$model),
    all(z$calibration_supported==FALSE),any(is.finite(z$rr)),file.exists(file.path(out,"COMPLETE.txt")))
  cat("PASS: 09--10 integration; original clinical estimates unchanged. Outputs:",out,"\n")
}
