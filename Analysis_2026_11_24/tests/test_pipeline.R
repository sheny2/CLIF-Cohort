#!/usr/bin/env Rscript
# Synthetic integration test: creates no real patient data and never touches inputs.
# Usage: HOSPITAL_BLOCK_HOURS=3 Rscript tests/test_pipeline.R (also test 6 and 12).
suppressPackageStartupMessages({library(data.table);library(arrow)})
file <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(normalizePath(gsub("~+~"," ",file,fixed=TRUE))))
work <- file.path("/tmp",paste0("clif_nov24_",Sys.getpid()));dir.create(work)
data <- file.path(work,"data");dir.create(data)
source(file.path(root,"variable_rules.R"))
stopifnot(identical(screen_measurements(c(0,29,30,37,45,46),"temperature"),c(NA_real_,NA_real_,30,37,45,NA_real_)))
stopifnot(identical(screen_measurements(c(.5,50,2),"fio2"),c(.5,.5,NA_real_)))
stopifnot(identical(map_support(c("imv","nippv","cpap","high_flow_nc","room_air","face_mask","other")),
 c("imv","nippv_cpap","nippv_cpap","hfno","room_air","conventional_oxygen",NA_character_)))
set.seed(152026)
n <- 350L;t0 <- as.POSIXct("2024-01-15 12:00:00",tz="UTC")+seq_len(n)*86400
patient <- sprintf("p%04d",seq_len(n));hospital <- sprintf("h%04d",seq_len(n))
dead <- runif(n)<.23; deathhours <- sample(c(5,14,30,60,120,250,450,650),n,TRUE)
dischargehours <- sample(5:16,n,TRUE)*24
dischargehours[dead & seq_len(n)%%2==0] <- pmin(dischargehours[dead&seq_len(n)%%2==0],deathhours[dead&seq_len(n)%%2==0])
dischargehours[c(7,8)] <- 12;dead[c(7,8)] <- FALSE
write <- function(x,f)write_parquet(x,file.path(data,paste0("clif_",f,".parquet")))
write(data.table(patient_id=patient,sex_category=sample(c("male","female"),n,TRUE),race_category="White",ethnicity_category="Not Hispanic",
 death_dttm=as.POSIXct(ifelse(dead,as.numeric(t0)+deathhours*3600,NA_real_),origin="1970-01-01",tz="UTC")),"patient")
sites <- rep(c("JHH","BMC","SMH","HCG","SUB"),length.out=n)
write(data.table(patient_id=patient,hospitalization_id=hospital,admission_dttm=t0-4*3600,
 discharge_dttm=t0+dischargehours*3600,age_at_admission=sample(25:90,n,TRUE),admission_type_category="ed",
 discharge_category=ifelse(dead&deathhours==dischargehours,"expired",ifelse(seq_len(n)%%17==0,"hospice","home")),hospital_id=sites),"hospitalization")
dx <- data.table(hospitalization_id=hospital,diagnosis_code="J18.9",diagnosis_code_format="ICD10CM",diagnosis_primary=TRUE,poa_present=TRUE)
chronic <- data.table(hospitalization_id=hospital,diagnosis_code=sample(c("J44","E1165","D80","C78","C90","I50","N18","I10"),n,TRUE),
 diagnosis_code_format="ICD10CM",diagnosis_primary=FALSE,poa_present=TRUE)
write(rbind(dx,chronic),"hospital_diagnosis")
write(data.table(patient_id=patient,start_dttm=t0-3600,code_status_category="full"),"code_status")
write(data.table(hospitalization_id=hospital,admin_dttm=t0-3600,med_category="norepinephrine",med_dose=0,
 med_dose_unit="mcg/kg/min",mar_action_category="stopped"),"medication_admin_continuous")
write(data.table(hospitalization_id=hospital,hospital_id=sites,in_dttm=t0-4*3600,out_dttm=t0+dischargehours*3600,
 location_category=sample(c("icu","ed","ward"),n,TRUE)),"adt")
vitals <- labs <- respiratory <- assessment <- list()
for(i in seq_len(n)) {
 hours <- seq(0,dischargehours[i],6)
 support <- c("imv","nippv","high_flow_nc")[(i%%3)+1]
 if(i==5)support<-"high_flow_nc"
 device <- ifelse(hours<48,support,ifelse(hours<96,"nasal_cannula","room_air"))
 if(i==5)device[]<-"high_flow_nc"
 respiratory[[i]] <- data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i]+hours*3600,
  device_category=device,fio2_set=ifelse(hours<48,.6,.3),lpm_set=40,tracheostomy=FALSE)
 sp <- if(i==5)98 else if(i==6)99 else 95
 temp_hours <- c(-1,hours)
 if(i==2)temp_hours<-c(2,hours[hours>0])
 if(i==3)temp_hours<-numeric()
 temp <- data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i]+temp_hours*3600,vital_category="temp_c",
 vital_value=36.8+rnorm(length(temp_hours),0,.7))
 if(i==1)temp<-rbind(temp,data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i],vital_category="temp_c",vital_value=0))
 vitals[[i]] <- rbind(temp,data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i]+hours*3600,vital_category="spo2",vital_value=sp),
  data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i]+hours*3600,vital_category="map",vital_value=70+runif(length(hours),0,30)),
  data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i]-3600,vital_category="weight_kg",vital_value=80))
 lh <- c(-1,if(dischargehours[i]>=24)seq(24,dischargehours[i],24) else numeric()); cats<-c("wbc","glucose_serum","albumin","platelet_count","bilirubin_total","creatinine","lactate","ph_arterial","pco2_arterial")
 lab <- CJ(h=lh,lab_category=cats)
 lab[,lab_value_numeric:=fcase(lab_category=="wbc",10+runif(.N,0,10),lab_category=="glucose_serum",100+runif(.N,0,150),
  lab_category=="albumin",2+runif(.N,0,2),lab_category=="platelet_count",100+runif(.N,0,200),lab_category=="bilirubin_total",runif(.N,.3,2),
  lab_category=="creatinine",runif(.N,.5,3),lab_category=="lactate",1,lab_category=="ph_arterial",7.3,default=45)]
 lab[,`:=`(hospitalization_id=hospital[i],lab_collect_dttm=t0[i]+h*3600)];lab[,h:=NULL]
 if(i==4)lab[lab_category=="glucose_serum",lab_value_numeric:=2007]
 labs[[i]]<-lab
 assessment[[i]]<-data.table(hospitalization_id=hospital[i],recorded_dttm=t0[i]+hours*3600,assessment_category="gcs_total",numerical_value=sample(10:15,length(hours),TRUE))
}
write(rbindlist(vitals),"vitals");write(rbindlist(labs),"labs");write(rbindlist(respiratory),"respiratory_support");write(rbindlist(assessment),"patient_assessments")
group <- sample(c("early","late","never"),n,TRUE)
init <- ifelse(group=="early",runif(n,.5,23),ifelse(group=="late",runif(n,30,150),800))
med <- data.table(hospitalization_id=hospital,admin_dttm=t0+init*3600,med_group="steroid",med_category="hydrocortisone",
 med_route_category="iv",med_dose=50,med_dose_unit="mg",mar_action_category="given")
med[7,mar_action_category:="not_given"]
med[8,admin_dttm:=t0[8]+6*3600]
write(med,"medication_admin_intermittent")
Sys.setenv(TREATMENT_GRACE_HOURS=24,CLIF_DATA_CUTOFF_DATE="2026-10-31",HOSPITAL_BLOCK_GAP_HOURS=6)
source(file.path(root,"analysis_preparation.R"))
timing_config <- analysis_timing(24)
step_hours <- timing_config$step_hours
M <- timing_config$M; tau <- timing_config$tau
run <- function(script,args) {
 log <- file.path(work,paste0(substr(script,1,2),".log"))
 status<-system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(root,script)),vapply(args,shQuote,character(1))),stdout=log,stderr=log)
 if(status!=0){cat(readLines(log),sep="\n");stop("Failed: ",script,"; fixture: ",work)}
}
cohort<-file.path(work,"cohort");profile<-file.path(work,"profile");candidate<-file.path(work,"candidate");treatment<-file.path(work,"treatment")
run("01_build_adult_pneumonia_cohort.R",c(data,cohort))
c<-as.data.table(read_parquet(file.path(cohort,"adult_pneumonia_cohort.parquet")))
stopifnot(!"p0005"%in%c$patient_id,nrow(c)==n-1)
design01 <- fread(file.path(cohort,"adult_pneumonia_cohort_design.csv"))
stopifnot(as.numeric(design01[parameter=="hfno_spo2_max",value])==97)
run("02_profile_glucocorticoid_sources.R",c(data,profile))
run("03_extract_candidate_steroid_administrations.R",c(data,cohort,candidate))
run("04_define_steroid_treatment.R",c(candidate,cohort,treatment,"24"))
eda<-file.path(work,"eda");eda_sensitive<-file.path(work,"eda_sensitive")
run("05_cohort_summary_eda.R",c(data,cohort,treatment,eda,"1","2026-10-31",eda_sensitive))
result<-file.path(work,"result")
run("07_prepare_and_run_five_methods.R",c(data,cohort,treatment,result,"2026-10-31"))
long<-as.data.table(read_parquet(file.path(result,"sensitive","analysis","pneumonia_ccw_long.parquet")))
primary <- as.data.table(read_parquet(file.path(result,"sensitive","analysis","pneumonia_ccw_cohort.parquet")))
stopifnot(primary[post_discharge_death_ignored & !grepl("hospice|expired",tolower(discharge_category)),all(delta==0L)],
          primary[event_type=="hospice_discharge",sum(delta==1L)]>0L)
stopifnot(setequal(list.dirs(result,full.names=FALSE,recursive=FALSE),c("aggregate","sensitive")))
overall <- fread(file.path(eda,"cohort_overall_summary.csv"))
stopifnot(as.integer(overall[parameter=="outcome_within_28d_n",value])==sum(primary$delta))
timing <- fread(file.path(eda,"outcome_and_treatment_six_hour_counts.csv"))
stopifnot(sum(timing[event_type=="In-hospital death or hospice discharge",n_reportable])==sum(primary$delta))
stopifnot(!any(grepl("^(ph|paco2|support_score)(_|$)",names(long))))
stopifnot(all(is.na(long$sf_ratio[long$spo2>97])))
estimates <- fread(file.path(result,"aggregate","analysis","five_method_results.csv"))
stopifnot(nrow(estimates)==5,all(is.finite(as.matrix(estimates[,.(risk_difference,risk_ratio,rmst_difference_days)]))))
fill<-as.data.table(read_parquet(file.path(result,"sensitive","analysis","adjacent_fill_rows_temperature.parquet")))
b2<-c[patient_id=="p0002",hospital_block_id]
stopifnot(fill[hospital_block_id==b2&m==1,future],fill[hospital_block_id==b2&m==1,distance_hours]==2)
b1<-c[patient_id=="p0001",hospital_block_id]
stopifnot(long[hospital_block_id==b1&m==1,temperature]>=30)
spec<-fread(file.path(result,"aggregate","analysis","ipcw_model_specification.csv"))
stopifnot(!any(grepl("support_score|paco2|^ph_|sf_ratio_model|sf_ratio_missing",spec$term)),
 any(grepl("^sf_",spec$term)),all(is.na(long$sf_ratio[long$spo2>97])))
states<-fread(file.path(eda,"respiratory_support_proportions_full_cohort.csv"))
stopifnot(all(abs(states[,sum(proportion),by=.(site,hours)]$V1-1)<1e-10))
diag<-readRDS(file.path(result,"sensitive","analysis","weight_diagnostics_patient_level.rds"))
stopifnot(primary[patient_id=="p0007",early_alive_discharge_untreated],
  !primary[patient_id=="p0008",early_alive_discharge_untreated])
component <- diag$site_components[[primary[patient_id=="p0007",site]]]
source(file.path(root,"06_generalized_five_methods.R"))
prepared <- g6_prepare_data(as.data.frame(long),"site","hospital_block_id","m","event_m","censor_m","init_m",NULL,
  diag$config$baseline,diag$config$tv,diag$config$horizon,diag$config$lag_initial)$data
uid <- as.character(unique(prepared$.uid[prepared$patient_id=="p0007"]))
i <- match(uid,component$uid)
stopifnot(all(component$R1[i,]),all(component$R0[i,]),
  all(component$raw_W1[i,component$adherence_end_m[i]:M]==component$raw_W1[i,component$adherence_end_m[i]]),
  all(component$raw_W0[i,component$adherence_end_m[i]:M]==component$raw_W0[i,component$adherence_end_m[i]]))
for(component in diag$site_components) {
 late<-which(component$S>diag$config$tau&component$S<=diag$config$horizon)
 for(i in late)stopifnot(component$R0[i,component$S[i]])
 # With no natural censoring, adherence weights freeze after grace.
 stopifnot(all(abs(component$raw_W0[,tau+1L]-component$raw_W0[,tau])<1e-10))
}
# Configurability check: same source/qualifying doses, fresh 96-hour analysis.
Sys.setenv(TREATMENT_GRACE_HOURS=96)
result96 <- file.path(work,"result96")
run("07_prepare_and_run_five_methods.R",c(data,cohort,treatment,result96,"2026-10-31"))
diag96 <- readRDS(file.path(result96,"sensitive","analysis","weight_diagnostics_patient_level.rds"))
tau96 <- as.integer(96 / step_hours)
stopifnot(diag96$config$tau==tau96)
late96 <- 0L
for(component in diag96$site_components) {
 late <- which(component$S>tau96 & component$S<=M)
 late96 <- late96+length(late)
 for(i in late)stopifnot(component$R0[i,component$S[i]])
 stopifnot(all(abs(component$raw_W0[,tau96+1L]-component$raw_W0[,tau96])<1e-10))
}
stopifnot(late96>0)
# Independently check interval dimensions, boundaries, and saved design.
for(run_root in c(result,result96)) {
 design <- readRDS(file.path(run_root,"sensitive","analysis","analysis_design.rds"))
 z <- as.data.table(read_parquet(file.path(run_root,"sensitive","analysis","pneumonia_ccw_long.parquet")))
 u <- as.data.table(read_parquet(file.path(run_root,"sensitive","analysis","pneumonia_ccw_cohort.parquet")))
 stopifnot(design$step_hours==step_hours,design$horizon_intervals==M,
   design$tau_intervals==design$grace_hours/step_hours,
   all(as.numeric(difftime(z$decision_dttm,z$index_dttm,units="hours"))==(z$m-1L)*step_hours),
   all(u$S==ifelse(is.finite(u$treatment_hours),pmin(M+1L,pmax(1L,ceiling(u$treatment_hours/step_hours))),M+1L)),
   all(u$T_event[u$delta==1L]==pmax(1L,ceiling(u$event_hours[u$delta==1L]/step_hours))))
 run("tests/test_local_meta.R",run_root)
}
cat("PASS: full 01--08 pipeline, composite endpoint agreement across EDA/methods, two-folder outputs, 24h/96h strategies, and weight diagnostics.\nFixture/output: ",work,"\n",sep="")
