#!/usr/bin/env Rscript
# Synthetic, nonclinical data only. Exercises every CLI stage in an isolated copy.
suppressPackageStartupMessages({library(arrow);library(data.table)})
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
folder <- dirname(dirname(normalizePath(self)))
tmp <- tempfile('proposal_e2e_',tmpdir='/tmp');dir.create(tmp);dir.create(file.path(tmp,'pipeline'));dir.create(file.path(tmp,'data'))
for(f in list.files(folder,pattern='[.]R$',full.names=TRUE)) file.copy(f,file.path(tmp,'pipeline',basename(f)))
set.seed(1941);n <- 240L;ids <- paste0('h',1:n);pids <- paste0('p',1:n)
t0 <- as.POSIXct('2025-01-01',tz='UTC')+seq_len(n)*86400
# Severity occurs two hours after admission. Survival/status remains independent of treatment in this fixture.
death <- ifelse(runif(n)<.23,runif(n,24,600),NA_real_)
disch <- t0+15*86400;disch[!is.na(death)] <- t0[!is.na(death)]+death[!is.na(death)]*3600
put <- function(x,name) write_parquet(x,file.path(tmp,'data',paste0('clif_',name,'.parquet')))
put(data.table(patient_id=pids,sex_category=sample(c('male','female'),n,TRUE),race_category='unknown',ethnicity_category='unknown',death_dttm=t0+death*3600),'patient')
put(data.table(hospitalization_id=ids,patient_id=pids,admission_dttm=t0,discharge_dttm=disch,
 age_at_admission=sample(30:85,n,TRUE),admission_type_category='ed',discharge_category=ifelse(is.na(death),'home','expired'),site_id=rep(c('JHH','SMH'),each=n/2)),'hospitalization')
put(data.table(hospitalization_id=ids,diagnosis_code='J18.9',diagnosis_code_format='ICD10CM',diagnosis_primary=1L,poa_present=1L),'hospital_diagnosis')
put(data.table(hospitalization_id=ids,in_dttm=t0,out_dttm=disch,location_category='icu',hospital_id=rep(c('JHH','SMH'),each=n/2)),'adt')
put(data.table(patient_id=pids,start_dttm=t0,code_status_category='full'),'code_status')
grid <- CJ(i=1:n,h=seq(0,360,by=12));grid[,t:=t0[i]+h*3600];grid <- grid[t<=disch[i]]
resp <- grid[,.(hospitalization_id=ids[i],recorded_dttm=t+2*3600,device_category=rep(c('imv','nippv','high_flow_nc'),length.out=.N),fio2_set=.6,lpm_set=40,tracheostomy=0L,peep_set=6)]
put(resp,'respiratory_support')
v <- rbindlist(lapply(c('temperature','spo2','map','weight_kg'),function(v) grid[,.(hospitalization_id=ids[i],recorded_dttm=t,vital_category=v,vital_value=switch(v,temperature=rnorm(.N,37,1),spo2=runif(.N,90,96),map=runif(.N,75,95),weight_kg=runif(.N,60,100)))]))
put(v,'vitals')
l <- rbindlist(lapply(c('po2_arterial','wbc','glucose_serum','albumin','pco2_arterial','platelet_count','bilirubin_total','creatinine','lactate'),function(v) grid[,.(hospitalization_id=ids[i],lab_collect_dttm=t,lab_category=v,lab_value_numeric=switch(v,po2_arterial=runif(.N,70,120),ph_arterial=runif(.N,7.25,7.45),wbc=runif(.N,5,20),glucose_serum=runif(.N,80,200),albumin=runif(.N,2,4),pco2_arterial=runif(.N,30,50),platelet_count=runif(.N,100,300),bilirubin_total=runif(.N,.3,1.5),creatinine=runif(.N,.5,2),lactate=runif(.N,.5,2)))]))
put(l,'labs')
put(grid[,.(hospitalization_id=ids[i],recorded_dttm=t,assessment_category='gcs_total',numerical_value=sample(12:15,.N,TRUE))],'patient_assessments')
treated <- which(runif(n)<.5);infusion <- treated[seq_along(treated)%%2==0];bolus <- setdiff(treated,infusion)
base <- data.table(hospitalization_id=ids,admin_dttm=t0,med_category='norepinephrine',med_dose=0,med_dose_unit='mcg/kg/min',mar_action_category='stopped',med_route_category='iv')
cont <- rbind(base,data.table(hospitalization_id=ids[infusion],admin_dttm=t0[infusion]+runif(length(infusion),4,20)*3600,med_category='hydrocortisone',med_dose=200/24,med_dose_unit='mg/hr',mar_action_category='started',med_route_category='iv'))
# Prior steroid users test baseline exclusion; incompatible regimens test censoring.
inter <- rbind(data.table(hospitalization_id=ids[bolus],admin_dttm=t0[bolus]+runif(length(bolus),4,20)*3600,med_category='hydrocortisone',med_dose=50,med_dose_unit='mg',mar_action_category='given',med_route_category='iv'),
 data.table(hospitalization_id=ids[1:5],admin_dttm=t0[1:5]+3600,med_category='prednisone',med_dose=5,med_dose_unit='mg',mar_action_category='given',med_route_category='enteral'),
 data.table(hospitalization_id=ids[6:10],admin_dttm=t0[6:10]+4*3600,med_category='dexamethasone',med_dose=6,med_dose_unit='mg',mar_action_category='given',med_route_category='iv'))
inter[,med_group:='steroid']
put(cont,'medication_admin_continuous');put(inter,'medication_admin_intermittent')
Sys.setenv(PROPOSAL_RUN_ID='synthetic',TREATMENT_GRACE_HOURS='24',PROPOSAL_MAX_MISSING='0.70',SUMMARY_SMALL_CELL_THRESHOLD='10')
status <- system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(tmp,'pipeline','run_pipeline.R')),shQuote(file.path(tmp,'data'))))
if(status!=0) stop('Synthetic pipeline failed. Inspect ',file.path(tmp,'pipeline','results','synthetic','logs'))
root <- file.path(tmp,'pipeline','results','synthetic')
stopifnot(file.exists(file.path(root,'COMPLETE.txt')))
co <- as.data.table(read_parquet(file.path(root,'sensitive','cohort','adult_pneumonia_cohort.parquet')))
qa <- fread(file.path(root,'shareable','ccw_results','covariate_availability_and_processing.csv'))
stopifnot(qa[variable=='glucose',role]=='baseline',qa[variable=='albumin',role]=='baseline',
          !qa[variable=='ph',usable])
long <- as.data.table(read_parquet(file.path(root,'sensitive','analysis','pneumonia_ccw_long.parquet')))
design <- readRDS(file.path(root,'sensitive','analysis','analysis_design.rds'))
stopifnot(!any(grepl('heart_failure|chronic_kidney|lactate|pressor|observed_prior',design$baseline_covariates)),
          !any(grepl('glucose|albumin',design$time_varying_covariates)),
          all(c('glucose_model','albumin_model') %in% design$baseline_covariates),
          all(long[site=='JHH',hospital_academic]==1L),all(long[site=='SMH',hospital_academic]==0L),
          all(long[,.(n=uniqueN(glucose_model)),by=hospital_block_id]$n==1L))
ans <- fread(file.path(root,'shareable','ccw_results','five_method_tidy_estimates.csv'))
stopifnot(all(is.finite(ans[estimand=='RR',estimate])))
cat('PASS: standalone proposal pipeline and baseline lab roles. Synthetic results:',root,'\n')
