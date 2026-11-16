#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow);library(dplyr)})
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"variable_rules.R"));source(file.path(root,"baseline_sf_helpers.R"))
work<-file.path("/tmp",paste0("sf_selection_",Sys.getpid()));dir.create(work)
for(d in c("data","aggregate","sensitive"))dir.create(file.path(work,d))
t0<-as.POSIXct("2024-01-01",tz="UTC")
x<-data.table(patient_id=LETTERS[1:8],hospital_block_id=LETTERS[1:8],source_site_id=rep(c("s1","s2"),4),
 admission_dttm=t0-24*3600,discharge_dttm=t0+5*86400,respiratory_failure_start=t0)
cw<-data.table(patient_id=LETTERS[1:8],hospital_block_id=LETTERS[1:8],hospitalization_id=LETTERS[1:8])
sp<-data.table(hospitalization_id=c("A","A","B","B","C","D","E","F","F","G","H"),
 recorded_dttm=t0+c(0,1,-1,0,1,-13,-1,-2,0,0,-12)*3600,
 vital_category="SpO2",vital_value=c(95,100,95,99,95,95,95,95,0,97,95))
# A blank category must never be interpreted as a saturation measurement.
sp<-rbind(sp,data.table(hospitalization_id="H",recorded_dttm=t0,vital_category=NA_character_,vital_value=99))
fi<-data.table(hospitalization_id=LETTERS[1:8],recorded_dttm=t0+c(0,0,0,0,1,0,0,-12)*3600,
 fio2_set=c(.5,.5,.5,.5,.5,.5,50,.5))
write_parquet(sp,file.path(work,"data","clif_vitals.parquet"))
write_parquet(fi,file.path(work,"data","clif_respiratory_support.parquet"))
y<-select_baseline_sf(file.path(work,"data"),x,cw,file.path(work,"aggregate"),file.path(work,"sensitive"))
stopifnot(setequal(y$cohort$patient_id,c("A","F","G","H")),
 y$cohort[patient_id=="G",baseline_sf]==194,
 y$audit[patient_id=="B",baseline_sf_exclusion]=="Latest valid baseline SpO2 >97%",
 y$audit[patient_id=="C",baseline_sf_exclusion]=="No usable prior SpO2 within 12h",
 y$audit[patient_id=="D",baseline_sf_exclusion]=="No usable prior SpO2 within 12h",
 y$audit[patient_id=="E",baseline_sf_exclusion]=="SpO2 eligible but no usable prior FiO2 within 12h")
encoded<-encode_baseline_sf(y$cohort)
stopifnot(all(is.finite(encoded$cohort$baseline_sf_model)),encoded$processing$missing==0L,
 nrow(x)==8,!"baseline_sf"%in%names(x))
cat("PASS: baseline selection, newer high saturation, no future eligibility, stale measurements, screening, percent FiO2 and inclusive 12h/97% boundaries.\n")
