#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow);library(data.table)})
f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
root <- dirname(dirname(normalizePath(gsub("~+~", " ", f, fixed=TRUE))))
work <- tempfile("steroid_actions_");dir.create(work)
for(d in c("data","cohort","profile","candidate","treatment"))dir.create(file.path(work,d))
t0 <- as.POSIXct("2024-01-01",tz="UTC")
c <- data.table(patient_id=paste0("p",1:8),hospital_block_id=paste0("b",1:8),
  admission_dttm=t0-86400,discharge_dttm=t0+10*86400,respiratory_failure_start=t0)
cw <- c[,.(patient_id,hospital_block_id,hospitalization_id=paste0("h",1:8))]
write_parquet(c,file.path(work,"cohort","adult_pneumonia_cohort.parquet"))
write_parquet(cw,file.path(work,"cohort","adult_pneumonia_hospitalization_crosswalk.parquet"))
med <- data.table(hospitalization_id=c("h1","h1","h2","h2","h3","h4","h5","h6","h7","h8","h8"),
  med_group=c(rep("steroid",10),"steroids"),
  admin_dttm=t0+c(1,20,1,30,2,3,4,5,6,6,1)*3600,
  med_category=c(rep("hydrocortisone",6),"methylprednisolone","fludrocortisone","hydrocortisone","prednisone","prednisone"),
  med_route_category=c(rep("iv",7),rep("enteral",4)),
  med_dose=c(50,50,50,50,50,49,32,1,50,40,40),
  med_dose_unit=c(rep("mg",6),"mg/kg/dose",rep("mg",4)),
  mar_action_category=c("not_given","given",NA,"given","not_given",rep("given",6)))
write_parquet(med,file.path(work,"data","clif_medication_admin_intermittent.parquet"))
run <- function(script,args) {
  log <- file.path(work,paste0(script,".log"))
  status <- system2(file.path(R.home("bin"),"Rscript"),
    c(shQuote(file.path(root,script)),vapply(args,shQuote,character(1))),stdout=log,stderr=log)
  if(status!=0){cat(readLines(log),sep="\n");stop(script," failed")}
}
run("02_profile_glucocorticoid_sources.R",file.path(work,c("data","profile")))
run("03_extract_candidate_steroid_administrations.R",file.path(work,c("data","cohort","candidate")))
run("04_define_steroid_treatment.R",c(file.path(work,c("candidate","cohort","treatment")),"24"))
x <- as.data.table(read_parquet(file.path(work,"treatment","steroid_administrations_with_treatment_eligibility.parquet")))
q <- x[qualifies_as_treatment==1L]
stopifnot(nrow(x)==10L,nrow(q)==3L,all(q$mar_action_clean=="given"),
  q[patient_id=="p1",hours_from_time_zero]==20,
  q[patient_id=="p2",hours_from_time_zero]==30,
  q[patient_id=="p8",hours_from_time_zero]==6,
  !any(x[mar_action_clean!="given"|is.na(mar_action_clean),qualifies_as_treatment]==1L),
  x[patient_id=="p5",qualification_reason]=="unsupported_or_missing_dose_unit",
  x[patient_id=="p7",qualification_reason]=="route_not_eligible_for_drug")
audit <- fread(file.path(work,"treatment","steroid_administration_action_audit.csv"))
stopifnot(sum(audit$qualifying_rows)==3L)
cat("PASS: 02--04 source checks, singular label, given-only initiation, late initiation, unknown action, dose/unit/route exclusions.\n")
