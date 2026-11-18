#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow);library(data.table)})
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root <- dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"hospital_block_helpers.R"));source(file.path(root,"history_eda_helpers.R"))
t0 <- as.POSIXct("2024-06-01",tz="UTC")
c <- data.table(patient_id=LETTERS[1:4],hospital_block_id=paste0("i",1:4),
  site=c("JHH","BMC","JHH","BMC"),S=c(1L,5L,113L,1L),admission_dttm=t0)
cw <- c[,.(patient_id,hospital_block_id,hospitalization_id=hospital_block_id)]
h <- data.table(patient_id=c$patient_id,hospitalization_id=c$hospital_block_id,
  admission_dttm=t0,discharge_dttm=t0+5*86400,admission_type_category="ed")
prior <- data.table(patient_id=c("A","A","A","A","A","B","B","D"),
  hospitalization_id=paste0("h",1:8),
  admission_dttm=t0+c(-10,-9+4/24,-400,2,-2,-50,-30,-365)*86400,
  discharge_dttm=t0+c(-9,-7,-360,3,-1/24,-40,NA,-364)*86400,
  admission_type_category=c("ed","osh","direct","ed","ed",NA,"ed","direct"))
h <- rbind(h,prior)
x <- build_history_eda(c,cw,h)
z <- x$patient[window=="Previous 365 days"&encounter_type=="All hospitalizations"]
stopifnot(z[patient_id=="A",observed_episode_count]==1L,z[patient_id=="A",total_los_days]==3,
  z[patient_id=="C",observed_episode_count]==0L,is.na(z[patient_id=="C",mean_los_days]),
  z[patient_id=="B",invalid_possible_prior_rows]==1L,z[patient_id=="D",observed_episode_count]==1L)
stopifnot(x$patient[patient_id=="A"&window=="All available prior history"&encounter_type=="All hospitalizations",observed_episode_count]==2L,
  x$patient[patient_id=="B"&window=="Previous 365 days"&encounter_type=="ED-origin hospitalizations",unknown_ed_episodes]==1L,
  x$overlap[patient_id=="A"&encounter_type=="All hospitalizations",observed_days_in_previous_year]==8)
args <- commandArgs(TRUE)
work <- if(length(args))args[1L] else tempfile("history_eda_")
dir.create(work,recursive=TRUE)
inp <- file.path(work,"sensitive","analysis");dir.create(inp,recursive=TRUE)
agg <- file.path(work,"aggregate","analysis");dir.create(agg,recursive=TRUE)
writeLines("Synthetic",file.path(agg,"COMPLETE.txt"))
write_parquet(c,file.path(inp,"pneumonia_ccw_cohort.parquet"))
write_parquet(cw,file.path(inp,"sampled_hospital_block_crosswalk.parquet"))
saveRDS(list(tau_intervals=4L),file.path(inp,"analysis_design.rds"))
hist <- file.path(work,"expanded.parquet");write_parquet(h,hist)
before <- tools::md5sum(file.path(agg,"COMPLETE.txt"));log <- file.path(work,"11.log")
status <- system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(root,"11_prior_utilization_eda.R")),shQuote(work),shQuote(hist)),stdout=log,stderr=log)
if(status!=0){cat(readLines(log),sep="\n");stop("11 integration failed")}
out <- file.path(work,"aggregate","prior_utilization_eda")
stopifnot(file.exists(file.path(out,"COMPLETE.txt")),length(list.files(out,pattern="[.]png$"))==5L,
  identical(before,tools::md5sum(file.path(agg,"COMPLETE.txt"))))
r <- fread(file.path(out,"history_availability_and_prevalence.csv"))
stopifnot(all(r[scope=="Overall",patients]==4L),!any(c("patient_id","hospitalization_id")%in%names(r)))
cat("PASS: merged encounters, index/future exclusion, annual boundary, overlapping-year duration, unknown ED, invalid dates, zeros, and 11 integration.\nOutputs:",out,"\n")
