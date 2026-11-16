#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"outcome_helpers.R"))
t0<-as.POSIXct("2024-01-01",tz="UTC")
x<-data.table(hospital_block_id=letters[1:9],admission_dttm=t0-86400,respiratory_failure_start=t0,
 discharge_dttm=t0+c(5,5,5,5,30,5,5,5,28)*86400,
 discharge_category=c("Home","Home","Expired","Hospice","Hospice","Expired","Home","Home","Hospice"),
 death_dttm=as.POSIXct(as.numeric(t0)+c(15,NA,NA,15,NA,15,4,-2,NA)*86400,origin="1970-01-01",tz="UTC"))
y<-derive_primary_outcome(x)
stopifnot(is.na(y$event_dttm[1]),is.na(y$event_dttm[2]),is.na(y$event_dttm[8]),
 y$event_hours[3]==120,y$event_hours[4]==120,y$event_hours[5]==720,
 y$event_hours[6]==120,y$event_hours[7]==96,y$event_hours[9]==672,
 y$event_type[4]=="hospice_discharge",y$event_type[7]=="in_hospital_death",
 sum(is.finite(y$event_hours)&y$event_hours<=672)==5L)
# The source table is not modified by the outcome helper.
stopifnot(!"event_dttm"%in%names(x))
cat("PASS: in-hospital death, hospice, expired fallback, excluded post-discharge death, and day-28 boundary.\n")
