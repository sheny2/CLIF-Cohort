#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr)})
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=3L) stop("Usage: 02_define_hydrocortisone.R DATA_DIR COHORT_DIR TREATMENT_DIR")
self <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(self)),"cape_cod_helpers.R"))
grace <- as.numeric(Sys.getenv("TREATMENT_GRACE_HOURS","24"))
co <- as.data.table(read_parquet(file.path(a[2],"adult_pneumonia_cohort.parquet")))
cw <- as.data.table(read_parquet(file.path(a[2],"adult_pneumonia_hospitalization_crosswalk.parquet")))
x <- cape_steroid_records(a[1],unique(as.character(cw$hospitalization_id)))
x <- merge(x,cw[,.(hospitalization_id,hospital_block_id,patient_id)],by="hospitalization_id")
x <- merge(x,co[,.(hospital_block_id,respiratory_failure_start,discharge_dttm)],by="hospital_block_id")
x[,hours_from_time_zero:=as.numeric(difftime(admin_dttm,respiratory_failure_start,units="hours"))]
x <- x[hours_from_time_zero>=0 & hours_from_time_zero<=672 & admin_dttm<=discharge_dttm]
x[,`:=`(med_category_clean=drug,med_route_clean=route)]
q <- x[qualifies_hydrocortisone==1L]
first <- q[order(admin_dttm),.SD[1L],by=hospital_block_id]
assign <- copy(co[,.(patient_id,hospital_block_id,respiratory_failure_start)])
assign[,`:=`(treatment_grace_hours=grace,ever_qualifying_steroid=as.integer(hospital_block_id %in% q$hospital_block_id))]
early <- first[hours_from_time_zero<=grace,.(hospital_block_id,treatment_admin_dttm=admin_dttm,
              treatment_hours_from_time_zero=hours_from_time_zero,treatment_med_category_clean=drug,
              treatment_dose_mg=dose_mg,treatment_med_route_clean=route)]
assign <- merge(assign,early,by="hospital_block_id",all.x=TRUE)
assign[,treated_within_grace:=as.integer(!is.na(treatment_admin_dttm))]
# The pragmatic policy permits subsequent usual-care steroids AFTER the grace
# window. Other systemic steroid regimens within the grace window are deviations.
dev <- x[qualifies_hydrocortisone==0L & hours_from_time_zero<=grace,
         .(protocol_deviation_dttm=min(admin_dttm)),by=hospital_block_id]
dir.create(a[3],recursive=TRUE,showWarnings=FALSE)
write_parquet(x,file.path(a[3],"all_systemic_steroid_records.parquet"))
write_parquet(q,file.path(a[3],"qualifying_steroid_administrations.parquet"))
write_parquet(assign,file.path(a[3],"steroid_treatment_assignment.parquet"))
write_parquet(dev,file.path(a[3],"protocol_deviations.parquet"))
# Patient-level drug/time provenance stays sensitive; only aggregate reports are exported by 05/07.
fwrite(x[,.(administrations=.N,patients=uniqueN(patient_id)),by=.(source_table,unit,qualifies_hydrocortisone)],
       file.path(a[3],"medication_source_audit_sensitive.csv"))
