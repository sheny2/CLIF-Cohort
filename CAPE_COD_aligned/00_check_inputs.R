#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2) stop("Usage: 00_check_inputs.R DATA_DIR OUTPUT_DIR")
self <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(self)),"cape_cod_helpers.R"))
mode <- Sys.getenv("CAPE_SEVERITY_MODE","pf")
hydro <- Sys.getenv("CAPE_HYDRO_MODE","pragmatic_iv")
if(!mode %in% c("pf","sf_proxy")) stop("Invalid CAPE_SEVERITY_MODE")
if(!hydro %in% c("pragmatic_iv","continuous_only")) stop("Invalid CAPE_HYDRO_MODE")
requirements <- list(
 patient=c("patient_id","death_dttm","sex_category"),
 hospitalization=c("hospitalization_id","patient_id","admission_dttm","discharge_dttm","age_at_admission","admission_type_category"),
 hospital_diagnosis=c("hospitalization_id","diagnosis_code","diagnosis_code_format","diagnosis_primary","poa_present"),
 adt=c("hospitalization_id","hospital_id","in_dttm","out_dttm","location_category"),
 respiratory_support=c("hospitalization_id","recorded_dttm","device_category","fio2_set","lpm_set","tracheostomy",Sys.getenv("CAPE_PEEP_COLUMN","peep_set")),
 vitals=c("hospitalization_id","recorded_dttm","vital_category","vital_value"),
 labs=c("hospitalization_id","lab_collect_dttm","lab_category","lab_value_numeric"),
 code_status=c("patient_id","start_dttm","code_status_category"),
 patient_assessments=c("hospitalization_id","recorded_dttm","assessment_category","numerical_value"),
 medication_admin_continuous=c("hospitalization_id","admin_dttm","med_category","med_dose","med_dose_unit","mar_action_category"),
 medication_admin_intermittent=c("hospitalization_id","admin_dttm","med_category","med_dose","med_dose_unit","med_route_category"))
audit <- rbindlist(lapply(names(requirements),function(nm){
 path <- file.path(a[1],paste0("clif_",nm,".parquet"))
 schema <- open_dataset(path)$schema$names
 missing <- setdiff(requirements[[nm]],schema)
 if(length(missing)) stop(nm," missing: ",paste(missing,collapse=", "))
 data.table(table=nm,column=schema)
}))
dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
fwrite(audit,file.path(a[2],"input_schema.csv"))
path <- file.path(a[1],"clif_labs.parquet")
cat <- open_dataset(path) |> distinct(lab_category) |> collect() |> as.data.table()
fwrite(cat,file.path(a[2],"lab_categories.csv"))
po <- strsplit(Sys.getenv("CAPE_PAO2_CATEGORIES","po2_arterial;pao2;pao2_arterial;arterial_po2"),";",fixed=TRUE)[[1]]
if(mode=="pf" && !any(tolower(cat$lab_category) %in% tolower(po))) stop(
 "No configured arterial PO2 category found. Review lab_categories.csv and set CAPE_PAO2_CATEGORIES, or explicitly use CAPE_SEVERITY_MODE=sf_proxy.")
writeLines(c(paste("severity_mode:",mode),paste("hydrocortisone_mode:",hydro),
 "PEEP is in cm H2O; PaO2 is assumed to be standardized to mmHg in CLIF. Verify source units.",
 "No supplemental RCT protocol has been operationalized beyond the documented observable criteria."),file.path(a[2],"input_contract.txt"))
cat("Schema and category preflight passed.\n")
