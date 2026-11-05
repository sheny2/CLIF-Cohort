#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
self <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
folder <- dirname(normalizePath(self))
data_dir <- if(length(args)) normalizePath(args[1],mustWork=TRUE) else
  normalizePath(file.path(folder,"..","..","rclif_2_1"),mustWork=TRUE)
run_id <- Sys.getenv("CAPE_RUN_ID",paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S")))
if(!grepl("^[A-Za-z0-9_-]+$",run_id)) stop("CAPE_RUN_ID may contain letters, numbers, underscores or hyphens only.")
root <- file.path(folder,"results",run_id)
if(dir.exists(root)) stop("Run exists: ",root,". Choose a fresh CAPE_RUN_ID.")
dir.create(root,recursive=TRUE)
Sys.setenv(TREATMENT_GRACE_HOURS=Sys.getenv("TREATMENT_GRACE_HOURS","24"),
           CLIF_DATA_CUTOFF_DATE=Sys.getenv("CLIF_DATA_CUTOFF_DATE","2026-10-31"),
           SUMMARY_SMALL_CELL_THRESHOLD=Sys.getenv("SUMMARY_SMALL_CELL_THRESHOLD","10"),
           CAPE_SEVERITY_MODE=Sys.getenv("CAPE_SEVERITY_MODE","pf"),
           CAPE_HYDRO_MODE=Sys.getenv("CAPE_HYDRO_MODE","pragmatic_iv"))
# Explicitly prevent inherited crosswalk overrides from another analysis.
Sys.unsetenv("PNEUMONIA_CROSSWALK_PATH")
grace <- suppressWarnings(as.numeric(Sys.getenv("TREATMENT_GRACE_HOURS")))
if(!is.finite(grace)||grace<6||grace>=672||grace%%6!=0) stop("Invalid grace period; use a multiple of six hours below 672.")
sensitive <- file.path(root,"sensitive");shareable <- file.path(root,"shareable")
cohort <- file.path(sensitive,"cohort"); treatment <- file.path(sensitive,"treatment")
dir.create(file.path(root,"logs"))
run <- function(script,arguments) {
 cat("Running ",script," ...\n",sep="")
 log <- file.path(root,"logs",paste0(script,".log"))
 status <- system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(folder,script)),shQuote(arguments)),stdout=log,stderr=log)
 if(status!=0) stop(script," failed. Inspect ",log,". Original analysis is unchanged.")
}
writeLines(c(paste("source:",data_dir),paste("grace_hours:",grace),paste("severity_mode:",Sys.getenv("CAPE_SEVERITY_MODE")),
 paste("hydro_mode:",Sys.getenv("CAPE_HYDRO_MODE")),paste("cutoff:",Sys.getenv("CLIF_DATA_CUTOFF_DATE")),
 "Target: early IV hydrocortisone initiation versus no systemic steroids within grace, subsequent usual care allowed.",
 "Time zero: first simultaneous ICU and observed respiratory eligibility, within 48h of admission.",
 "Pragmatic initiation estimand; not exact trial assignment or sustained dose/taper replication."),file.path(root,"run_design.txt"))
files <- list.files(folder,pattern="[.]R$",full.names=TRUE)
write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),file.path(root,"source_versions.csv"),row.names=FALSE)
run("00_check_inputs.R",c(data_dir,file.path(sensitive,"input_checks")))
run("01_build_adult_pneumonia_cohort.R",c(data_dir,cohort,"poa_any","6",Sys.getenv("COHORT_SELECTION_SEED","32284"),"hopkins_unspecified"))
run("02_define_hydrocortisone.R",c(data_dir,cohort,treatment))
run("05_cohort_summary_eda.R",c(data_dir,cohort,treatment,file.path(shareable,"cohort_eda"),Sys.getenv("SUMMARY_SMALL_CELL_THRESHOLD"),Sys.getenv("CLIF_DATA_CUTOFF_DATE")))
run("07_prepare_and_run_five_methods.R",c(data_dir,cohort,treatment,file.path(sensitive,"analysis"),file.path(shareable,"ccw_results"),Sys.getenv("CCW_MAX_PATIENTS","0"),Sys.getenv("CLIF_DATA_CUTOFF_DATE")))
writeLines("All pipeline stages completed.",file.path(root,"COMPLETE.txt"))
cat("Completed. Results: ",root,"\n",sep="")
