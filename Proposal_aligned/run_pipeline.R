#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
folder <- dirname(normalizePath(self))
data_dir <- if(length(args)) normalizePath(args[1],mustWork=TRUE) else normalizePath(file.path(folder,'..','..','rclif_2_1'),mustWork=TRUE)
run_id <- Sys.getenv('PROPOSAL_RUN_ID',paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S')))
if(!grepl('^[A-Za-z0-9_-]+$',run_id)) stop('Invalid PROPOSAL_RUN_ID')
root <- file.path(folder,'results',run_id)
if(dir.exists(root)) stop('Run already exists: ',root)
Sys.setenv(TREATMENT_GRACE_HOURS=Sys.getenv('TREATMENT_GRACE_HOURS','24'),
 PROPOSAL_MAX_MISSING=Sys.getenv('PROPOSAL_MAX_MISSING','0.70'),
 CLIF_DATA_CUTOFF_DATE=Sys.getenv('CLIF_DATA_CUTOFF_DATE','2026-10-31'),
 SUMMARY_SMALL_CELL_THRESHOLD=Sys.getenv('SUMMARY_SMALL_CELL_THRESHOLD','10'))
grace <- as.numeric(Sys.getenv('TREATMENT_GRACE_HOURS')); miss <- as.numeric(Sys.getenv('PROPOSAL_MAX_MISSING'))
if(!is.finite(grace)||grace<6||grace>=672||grace%%6!=0) stop('Grace must be a multiple of 6 hours below 672.')
if(!is.finite(miss)||miss<0||miss>=1) stop('PROPOSAL_MAX_MISSING must be in [0,1).')
Sys.unsetenv('PNEUMONIA_CROSSWALK_PATH')
dir.create(file.path(root,'logs'),recursive=TRUE)
cohort <- file.path(root,'sensitive','cohort');candidate <- file.path(root,'sensitive','steroid_candidates');treatment <- file.path(root,'sensitive','treatment')
run <- function(script,arguments) {
 cat('Running ',script,' ...\n',sep='')
 log <- file.path(root,'logs',paste0(script,'.log'))
 status <- system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(folder,script)),shQuote(arguments)),stdout=log,stderr=log)
 if(status!=0) stop(script,' failed; inspect ',log)
}
writeLines(c(paste('data_source:',data_dir),paste('grace_hours:',grace),paste('maximum_acute_missing_fraction:',miss),
 'Static: age, recorded sex, five proposal comorbidities, hospital type, time-zero location, year/season; site-specific nuisance models.',
 'Glucose/albumin: nearest valid record at/before time zero within the hospital block; fixed baseline covariates.',
 'Physiology: LOCF from admission forward; no backfill and no expiry after a valid measurement.',
 'Excluded clinical adjustments: heart failure, CKD, lactate, vasopressors and historical NCO proxies.',
 'Remaining missing values in retained acute features: reference-median single imputation plus missing indicator.',
 'Hospital type is constant within hospital and absorbed by hospital-specific models.',
 'Original cohort/pooled systemic steroid definition retained; not the CAPE COD-aligned analysis.'),file.path(root,'run_design.txt'))
files <- list.files(folder,pattern='[.]R$',full.names=TRUE)
write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),file.path(root,'source_versions.csv'),row.names=FALSE)
run('00_check_clif_cohort_inputs.R',c(data_dir,file.path(root,'sensitive','input_checks')))
run('01_build_adult_pneumonia_cohort.R',c(data_dir,cohort,'poa_any','6',Sys.getenv('COHORT_SELECTION_SEED','32284'),'hopkins_unspecified'))
run('03_extract_candidate_steroid_administrations.R',c(data_dir,cohort,candidate))
run('04_define_steroid_treatment.R',c(candidate,cohort,treatment,as.character(grace)))
run('05_cohort_summary_eda.R',c(data_dir,cohort,treatment,file.path(root,'shareable','cohort_eda'),Sys.getenv('SUMMARY_SMALL_CELL_THRESHOLD'),Sys.getenv('CLIF_DATA_CUTOFF_DATE')))
run('07_prepare_and_run_five_methods.R',c(data_dir,cohort,treatment,file.path(root,'sensitive','analysis'),file.path(root,'shareable','ccw_results'),Sys.getenv('CCW_MAX_PATIENTS','0'),Sys.getenv('CLIF_DATA_CUTOFF_DATE')))
writeLines('All stages completed.',file.path(root,'COMPLETE.txt'))
cat('Completed. Results: ',root,'\n',sep='')
