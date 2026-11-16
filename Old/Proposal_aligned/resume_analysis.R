#!/usr/bin/env Rscript
# Rerun only 07 with the saved run's settings and unchanged cohort/treatment.
a <- commandArgs(trailingOnly=TRUE)
if(!length(a)||length(a)>2) stop('Usage: Rscript Proposal_aligned/resume_analysis.R EXISTING_RUN_DIR [DATA_DIR]')
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
folder <- dirname(normalizePath(self));root <- normalizePath(a[1],mustWork=TRUE)
lines <- readLines(file.path(root,'run_design.txt'))
setting <- function(key) {
 prefix <- paste0(key,':')
 x <- lines[startsWith(lines,prefix)]
 if(length(x)!=1L) stop('Missing/ambiguous saved setting: ',key)
 trimws(substring(x,nchar(prefix)+1L))
}
data_dir <- if(length(a)==2) a[2] else setting('data_source')
Sys.setenv(TREATMENT_GRACE_HOURS=setting('grace_hours'),
           PROPOSAL_MAX_MISSING=setting('maximum_acute_missing_fraction'))
Sys.unsetenv('PNEUMONIA_CROSSWALK_PATH')
# The original runner used this documented default unless explicitly overridden.
cutoff_file <- file.path(root,'shareable','cohort_eda','cohort_overall_summary.csv')
cutoff <- if(file.exists(cutoff_file)) {
 x <- read.csv(cutoff_file,stringsAsFactors=FALSE)
 if(all(c('measure','value') %in% names(x))) x$value[x$measure=='data_cutoff'] else character()
} else character()
# Most original EDA files name the first column "metric"; use either label.
if(length(cutoff)!=1L && file.exists(cutoff_file)) {
 x <- read.csv(cutoff_file,stringsAsFactors=FALSE)
 cutoff <- as.character(x[[2]][as.character(x[[1]])=='data_cutoff'])
}
if(length(cutoff)!=1L || is.na(as.Date(cutoff))) stop('Cannot recover original cutoff from cohort_overall_summary.csv; inspect original run.')
threshold <- as.character(x[[2]][as.character(x[[1]])=='minimum_reportable_cell'])
if(length(threshold)==1L && is.finite(as.numeric(threshold))) Sys.setenv(SUMMARY_SMALL_CELL_THRESHOLD=threshold)
log <- file.path(root,'logs',paste0('07_resumed_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.log'))
files <- list.files(folder,pattern='[.]R$',full.names=TRUE)
write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),file.path(root,'analysis_resume_source_versions.csv'),row.names=FALSE)
if(file.exists(file.path(root,'COMPLETE.txt'))) unlink(file.path(root,'COMPLETE.txt'))
cat('Resuming 07: grace=',setting('grace_hours'),' hours; missingness cutoff=',setting('maximum_acute_missing_fraction'),'\nLog: ',log,'\n',sep='')
status <- system2(file.path(R.home('bin'),'Rscript'),shQuote(c(file.path(folder,'07_prepare_and_run_five_methods.R'),
 data_dir,file.path(root,'sensitive','cohort'),file.path(root,'sensitive','treatment'),
 file.path(root,'sensitive','analysis'),file.path(root,'shareable','ccw_results'),'0',cutoff)),stdout=log,stderr=log)
if(status!=0L) {
 cat(paste(tail(readLines(log,warn=FALSE),20),collapse='\n'),'\n')
 stop('Analysis failed; inspect ',log)
}
writeLines('All stages completed; 07 resumed using saved cohort/treatment and saved design.',file.path(root,'COMPLETE.txt'))
cat('Analysis complete. Results: ',file.path(root,'shareable','ccw_results'),'\n',sep='')
