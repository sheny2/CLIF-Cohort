suppressPackageStartupMessages({library(arrow);library(data.table)})
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1])
source_folder <- dirname(normalizePath(self))
t <- tempfile();dir.create(t);dir.create(file.path(t,'historical_proxy_sensitivity'))
for(f in c('run_analysis.R','history_helpers.R')) file.copy(file.path(source_folder,f),file.path(t,'historical_proxy_sensitivity',f))
dir.create(file.path(t,'data'));dir.create(file.path(t,'cohort'))
tm <- function(x) as.POSIXct(x,tz='UTC')
put <- function(x,n) write_parquet(x,file.path(t,'data',paste0('clif_',n,'.parquet')))
write_parquet(data.table(hospital_block_id=c('b1','b2'),patient_id=c('p1','p2'),admission_dttm=tm(c('2021-06-01','2022-06-01'))),file.path(t,'cohort','adult_pneumonia_cohort.parquet'))
put(data.table(hospitalization_id='h1',patient_id='p1',admission_dttm=tm('2021-02-01'),discharge_dttm=tm('2021-02-03'),admission_type_category='ed'),'hospitalization')
put(data.table(hospitalization_id='h1',in_dttm=tm('2021-02-01'),location_category='icu'),'adt')
put(data.table(hospitalization_id='h1',recorded_dttm=tm('2021-02-01'),device_category='imv'),'respiratory_support')
put(data.table(hospitalization_id='h1',diagnosis_code='J96.01',diagnosis_code_format='ICD10CM'),'hospital_diagnosis')
writeLines(c('rm(list=ls())','args <- commandArgs(trailingOnly=TRUE)','grace_hours <- 24L','stopifnot(grace_hours == 96L)',
 'long <- data.table(hospital_block_id=c("b1","b2"),site="s1",m=1L,age_z=c(0,1))',
 'baseline_covariates <- c("age_z")',
 'baseline_covariates <- baseline_covariates[vapply(long[, ..baseline_covariates], function(x) uniqueN(x)>1, logical(1))]',
 'dir.create(args[4],recursive=TRUE);dir.create(args[5],recursive=TRUE)',
 'fwrite(data.table(method="fed_ccw",estimand="RR",estimate=length(baseline_covariates)),file.path(args[5],"five_method_tidy_estimates.csv"))'),file.path(t,'07_prepare_and_run_five_methods.R'))
writeLines('# fixture',file.path(t,'06_generalized_five_methods.R'))
Sys.setenv(HISTORY_MODE='observed_extract',HISTORY_COVERAGE_START='2021-01-01',TREATMENT_GRACE_HOURS='96')
result <- system2(file.path(R.home('bin'),'Rscript'),c(file.path(t,'historical_proxy_sensitivity','run_analysis.R'),file.path(t,'data'),file.path(t,'cohort'),file.path(t,'treatment'),file.path(t,'data'),file.path(t,'out')),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(result,'status'))) {cat(paste(result,collapse='\n'));stop('integration failed')}
z <- fread(file.path(t,'out','shareable','reference_vs_augmented.csv'))
stopifnot(z$estimate_augmented > z$estimate_reference)
p <- as.data.table(read_parquet(file.path(t,'out','sensitive','historical_proxies.parquet')))
stopifnot(p[hospital_block_id=='b1',observed_prior_imv]==1L,p[hospital_block_id=='b1',prior_history_short_lookback]==1L)
unlink(t,recursive=TRUE)
cat('PASS: synthetic Parquet wrapper integration (lightweight parent fixture).\n')
