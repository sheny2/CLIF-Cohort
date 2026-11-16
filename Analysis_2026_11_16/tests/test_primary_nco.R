#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow)})
f<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);root<-dirname(dirname(normalizePath(f)))
source(file.path(root,'hospital_block_helpers.R'));source(file.path(root,'nco_history_helpers.R'))
t0<-as.POSIXct('2024-06-01',tz='UTC')
c<-data.table(patient_id=LETTERS[1:5],hospital_block_id=paste0('index',1:5),site='s',admission_dttm=t0)
cw<-data.table(patient_id=LETTERS[1:5],hospital_block_id=c$hospital_block_id,hospitalization_id=paste0('i',1:5))
h<-data.table(patient_id=LETTERS[1:5],hospitalization_id=cw$hospitalization_id,admission_dttm=t0,discharge_dttm=t0+86400)
h<-rbind(h,data.table(patient_id=c('A','B','C','D','E'),hospitalization_id=paste0('h',1:5),
 admission_dttm=t0+c(-365,-366,5,-3,-100)*86400,
 discharge_dttm=as.POSIXct(as.numeric(t0)+c(-360,-360,6,-1/24,NA)*86400,origin='1970-01-01',tz='UTC')))
x<-build_prior_admission_nco(c,cw,h)
stopifnot(identical(x$cohort$prior_year_admission,c(1L,0L,0L,0L,NA_integer_)))
y<-build_prior_admission_nco(c,cw,h,as.Date('2024-01-01'))
stopifnot(!any(y$cohort$nco_eligible))
# Identical clones cancel at the patient influence level, rather than falsely
# treating the two copies as independent samples.
z<-nco_weighted_contrast(c(0,1,0,1),rep(1,4),rep(1,4))
stopifnot(z$difference==0,z$se==0,z$ess_early==4)
cat('PASS: prior-year boundary, future exclusion, index-connected episode exclusion, missing-date handling, coverage and clone covariance.\n')
a<-commandArgs(TRUE)
if(length(a)) {
 run<-a[1];data<-file.path(dirname(run),'data');work<-tempfile('nco_expanded_');dir.create(work)
 h<-as.data.table(read_parquet(file.path(data,'clif_hospitalization.parquet')))[,.(patient_id,hospitalization_id,admission_dttm,discharge_dttm)]
 prior<-copy(h)[seq_len(.N)%%2==0]
 prior[,`:=`(hospitalization_id=paste0('prior_',hospitalization_id),admission_dttm=admission_dttm-100*86400,discharge_dttm=admission_dttm-95*86400)]
 path<-file.path(work,'expanded.parquet');write_parquet(rbind(h,prior),path)
 original<-file.path(run,'aggregate','analysis','five_method_results.csv');before<-tools::md5sum(original)
 log<-file.path(work,'test.log')
 status<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(root,'09_prior_hospitalization_nco_diagnostic.R')),shQuote(run),shQuote(path),'2020-01-01'),stdout=log,stderr=log)
 if(status!=0){cat(readLines(log),sep='\n');stop('09 integration failed')}
 stopifnot(identical(before,tools::md5sum(original)))
 out<-fread(file.path(run,'aggregate','nco_prior_hospitalization','primary_nco_weighted_contrasts.csv'))
 stopifnot(any(is.finite(out$difference)),any(out$assessment=='Landmark survivor comparison'))
 cat('PASS: 09 reads expanded history, reuses saved CCW/landmark weights, and leaves 07 estimates unchanged.\n')
}
