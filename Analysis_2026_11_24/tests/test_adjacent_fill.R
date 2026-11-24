#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow)})
f<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(gsub('~+~',' ',f,fixed=TRUE))))
source(file.path(root,'variable_rules.R'))
stopifnot(identical(screen_measurements(c(29,30,45,46),'temperature'),c(NA_real_,30,45,NA_real_)),
 identical(screen_measurements(c(0,500,501),'wbc'),c(0,500,NA_real_)),
 identical(screen_measurements(c(0,2000,2007),'glucose'),c(0,2000,NA_real_)))
t0<-as.POSIXct('2024-01-01',tz='UTC')
grid<-data.table(hospital_block_id=c('A','B','C','D'),m=1L,site='s',decision_dttm=t0,
 discharge_dttm=t0+c(30,30,30,-1)*3600)
a<-grid[,.(hospital_block_id,m,decision_dttm,measurement_seconds=NA_real_,value=NA_real_,age_hours=NA_real_)]
x<-data.table(hospital_block_id=c('A','A','B','B','D'),
 measurement_seconds=as.numeric(t0)+c(-20,2,-20,20,-2)*3600,value=c(36,38,35,39,37))
x[,recorded_dttm:=as.POSIXct(measurement_seconds,origin='1970-01-01',tz='UTC')]
work<-tempfile();dir.create(work)
y<-fill_adjacent_measurements(a,x,grid,'temperature',4,work,work)
stopifnot(y[hospital_block_id=='A',value]==38,y[hospital_block_id=='B',value]==35,
 is.na(y[hospital_block_id=='C',value]),is.na(y[hospital_block_id=='D',value]))
audit<-fread(file.path(work,'adjacent_fill_audit_temperature.csv'))
stopifnot(audit$filled_adjacent==2L,audit$filled_future==1L)
cat('PASS: Chad bounds, future nearest fills, earlier tie-break, no cross-block borrowing, no new post-discharge fill, and audit counts.\n')
