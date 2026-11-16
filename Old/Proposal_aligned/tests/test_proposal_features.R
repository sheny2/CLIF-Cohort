suppressPackageStartupMessages(library(data.table))
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(self))),'proposal_covariate_helpers.R'))
t0 <- as.POSIXct('2025-01-01',tz='UTC')
g <- data.table(hospital_block_id='h',m=1:4,admission_dttm=t0-86400,index_dttm=t0,decision_dttm=t0+c(0,6,12,24)*3600)
a <- data.table(hospital_block_id='h',t=t0+c(-2,1,8)*3600,v=c(80,200,500))
b <- proposal_latest(g,a,'t','v','glucose',TRUE)
stopifnot(all(b$glucose==80))
l <- proposal_latest(g,a,'t','v','glucose',FALSE)
stopifnot(identical(l$glucose,c(80,200,500,500)))
a$v <- c(NA,37,39)
t <- proposal_latest(g,a,'t','v','temperature')
stopifnot(is.na(t$temperature[1]),identical(t$temperature[2:4],c(37,39,39)))
# A bad recent observation must not hide an earlier valid reading.
a$v <- c(36,100,39)
t <- proposal_latest(g,a,'t','v','temperature')
stopifnot(t$temperature[2]==36)
# No cross-block contamination; no hospital data from before this admission.
a <- data.table(hospital_block_id=c('h','other'),t=c(t0-2*86400,t0),v=c(36,39))
t <- proposal_latest(g,a,'t','v','temperature')
stopifnot(all(is.na(t$temperature)))
stopifnot(identical(proposal_location(c('ED','ward','IMC','ICU','medical unknown','Missing')),c('ED','Ward','IMC','ICU','Unknown','Unknown')))
# Exclusion threshold is strictly >70%, so 70% is retained, 80% is dropped.
x <- CJ(hospital_block_id=as.character(1:10),m=1:2)
x[,S:=10L]
for(v in c('glucose','albumin','temperature','wbc','ph','paco2','sf_ratio','support_score','nonresp_sofa_proxy')) x[,(v):=NA_real_]
x[hospital_block_id %in% c('1','2','3'),glucose:=c(100,100,110,110,120,120)]
x[hospital_block_id %in% c('1','2'),albumin:=c(2,2,3,3)]
x[,temperature:=36+m/10+as.numeric(hospital_block_id)/10]
p <- proposal_process(x,2,.70)
stopifnot(p$qa[variable=='glucose',usable],!p$qa[variable=='albumin',usable],
          !p$qa[variable=='ph',usable],all(is.finite(p$data$glucose_model)))
# All acute features absent: valid empty term vectors, never a bogus '_model'.
for(v in c('glucose','albumin','temperature')) x[,(v):=NA_real_]
p <- proposal_process(x,2,.70)
stopifnot(length(p$baseline)==0,length(p$tv)==0)
cat('PASS: pre-time-zero baseline labs, indefinite LOCF, no backfill, plausibility filtering, block boundaries, location categories and 70% cutoff.\n')
