# Nonclinical stress test of matching with repeated timestamps and many queries.
suppressPackageStartupMessages(library(data.table))
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
folder <- dirname(dirname(normalizePath(self)))
source(file.path(folder,'proposal_covariate_helpers.R'))
set.seed(7);blocks<-1000L;per_block<-100L
sid <- rep(as.character(seq_len(blocks)),each=per_block)
st <- as.POSIXct('2025-01-01',tz='UTC')+sample(0:2000,length(sid),TRUE)*3600
qid <- rep(as.character(seq_len(blocks)),each=112)
qt <- as.POSIXct('2025-01-01',tz='UTC')+sample(-10:2100,length(qid),TRUE)*3600
at <- proposal_last_indices(sid,st,qid,qt)
for(i in sample(seq_along(qid),40)) {
 rows <- which(sid==qid[i]&st<=qt[i]);expected<-NA_integer_
 if(length(rows)) expected<-tail(rows[st[rows]==max(st[rows])],1)
 stopifnot(identical(at[i],expected))
}
stopifnot(length(at)==112000L,all(st[at[!is.na(at)]]<=qt[!is.na(at)]))
cat('PASS: 100,000 source records / 112,000 queries, duplicate timestamps, no future matching.\n')
