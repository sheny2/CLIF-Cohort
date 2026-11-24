#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"hospital_block_helpers.R"))
setDTthreads(1)
t0 <- as.POSIXct("2024-01-01",tz="UTC")
# Nested overlap, exact-six-hour gap, >six-hour gap, interleaved patients,
# singleton patient and shuffled source order.
x <- data.table(patient_id=c(rep("p1",5),"p2","p3","p3"),hospitalization_id=paste0("h",1:8),
 admission_dttm=t0+c(0,1,15,26,40,0,0,20)*3600,
 discharge_dttm=t0+c(20,2,16,30,50,2,10,25)*3600)
y <- assign_patient_blocks_base(x[c(8,2,6,4,1,7,5,3)],6)
stopifnot(identical(y[patient_id=="p1",hospital_block_number],c(1L,1L,1L,1L,2L)),
 identical(y[patient_id=="p1",component_order],c(1L,2L,3L,4L,1L)),
 identical(y[patient_id=="p3",hospital_block_number],c(1L,2L)))
set.seed(15)
x <- rbindlist(lapply(seq_len(200),function(i) {
 start <- cumsum(sample(0:30,30,TRUE))
 data.table(patient_id=paste0("p",i),hospitalization_id=paste0("h",i,"_",seq_along(start)),
 admission_dttm=t0+start*3600,discharge_dttm=t0+(start+sample(1:60,30,TRUE))*3600)
}))
reference <- copy(x)
setorder(reference,patient_id,admission_dttm,discharge_dttm,hospitalization_id)
reference[,c("hospital_block_number","gap_hours_from_running_block_end","running_block_end_dttm"):={
 a<-assign_blocks_running_end(admission_dttm,discharge_dttm,6);list(a$block,a$gap,a$running_end)
},by=patient_id]
reference[,hospital_block_id:=paste0(patient_id,"_block_",hospital_block_number)]
reference[,component_order:=seq_len(.N),by=.(patient_id,hospital_block_id)]
result <- assign_patient_blocks_base(x[sample(.N)],6)
stopifnot(isTRUE(all.equal(as.data.frame(result),as.data.frame(reference),check.attributes=FALSE)))
cat("PASS: base-R block assignment matches prior grouped implementation on 6,000 encounters and boundary cases.\n")
