#!/usr/bin/env Rscript
f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
root <- dirname(dirname(normalizePath(gsub("~+~", " ", f, fixed=TRUE))))
source(file.path(root,"analysis_preparation.R"))
for(step in c(1,2,3,4,6,8,12,24)) {
 x <- analysis_timing(24,step)
 stopifnot(x$step_hours==step,x$tau*step==24,x$M*step==672)
}
for(step in c("", "bad", NA, Inf, 0, -3, 1.5, 5, 7, 48, 672))
 stopifnot(inherits(try(analysis_timing(24,step),silent=TRUE),"try-error"))
for(grace in c(NA,Inf,0,5,25,672))
 stopifnot(inherits(try(analysis_timing(grace,6),silent=TRUE),"try-error"))
stopifnot(analysis_timing(96,12)$tau==8L)
cat("PASS: 1/2/3/4/6/8/12/24-hour intervals, 96-hour grace, invalid timing rejected.\n")
