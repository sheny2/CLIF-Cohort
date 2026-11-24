#!/usr/bin/env Rscript
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root <- dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"06_generalized_five_methods.R"))
M <- 8L;tau <- 4L
ids <- c("discharge12","discharge30","treated_discharge12","death20","hospice12","late_start","no_start","discharge24","discharge25","admin_censor")
S <- c(Inf,Inf,1,Inf,Inf,6,Inf,Inf,Inf,Inf)
end <- c(2,5,2,8,8,8,8,4,5,8)
exempt <- c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE)
terminal <- c(8,8,8,4,2,8,8,8,8,3)
d <- do.call(rbind,lapply(seq_along(ids),function(i)data.frame(site="A",id=ids[i],m=seq_len(terminal[i]),
  event=as.integer(seq_len(terminal[i])==terminal[i]&i%in%c(4,5)),
  censor=as.integer(i==10L&seq_len(terminal[i])==terminal[i]),
  init=as.integer(seq_len(terminal[i])==S[i]),adherence_end_m=end[i],early_alive_discharge_untreated=exempt[i])))
p <- g6_prepare_data(d,"site","id","m","event","censor","init",NULL,character(),character(),M,0)
learner <- list(fit=function(formula,data)list(p=if(all.vars(formula)[1L]==".init") .3 else
  if(".trt"%in%all.vars(formula)) .2 else .1),predict=function(model,newdata)rep(model$p,nrow(newdata)))
cfg <- list(tau=tau,horizon=M,terms=list(initiation_denominator=".time",initiation_numerator=".time",
  censoring_denominator=c(".time",".trt"),censoring_numerator=".time"),learners=list(initiation=learner,censoring=learner),trunc=c(0,1))
x <- g6_local_components(p$data,cfg,rep(.1,M))
at <- function(id)match(as.character(unique(p$data$.uid[p$data$.id==id])),x$uid)
i <- at("discharge12")
stopifnot(all(x$R1[i,]),all(x$R0[i,]),all(is.na(x$Hden[i,3:M])),
  all(x$raw_W1[i,2:M]==x$raw_W1[i,2]),all(x$raw_W0[i,2:M]==x$raw_W0[i,2]))
stopifnot(all(x$R1[at("discharge24"),]),all(x$R0[at("discharge24"),]),
  !x$R1[at("discharge25"),5],all(x$R0[at("discharge25"),]),
  !x$R1[at("discharge30"),5],x$R1[at("discharge30"),4],
  all(x$R1[at("treated_discharge12"),]),!any(x$R0[at("treated_discharge12"),]),
  x$D1[at("death20"),4],x$D0[at("death20"),4],
  x$D1[at("hospice12"),2],x$D0[at("hospice12"),2],
  all(x$R0[at("late_start"),]),x$R1[at("no_start"),4],!x$R1[at("no_start"),5])
stopifnot(is.finite(g6_central_ccw(list(g6_local_stats(x)))$RR))
cat("PASS: both clones retained for untreated early discharge; weights frozen; initiation history, deadline, deaths, hospice, and late initiation handled.\n")
