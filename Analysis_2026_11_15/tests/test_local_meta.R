#!/usr/bin/env Rscript
# Optional argument: a completed synthetic 07 run root for integration testing.
suppressPackageStartupMessages({library(data.table);library(metafor)})
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE))))
source(file.path(root,"local_meta_helpers.R"))
x<-data.table(site=c("A","B","C"),estimand="risk_difference",estimate=c(.1,.2,.3),se=c(.05,.1,.2),site_valid=TRUE,site_status="ok")
p<-pool_local_contrasts(x)
w<-1/x$se^2
stopifnot(abs(p$summary[meta_method=="fixed_effect",estimate]-sum(w*x$estimate)/sum(w))<1e-12,
 all(abs(p$weights[included==TRUE,sum(weight_percent),by=meta_method]$V1-100)<1e-10))
y<-copy(x);y[,`:=`(estimand="log_risk_ratio",estimate=log(c(.8,1.2,1.5)))]
p<-pool_local_contrasts(y)
stopifnot(abs(p$summary[meta_method=="fixed_effect",estimate]-exp(sum(w*y$estimate)/sum(w)))<1e-12)
x[site=="B",se:=0];x[site=="C",site_valid:=FALSE]
p<-pool_local_contrasts(x)
stopifnot(all(p$summary$k==1),all(is.na(p$summary$estimate)),all(p$weights[site!="A",included]==FALSE))
cat("PASS: fixed-effect arithmetic, log-ratio scale, weights, exclusions and insufficient-site handling.\n")
a<-commandArgs(TRUE)
if(length(a)) {
 run_root<-a[1];original<-file.path(run_root,"aggregate","analysis","five_method_results.csv")
 before<-tools::md5sum(original)
 log<-tempfile(fileext=".log")
 status<-system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(root,"08_local_ccw_meta_analysis.R")),shQuote(run_root)),stdout=log,stderr=log)
 if(status!=0){cat(readLines(log),sep="\n");stop("08 integration failed")}
 stopifnot(identical(before,tools::md5sum(original)))
 out<-file.path(run_root,"aggregate","local_ccw_meta")
 sites<-fread(file.path(out,"local_ccw_site_results.csv"))
 six<-fread(file.path(out,"six_method_results.csv"));meta<-fread(file.path(out,"meta_analysis_results.csv"))
 stopifnot(nrow(sites)==5L,nrow(six)==6L,nrow(meta)==8L,
  isTRUE(all.equal(fread(original),six[method_key!="local_ccw_meta"],check.attributes=FALSE)))
 saved<-readRDS(file.path(run_root,"sensitive","local_ccw_meta","local_ccw_meta_fits.rds"))
 for(s in names(saved$sites)) {
  c<-saved$sites[[s]]$component
  # Independently reconstruct the weighted KM curves from local clone arrays.
  curve<-function(W,R,D) {
   risk<-colSums(W*R);event<-colSums(W*D)
   c(1,cumprod(1-ifelse(risk>0,event/risk,0)))
  }
  s1<-curve(c$W1,c$R1,c$D1);s0<-curve(c$W0,c$R0,c$D0)
  rd<-tail(s0,1)-tail(s1,1)
  rmst<-sum(head(s1,-1)-head(s0,-1))*saved$design$step_hours/24
  stopifnot(abs(sites[site==s,risk_difference]-rd)<1e-12,abs(sites[site==s,rmst_difference_days]-rmst)<1e-12)
 }
 weights<-fread(file.path(out,"meta_analysis_site_weights_and_exclusions.csv"))
 sums<-weights[included==TRUE & is.finite(weight_percent),sum(weight_percent),by=.(estimand,meta_method)]
 stopifnot(all(abs(sums$V1-100)<1e-8),file.exists(file.path(out,"six_method_effect_comparison.png")))
 site_forest<-fread(file.path(out,"local_ccw_site_forest_data.csv"))
 stopifnot(all(grepl("^Overall CCW MA",site_forest[,tail(method,1),by=estimand]$V1)),
  nrow(site_forest)==12L,file.exists(file.path(out,"local_ccw_site_risk_ratio_forest.png")))
 cat("PASS: local weighted curves, six-method output, site forests with overall CCW MA last, and unchanged 07 results.\n")
}
