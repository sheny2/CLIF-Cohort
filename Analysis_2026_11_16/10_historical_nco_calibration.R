#!/usr/bin/env Rscript
# Run after 09: Rscript 10_historical_nco_calibration.R RUN_ROOT EXPANDED_HOSPITALIZATION_PARQUET [SYSTEM_COVERAGE_START_DATE]
# Exploratory transport of historical log-prevalence-ratio systematic error
# to the primary CCW log RR. Does not identify a corrected causal effect.
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2)})
setDTthreads(1L)
a<-commandArgs(TRUE)
if(length(a)<2L)stop("Usage: 10 RUN_ROOT EXPANDED_HOSPITALIZATION_PARQUET [SYSTEM_COVERAGE_START_DATE]")
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
root<-dirname(normalizePath(gsub("~+~"," ",f,fixed=TRUE)))
for(h in c("06_generalized_five_methods.R","hospital_block_helpers.R","site_helpers.R",
           "nco_history_helpers.R","nco_panel_helpers.R","nco_calibration_helpers.R"))source(file.path(root,h))
run<-normalizePath(a[1]);expanded<-normalizePath(a[2])
coverage<-if(length(a)>=3L)as.Date(a[3]) else as.Date(NA)
if(length(a)>=3L&&is.na(coverage))stop("Coverage start must be YYYY-MM-DD")
inp<-file.path(run,"sensitive","analysis")
for(p in c(file.path(run,"aggregate","analysis","COMPLETE.txt"),
           file.path(run,"aggregate","nco_prior_hospitalization","COMPLETE.txt")))
  if(!file.exists(p))stop("Run updated 07 and 09 first; missing: ",p)
out<-file.path(run,"aggregate","nco_calibration")
sensitive<-file.path(run,"sensitive","nco_calibration")
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
unlink(file.path(out,"COMPLETE.txt"))
# Prevent obsolete successful-fit artifacts surviving a later unsuccessful run.
unlink(file.path(out,c("nn_mnn_bias_distributions.png","fitted_bias_density_plot_data.csv",
  "leave_one_control_out.csv","nco_mnn_membership.csv","exploratory_ccw_rr_forest.png")))
write<-function(x,f) {
  if(any(c("patient_id","hospitalization_id","hospital_block_id","uid")%in%names(x)))stop("Identifier in aggregate output")
  fwrite(x,file.path(out,f),na="")
}
minimum<-as.integer(Sys.getenv("NCO_MIN_CELL_N","10"))
if(!is.finite(minimum)||minimum<2L)stop("NCO_MIN_CELL_N must be >=2")
c<-as.data.table(read_parquet(file.path(inp,"pneumonia_ccw_cohort.parquet")))
cw<-as.data.table(read_parquet(file.path(inp,"sampled_hospital_block_crosswalk.parquet")))
assert_mapped_sites(c$site)
if(anyDuplicated(c$patient_id)||anyDuplicated(c$hospital_block_id))stop("Expected one block per patient")
ds<-open_dataset(expanded)
needed<-c("patient_id","hospitalization_id","admission_dttm","discharge_dttm")
if(!all(needed%in%ds$schema$names))stop("Expanded hospitalization table lacks required columns")
optional<-intersect(c("admission_type_category","discharge_category"),ds$schema$names)
ids<-as.character(c$patient_id)
hosp<-ds |> filter(.data$patient_id%in%ids) |> select(all_of(c(needed,optional))) |> collect() |> as.data.table()
panel<-build_historical_nco_panel(c,cw,hosp,coverage);c<-panel$cohort
write(panel$definitions,"nco_definitions.csv");write(panel$source_summary,"historical_source_audit.csv")
if(ncol(panel$mapping))write(panel$mapping,"historical_category_mapping_audit.csv")
write(data.table(nco="prior_icu_admission",status="Not constructed: expanded hospitalization alone does not establish ICU location history"),"unavailable_candidates.csv")
write_parquet(c,file.path(sensitive,"historical_nco_panel.parquet"))
write_parquet(panel$history,file.path(sensitive,"qualifying_historical_episodes.parquet"))
definitions<-panel$definitions
availability<-rbindlist(lapply(definitions$nco,function(nm)c[,.(nco=nm,patients=.N,
  observed=sum(!is.na(get(nm))),positive=sum(get(nm)==1L,na.rm=TRUE),unknown=sum(is.na(get(nm))),
  unverified_coverage=sum(is.na(coverage_verified))),by=site]))
write(availability,"nco_availability_by_site.csv")
weights<-readRDS(file.path(inp,"weight_diagnostics_patient_level.rds"));cfg<-weights$config
long<-as.data.frame(read_parquet(file.path(inp,"pneumonia_ccw_long.parquet")))
assert_mapped_sites(long$site)
p<-g6_prepare_data(long,"site","hospital_block_id","m","event_m","censor_m","init_m",NULL,
                   cfg$baseline,cfg$tv,cfg$horizon,cfg$lag_initial)
map<-unique(as.data.table(p$data)[,.(uid=as.character(.uid),hospital_block_id=as.character(.id))])
if(anyDuplicated(map$uid))stop("Nonunique saved unit map")
map<-merge(map,c[,c("hospital_block_id",definitions$nco),with=FALSE],by="hospital_block_id")
if(!all(unlist(lapply(weights$site_components,`[[`,"uid"))%in%map$uid))stop("Saved CCW units do not match cohort")
contrasts<-list()
for(m in seq_len(cfg$tau))for(version in c("raw","clipped")) {
  arms<-rbindlist(lapply(weights$site_components,function(z) {
    rows<-map[match(z$uid,map$uid)]
    w1<-if(version=="raw")z$raw_W1[,m] else z$W1[,m]
    w0<-if(version=="raw")z$raw_W0[,m] else z$W0[,m]
    rows[,`:=`(w1=w1*z$R1[,m],w0=w0*z$R0[,m])]
    rows
  }))
  for(nm in definitions$nco)contrasts[[length(contrasts)+1L]]<-cbind(
    data.table(nco=nm,interval=m,decision_hours=(m-1L)*6,version=version),
    historical_log_pr(arms[[nm]],arms$w1,arms$w0))
}
contrasts<-rbindlist(contrasts)
write(contrasts,"nco_log_prevalence_ratios_by_interval.csv")
final<-merge(contrasts[interval==cfg$tau&version=="clipped"],definitions,by="nco",sort=FALSE)
final[,exclusion_reason:=fcase(!field_available,"Historical field unavailable",
  !is.finite(log_pr)|!is.finite(se_log_pr)|se_log_pr<=0,"Nonestimable log prevalence ratio/SE",
  positive_early<minimum|positive_control<minimum|negative_early<minimum|negative_control<minimum,
  "Sparse positive/negative observations in one or both arms",default="")]
# Exact duplicate patient-level controls do not create additional evidence.
seen<-character()
for(i in seq_len(nrow(final)))if(final$exclusion_reason[i]=="") {
  signature<-paste(ifelse(is.na(c[[final$nco[i]]]),"?",c[[final$nco[i]]]),collapse="")
  if(signature%in%seen)final$exclusion_reason[i]<-"Duplicates another included patient-level NCO"
  else seen<-c(seen,signature)
}
final[,included:=exclusion_reason==""]
write(final,"nco_calibration_inputs.csv")
ymat<-as.matrix(c[,definitions$nco,with=FALSE])
corr<-suppressWarnings(cor(ymat,use="pairwise.complete.obs"))
write(data.table(nco=rownames(corr),as.data.table(corr)),"nco_patient_level_correlations.csv")
primary<-fread(file.path(run,"aggregate","analysis","five_method_results.csv"))[method_key=="fed_ccw"]
if(nrow(primary)!=1L||!is.finite(primary$risk_ratio)||primary$risk_ratio<=0||
   !is.finite(primary$log_risk_ratio_se)||primary$log_risk_ratio_se<=0)stop("Primary CCW RR/SE unavailable")
theta<-log(primary$risk_ratio);primary_se<-primary$log_risk_ratio_se
inputs<-final[included==TRUE]
fits<-list();fit_rows<-list();calibrated<-list();sensitivity<-list()
fit_one<-function(d,model)tryCatch(
  if(model=="NN")fit_nn_calibration(d$log_pr,d$se_log_pr) else fit_mnn_calibration(d$log_pr,d$se_log_pr),error=identity)
for(model in c("NN","MNN")) {
  fit<-fit_one(inputs,model)
  if(inherits(fit,"error")) {
    fit_rows[[model]]<-data.table(model=model,n_controls=nrow(inputs),converged=FALSE,status=conditionMessage(fit))
  } else {
    fits[[model]]<-fit
    fit_rows[[model]]<-data.table(model=model,n_controls=nrow(inputs),converged=fit$converged,
      mu_bias=fit$mu_valid,sigma_bias=fit$sigma_valid,pi_valid=fit$pi_valid,
      mu_other=if(model=="MNN")fit$mu_invalid else NA_real_,sigma_other=if(model=="MNN")fit$sigma_invalid else NA_real_,
      log_likelihood=fit$logLik,small_panel=nrow(inputs)<20L,
      weak_component_separation=if(model=="MNN")abs(fit$mu_valid-fit$mu_invalid)<
        sqrt(median(inputs$se_log_pr^2)+(fit$sigma_valid^2+fit$sigma_invalid^2)/2) else NA,
      mixture_boundary=model=="MNN"&&(fit$pi_valid<.51||fit$pi_valid>.99),status="Exploratory fit only")
    z<-calibrate_primary(theta,primary_se,fit,model)
    z$rr<-exp(z$estimate_logRR);z$rr_lower<-exp(z$ci_low);z$rr_upper<-exp(z$ci_high)
    z$p_value<-NULL # No calibrated hypothesis-test claim for this small correlated panel.
    z$calibration_supported<-FALSE
    z$interpretation<-"Assumption-dependent sensitivity result; historical prevalence bias transport not established"
    z$interval_type<-"Plug-in normal interval: primary sampling SE plus between-control bias SD; calibration-parameter uncertainty and dependence omitted"
    calibrated[[model]]<-as.data.table(z)
  }
  for(omit in inputs$nco) {
    loo<-fit_one(inputs[nco!=omit],model)
    sensitivity[[length(sensitivity)+1L]]<-data.table(model=model,omitted=omit,n_controls=nrow(inputs)-1L,
      mu_bias=if(inherits(loo,"error"))NA_real_ else loo$mu_valid,
      sigma_bias=if(inherits(loo,"error"))NA_real_ else loo$sigma_valid,
      rr=if(inherits(loo,"error"))NA_real_ else exp(theta-loo$mu_valid),
      status=if(inherits(loo,"error"))conditionMessage(loo) else "Exploratory")
  }
}
write(rbindlist(fit_rows,fill=TRUE),"nn_mnn_fit_diagnostics.csv")
write(data.table(method="Uncalibrated federated CCW",estimate_logRR=theta,se_logRR=primary_se,
  rr=primary$risk_ratio,rr_lower=primary$risk_ratio_lower,rr_upper=primary$risk_ratio_upper),"primary_ccw_reference.csv")
if(length(calibrated)) {
  write(rbindlist(calibrated,fill=TRUE),"exploratory_calibrated_ccw_rr.csv")
} else write(data.table(status="No usable NN/MNN calibration fit; see fit diagnostics"),"exploratory_calibrated_ccw_rr.csv")
if(length(sensitivity))write(rbindlist(sensitivity),"leave_one_control_out.csv")
if(!is.null(fits$MNN)) {
  z<-fits$MNN
  lv<-log(z$pi_valid)+dnorm(inputs$log_pr,z$mu_valid,sqrt(inputs$se_log_pr^2+z$sigma_valid^2),log=TRUE)
  li<-log1p(-z$pi_valid)+dnorm(inputs$log_pr,z$mu_invalid,sqrt(inputs$se_log_pr^2+z$sigma_invalid^2),log=TRUE)
  write(data.table(nco=inputs$nco,probability_majority_component=plogis(lv-li),
    interpretation="Model membership only; does not validate the NCO"),"nco_mnn_membership.csv")
}
saveRDS(list(fits=fits,inputs=inputs,primary=primary),file.path(sensitive,"calibration_fits.rds"))
if(length(calibrated)) {
  comparison<-rbindlist(list(data.table(method="Uncalibrated CCW",rr=primary$risk_ratio,
    rr_lower=primary$risk_ratio_lower,rr_upper=primary$risk_ratio_upper),
    rbindlist(calibrated)[,.(method=paste0(method," exploratory calibration"),rr,rr_lower,rr_upper)]))
  pl<-ggplot(comparison,aes(rr,method))+geom_vline(xintercept=1,linetype=3)+
    geom_errorbar(aes(xmin=rr_lower,xmax=rr_upper),orientation="y",width=.15)+geom_point()+scale_x_log10()+
    theme_minimal()+labs(title="CCW RR: assumption-dependent NCO calibration",x="Risk ratio",y=NULL,
      caption="Calibrated bars are plug-in sensitivity intervals; historical-to-mortality bias transport is unverified.")
  ggsave(file.path(out,"exploratory_ccw_rr_forest.png"),pl,width=11,height=5,dpi=180)
}
pl<-ggplot(final[is.finite(prevalence_ratio)],aes(prevalence_ratio,reorder(label,prevalence_ratio),color=included))+
  geom_vline(xintercept=1,linetype=3)+geom_errorbar(aes(xmin=pr_lower,xmax=pr_upper),orientation="y",width=.15)+
  geom_point()+scale_x_log10()+theme_minimal()+labs(title="Historical NCOs: final grace-period CCW risk set",
    subtitle="Historical prevalence ratios, not prospective event risk ratios",x="Weighted prevalence ratio, 95% fixed-weight interval",y=NULL)
ggsave(file.path(out,"historical_nco_forest.png"),pl,width=11,height=6,dpi=180)
if(length(fits)) {
  span<-range(c(inputs$log_pr,unlist(lapply(fits,function(z)c(z$mu_valid-4*max(z$sigma_valid,.05),z$mu_valid+4*max(z$sigma_valid,.05))))))
  x<-seq(span[1],span[2],length.out=600)
  density<-rbindlist(lapply(names(fits),function(nm) {
    z<-fits[[nm]];d<-data.table(model=nm,log_bias=x,density=dnorm(x,z$mu_valid,max(z$sigma_valid,.001)),component="Assumed valid bias distribution")
    if(nm=="MNN")d<-rbind(d,data.table(model=nm,log_bias=x,density=z$pi_valid*dnorm(x,z$mu_valid,max(z$sigma_valid,.001))+
      (1-z$pi_valid)*dnorm(x,z$mu_invalid,max(z$sigma_invalid,.001)),component="Full mixture"))
    d
  }))
  write(density,"fitted_bias_density_plot_data.csv")
  pl<-ggplot(density,aes(log_bias,density,color=component))+geom_line()+facet_wrap(~model,scales="free_y")+
    geom_vline(xintercept=0,linetype=3)+geom_rug(data=inputs,aes(x=log_pr),inherit.aes=FALSE)+theme_minimal()+
    labs(title="Exploratory historical NCO systematic-error distributions",x="Log prevalence-ratio bias",y="Density",
      caption="Rug: noisy control estimates. Near-zero SD rendered with floor 0.001. Correlated small panel; mortality-bias transport unverified.")
  ggsave(file.path(out,"nn_mnn_bias_distributions.png"),pl,width=12,height=6,dpi=180)
}
writeLines(c(
  "This is exploratory empirical bias calibration, not identified causal correction or validated CI calibration.",
  "Controls: any/recurrent admission, stay >=7 days, ED-origin admission, SNF/rehab discharge in the preceding 365 days.",
  "Only completed episodes starting in the window count; index-connected and future episodes excluded; <=6-hour gaps merged.",
  "ED uses first component admission category; SNF/rehab uses last component discharge category. Inspect source-category mapping audit.",
  "Prior ICU admission is unavailable without historical ADT. Missing optional fields are not converted to negative controls.",
  paste("System coverage start:",coverage,"; absent date means capture is unverified, not complete."),
  "Reuses unchanged 07 confounders and weights. Fits final-grace-interval CLIPPED historical prevalence log ratios only; raw/all-grace-interval contrasts are exported for sensitivity.",
  "Historical prevalence ratios conditional on adherence/survival are not day-28 event RRs. Selection can generate associations without residual confounding.",
  "Transport assumes primary log-RR systematic error is exchangeable with valid historical log-prevalence-ratio errors and invariant with true effect size.",
  "NN assumes all included controls are valid. MNN assumes a strict majority form the valid component, not the component nearest zero.",
  "NN requires >=3; exploratory MNN >=4 distinct finite nonsparse controls. These are computational guards, not sample-size adequacy claims.",
  "Five overlapping hospitalization controls do not provide five independent draws; mixtures and dispersions can be unstable. No extra windows are manufactured to inflate panel size.",
  "Failed/nonconverged MNN is reported as failed, never relabeled NN. Exact duplicate controls excluded; pairwise correlations and leave-one-out fits exported.",
  "Correction: calibrated logRR=primary logRR-mu_valid; RR=exp(calibrated logRR). Plug-in variance=primary SE^2+sigma_valid^2.",
  "Intervals omit fitted-weight uncertainty, calibration-parameter uncertainty, cross-control/primary dependence, and clinical transport uncertainty. No formal calibrated p-values.",
  "calibration_supported is FALSE pending clinical justification, sufficient controls and validation. Inspect results as sensitivity analyses; do not select controls to obtain benefit.",
  "No RD/RMST calibration. No changes to 07/08/09 outputs. Original patient-level data stay under sensitive/."
),file.path(out,"METHOD_NOTES.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"))
print(rbindlist(fit_rows,fill=TRUE));cat("\nResults:",out,"\n")
