# Helpers for 08 only; no change to the running 07 pipeline.
local_ccw_prepare <- function(long,design) {
 x <- data.table::copy(data.table::as.data.table(long))
 tau <- design$tau_intervals
 decision <- x$m<=tau & x$m<=x$S & x$m<=x$adherence_end_m
 audit <- list()
 numeric <- sub("_model$","",grep("_model$",design$time_varying_covariates,value=TRUE))
 for(v in numeric) {
  if(!v%in%names(x))stop("Missing screened raw covariate for local imputation: ",v)
  z <- screen_measurements(x[[v]],v);observed<-is.finite(z)
  pool <- z[decision & observed]
  limits <- c(NA_real_,NA_real_) # No additional winsorization
  if(length(pool)) {
   fill <- median(z[decision&observed]);z[!observed] <- fill
  } else {fill<-0;z[]<-0}
  if(v%in%c("wbc","glucose"))z<-log1p(z)
  center <- mean(z[decision]);scale<-sd(z[decision])
  if(!is.finite(scale)||scale<=0)scale<-1
  x[,(paste0(v,"_missing")):=as.integer(!observed)]
  x[,(paste0(v,"_model")):=if(v=="nonresp_sofa_proxy")z else (z-center)/scale]
  audit[[v]] <- data.table::data.table(variable=v,decision_rows=sum(decision),
   observed_decision_rows=length(pool),imputation_value=fill,winsor_lower=limits[1],winsor_upper=limits[2],
   transformed_center=center,transformed_scale=scale)
 }
 # Local constants cannot distinguish initiation decisions; omit their terms.
 informative <- function(v)length(unique(x[[v]][decision]))>1L
 baseline <- design$baseline_covariates[vapply(design$baseline_covariates,informative,logical(1))]
 tv <- design$time_varying_covariates[vapply(design$time_varying_covariates,informative,logical(1))]
 prepared <- g6_prepare_data(data=as.data.frame(x),site_col="site",id_col="hospital_block_id",time_col="m",
   event_col="event_m",censor_col="censor_m",initiation_col="init_m",treatment_col=NULL,
   baseline_covariates=baseline,time_varying_covariates=tv,horizon=design$horizon_intervals,lag_initial=0)
 terms <- list(initiation_denominator=c(".time",baseline,tv,prepared$lag_names),initiation_numerator=".time",
   censoring_denominator=c(".time",baseline,".trt"),censoring_numerator=".time")
 cfg <- list(tau=tau,horizon=design$horizon_intervals,baseline=baseline,tv=tv,lag_names=prepared$lag_names,
   terms=terms,learners=list(initiation=g6_default_binary_learner(),censoring=g6_default_binary_learner()),
   trunc=c(.01,.99),lag_initial=0,discharge_rule=G6_DISCHARGE_RULE)
 list(data=prepared$data,cfg=cfg,processing=data.table::rbindlist(audit,fill=TRUE))
}

meta_contrasts <- function(site_results) {
 x <- data.table::as.data.table(site_results)
 data.table::rbindlist(list(
  x[,.(site,estimand="risk_difference",estimate=risk_difference,se=risk_difference_se,site_valid,site_status)],
  x[,.(site,estimand="log_risk_ratio",estimate=log(risk_ratio),se=log_risk_ratio_se,site_valid,site_status)],
  x[,.(site,estimand="log_odds_ratio",estimate=log(odds_ratio),se=log_odds_ratio_se,site_valid,site_status)],
  x[,.(site,estimand="rmst_difference_days",estimate=rmst_difference_days,se=rmst_difference_se_days,site_valid,site_status)]))
}

pool_local_contrasts <- function(contrasts) {
 summaries<-weights<-list();models<-list();warnings<-character()
 for(e in unique(contrasts$estimand)) {
  z<-data.table::copy(contrasts[estimand==e])
  z[,included:=site_valid & is.finite(estimate)&is.finite(se)&se>0]
  z[,exclusion_reason:=data.table::fcase(!site_valid,site_status,!is.finite(estimate),"Nonfinite contrast (e.g. zero arm risk for a ratio)",
    !is.finite(se)|se<=0,"Nonpositive/nonfinite variance",default="")]
  usable<-z[included==TRUE]
  ratio<-e%in%c("log_risk_ratio","log_odds_ratio")
  for(method in c("REML_HK","fixed_effect")) {
   model<-NULL;failure<-""
   if(nrow(usable)>=2) {
    model<-tryCatch(withCallingHandlers(
     metafor::rma.uni(yi=usable$estimate,sei=usable$se,slab=usable$site,
       method=if(method=="REML_HK")"REML" else "FE",test=if(method=="REML_HK")"knha" else "z"),
     warning=function(w){warnings<<-c(warnings,paste(e,method,conditionMessage(w)));invokeRestart("muffleWarning")}),
     error=function(err){failure<<-conditionMessage(err);NULL})
   } else failure<-"Fewer than two sites with valid estimates and positive variances"
   estimate<-se<-lo<-hi<-pi_lo<-pi_hi<-tau2<-i2<-q<-qp<-NA_real_
   if(!is.null(model)) {
    estimate<-as.numeric(model$b);se<-as.numeric(model$se);lo<-model$ci.lb;hi<-model$ci.ub
    tau2<-model$tau2;i2<-model$I2;q<-model$QE;qp<-model$QEp
    if(method=="REML_HK") {
     pred<-tryCatch(predict(model),error=function(err)NULL)
     if(!is.null(pred)){pi_lo<-pred$pi.lb;pi_hi<-pred$pi.ub}
    }
    models[[paste(e,method)]]<-model
   }
   summaries[[paste(e,method)]]<-data.table::data.table(estimand=e,meta_method=method,k=nrow(usable),
    analysis_scale_estimate=estimate,analysis_scale_se=se,analysis_scale_lower=lo,analysis_scale_upper=hi,
    estimate=if(ratio)exp(estimate) else estimate,lower=if(ratio)exp(lo) else lo,upper=if(ratio)exp(hi) else hi,
    prediction_lower=if(ratio)exp(pi_lo) else pi_lo,prediction_upper=if(ratio)exp(pi_hi) else pi_hi,
    tau_squared=tau2,I_squared_percent=i2,Q=q,Q_pvalue=qp,status=if(is.null(model))failure else "ok")
   z[,`:=`(meta_method=method,weight_percent=NA_real_)]
   if(!is.null(model))z[included==TRUE,weight_percent:=as.numeric(weights(model))]
   weights[[paste(e,method)]]<-data.table::copy(z)
  }
 }
 list(summary=data.table::rbindlist(summaries),weights=data.table::rbindlist(weights),models=models,warnings=unique(warnings))
}
