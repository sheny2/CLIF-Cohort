#!/usr/bin/env Rscript
# Two pooled weighted Cox analyses using script-07 inputs; existing files are
# not changed. Usage:
# Rscript 09_ccw_and_landmark_cox.R /path/to/07/sensitive/ccw_analysis [NEW_OUTPUT_DIR]
# Default input: pipeline_output/sensitive/ccw_analysis.
# The saved grace period, six-hour outcome grid, 1st/99th percentile clipping,
# CCW risk sets, and original landmark propensity model are preserved.
# Cox models adjust for strategy and FIXED time-zero covariates, stratify by
# hospital, and cluster robust SEs on the original subject across clones/rows.
# Current/lagged longitudinal covariates remain in weight models only.
# Important: the preserved landmark propensity model itself uses end-of-grace
# covariates that can be post-treatment. This script changes the outcome model,
# not that design choice; it does not cure survivor selection or confounding.

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(survival); library(ggplot2)
})
a <- commandArgs(TRUE)
sf <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(normalizePath(gsub("~+~"," ",sf,fixed=TRUE)))
source(file.path(root,"06_generalized_five_methods.R"))
input <- if(length(a)) a[1] else file.path(root,"pipeline_output/sensitive/ccw_analysis")
out <- if(length(a)>=2) a[2] else file.path(dirname(input),paste0("cox_analysis_",format(Sys.time(),"%Y%m%d_%H%M%S")))
if(normalizePath(input,mustWork=FALSE)==normalizePath(out,mustWork=FALSE)) stop("Use a separate output directory.")
agg <- file.path(out,"aggregate"); sensitive <- file.path(out,"sensitive")
dir.create(agg,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
write_agg <- function(x,f) {
  if(any(c("patient_id","hospital_block_id",".uid","subject")%in%names(x))) stop("Identifier in aggregate output")
  fwrite(x,file.path(agg,f),na="")
}
design <- readRDS(file.path(input,"analysis_design.rds"))
long <- as.data.frame(read_parquet(file.path(input,"pneumonia_ccw_long.parquet")))
M <- design$horizon_intervals; tau <- design$tau_intervals; step <- design$step_hours
prep <- g6_prepare_data(long,"site","hospital_block_id","m","event_m","censor_m",
  "init_m",NULL,design$baseline_covariates,design$time_varying_covariates,M,0)
d <- prep$data;rm(long);invisible(gc())
terms <- list(initiation_denominator=c(".time",design$baseline_covariates,
  design$time_varying_covariates,prep$lag_names),initiation_numerator=".time",
  censoring_denominator=c(".time",design$baseline_covariates,".trt"),censoring_numerator=".time")
learner <- g6_default_binary_learner()
cfg <- list(tau=tau,horizon=M,trunc=c(.01,.99),terms=terms,
            learners=list(initiation=learner,censoring=learner))
warnings <- data.table(method=character(),stage=character(),message=character())
capture <- function(expr,method,stage) withCallingHandlers(expr,warning=function(w) {
  warnings <<- rbind(warnings,data.table(method=method,stage=stage,message=conditionMessage(w)))
  invokeRestart("muffleWarning")
})

# Fix all extra outcome-model covariates at m=1 (time zero). Do not carry
# evolving post-initiation measurements into a conventional adjusted Cox model.
base <- as.data.table(d[d$.time==1,,drop=FALSE])
base[, .uid := as.character(.uid)]
tv <- design$time_varying_covariates
setnames(base,tv,paste0("t0_",tv))
candidate <- c(design$baseline_covariates,paste0("t0_",tv))
base <- base[,c(".uid",candidate),with=FALSE]
retained <- character(); audit <- list()
for(v in candidate) {
  duplicate <- retained[vapply(retained,function(u) identical(base[[u]],base[[v]]),logical(1))]
  reason <- if(uniqueN(base[[v]])<=1L) "omitted: constant" else
    if(length(duplicate)) paste0("omitted: identical to ",duplicate[1]) else "included"
  if(reason=="included") retained <- c(retained,v)
  audit[[v]] <- data.table(variable=v,status=reason,measurement_window="fixed at time zero")
}
write_agg(rbindlist(audit),"cox_covariate_specification.csv")
common <- g6_common_init_numerator(d,M)
site_data <- split(d,d$.site)
land <- lapply(site_data,function(z) z[z$.time==tau & z$.event==0L & z$.censor==0L,,drop=FALSE])
land_n <- sum(vapply(land,nrow,integer(1)))
if(!land_n) stop("No landmark survivors")
p_num <- g6_clamp_prob(sum(vapply(land,function(z) sum(z$.S<=tau),numeric(1)))/land_n)

make_frame <- function(ds,x,uid) {
  row <- match(as.character(ds$.uid),uid)
  cell <- cbind(row,ds$.time)
  rbindlist(lapply(0:1,function(arm) {
    R <- if(arm==1) x$R1 else x$R0
    W <- if(arm==1) x$W1 else x$W0
    D <- if(arm==1) x$D1 else x$D0
    keep <- !is.na(row) & R[cell]
    data.table(.uid=as.character(ds$.uid[keep]),site=as.character(ds$.site[keep]),
      strategy=arm,start=(ds$.time[keep]-1)*step/24,stop=ds$.time[keep]*step/24,
      event=as.integer(D[cell][keep]),weight=W[cell][keep])
  }))
}
ccw_frames <- landmark_frames <- components <- list()
for(site in names(site_data)) {
  message("Reconstructing preserved weights and risk sets: ",site)
  ds <- site_data[[site]]
  x <- capture(g6_local_components(ds,cfg,common),"CCW",site)
  components[[site]] <- x
  ccw_frames[[site]] <- make_frame(ds,x,x$uid)
  if(nrow(land[[site]])) {
    ls <- capture(g6_landmark_local_stats(ds,cfg,p_num),"Landmark",site)
    landmark_frames[[site]] <- make_frame(ds,ls$local,as.character(land[[site]]$.uid))
  }
}
# Verify reconstructed weights/risk sets against the saved original estimates.
saved_file <- file.path(input,"five_method_results.rds")
if(file.exists(saved_file)) {
  saved <- readRDS(saved_file)$fits
  ccw_check <- g6_central_ccw(lapply(components,g6_local_stats))
  lm_check <- capture(g6_run_landmark_ipw_general(d,cfg),"Landmark","reproduction")
  check <- rbindlist(lapply(c("fed_ccw","landmark_ipw"),function(method) {
    current <- if(method=="fed_ccw") ccw_check else lm_check
    if(is.null(saved[[method]])) return(NULL)
    data.table(method=method,quantity=c("RD","RR"),
      original=c(saved[[method]]$RD,saved[[method]]$RR),reconstructed=c(current$RD,current$RR))
  }))
  if(nrow(check)) {
    check[,`:=`(absolute_difference=abs(original-reconstructed),passed=abs(original-reconstructed)<1e-8)]
    write_agg(check,"original_07_reproduction_check.csv")
    if(!all(check$passed)) stop("07 reproduction check failed; investigate before fitting Cox models.")
  }
}
rm(components);invisible(gc())
frames <- list(CCW=rbindlist(ccw_frames),Landmark=rbindlist(landmark_frames))
rm(ccw_frames,landmark_frames)
fits <- results <- counts <- list()
for(method in names(frames)) {
  message("Fitting weighted adjusted Cox model: ",method)
  z <- merge(frames[[method]],base,by=".uid",all.x=TRUE,sort=FALSE)
  if(any(!is.finite(z$weight)|z$weight<=0) || any(z$stop<=z$start)) stop("Invalid Cox rows")
  # .uid is the original patient/block key, shared by both clones.
  z[,subject:=.uid]
  active <- retained[vapply(retained,function(v) uniqueN(z[[v]])>1L,logical(1))]
  form <- reformulate(c("strategy",active,"strata(site)"),response="Surv(start, stop, event)")
  fit <- capture(coxph(form,data=z,weights=weight,cluster=subject,robust=TRUE,
    ties="breslow",x=TRUE,y=TRUE,model=TRUE),method,"cox_fit")
  if(!is.finite(coef(fit)["strategy"])) stop("Strategy coefficient not estimable for ",method)
  sm <- summary(fit)$coefficients
  beta <- coef(fit)["strategy"]; se <- sqrt(vcov(fit)["strategy","strategy"])
  results[[method]] <- data.table(method=paste0(method," + weighted adjusted Cox"),
    patients=uniqueN(z$subject),grace_hours=tau*step,hazard_ratio=unname(exp(beta)),
    ci_lower=unname(exp(beta-1.96*se)),ci_upper=unname(exp(beta+1.96*se)),
    log_hr_robust_se=unname(se),p_value=unname(2*pnorm(-abs(beta/se))),
    outcome_covariates=length(active),ties="breslow",time_scale="days since original time zero; six-hour interval endpoints")
  coef_table <- data.table(term=rownames(sm),as.data.table(sm))
  write_agg(coef_table,paste0(tolower(method),"_cox_coefficients.csv"))
  counts[[method]] <- z[,.(patients=uniqueN(subject),patient_intervals=.N,
    events=sum(event),min_weight=min(weight),max_weight=max(weight)),by=.(site,strategy)][,method:=method]
  ph <- tryCatch(capture(cox.zph(fit),method,"PH_assessment"),error=function(e)e)
  if(inherits(ph,"error")) {
    warnings <- rbind(warnings,data.table(method=method,stage="PH_assessment",message=conditionMessage(ph)))
  } else {
    write_agg(data.table(term=rownames(ph$table),as.data.table(ph$table)),
              paste0(tolower(method),"_ph_assumption_check.csv"))
    if("strategy" %in% colnames(ph$y)) {
      png(file.path(agg,paste0(tolower(method),"_strategy_schoenfeld.png")),width=1100,height=750,res=130)
      plot(ph,var=match("strategy",colnames(ph$y)),resid=FALSE)
      title(main=paste(method,"strategy: scaled Schoenfeld diagnostic"))
      dev.off()
    }
  }
  fits[[method]] <- fit
}
result <- rbindlist(results)
write_agg(result,"weighted_cox_hazard_ratios.csv")
write_agg(rbindlist(counts),"cox_risk_set_summary.csv")
write_agg(unique(warnings),"model_warnings.csv")
saveRDS(list(fits=fits,design=design),file.path(sensitive,"weighted_cox_models.rds"))
p <- ggplot(result,aes(hazard_ratio,method))+
  geom_vline(xintercept=1,linetype=2,color="grey50")+
  geom_segment(aes(x=ci_lower,xend=ci_upper,yend=method))+geom_point(size=3)+scale_x_log10()+
  labs(title="Weighted Cox strategy hazard ratios",x="Hazard ratio (95% robust CI)",y=NULL)+theme_minimal(base_size=12)
ggsave(file.path(agg,"weighted_cox_hazard_ratios.png"),p,width=10,height=4.5,dpi=180)
writeLines(c(
  "Two methods: pooled weighted adjusted Cox following CCW, and pooled weighted adjusted landmark Cox.",
  "Weight models, numerator construction, clipping, risk sets and six-hour outcome grid are preserved from 07.",
  "Events are placed at their six-hour interval endpoints; this is not an exact-event-time Cox analysis.",
  "CCW strategy is the assigned clone arm, not observed time-varying treatment. Clone rows share original-subject clustering.",
  "Landmark begins follow-up after the grace period among survivors; treatment group is fixed by initiation through the landmark.",
  "Landmark controls may initiate later, matching 07; its contrast is not identical to CCW's never-initiate control strategy.",
  "Both Cox models adjust for fixed time-zero demographics/comorbidities and first-interval clinical covariates; hospital baseline hazards are stratified.",
  "Longitudinal current/lagged covariates remain in the weights, not as evolving covariates in the outcome model.",
  "The baseline-adjusted HR is conditional on included baseline covariates and site. It is not a risk ratio or RMST contrast, or a purely marginal HR.",
  "Original landmark weights use end-of-grace measurements that may follow treatment. That limitation is preserved and not corrected here.",
  "Robust SEs cluster on original subjects, including both clones; estimated weight-model and truncation uncertainty is not propagated.",
  "Schoenfeld diagnostics are exploratory usual cox.zph diagnostics, not bootstrap tests accounting for estimated weights/cloning.",
  "Changing outcome models does not resolve residual confounding, positivity issues, or guarantee a beneficial estimate.",
  "Aggregate summaries contain no patient identifiers. Keep sensitive fitted model objects on the server."
),file.path(out,"START_HERE.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
writeLines("Complete",file.path(out,"COMPLETE.txt"))
print(result)
message("Cox results: ",normalizePath(out))
