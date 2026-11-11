#!/usr/bin/env Rscript
# Early versus later/never initiation; same primary cohort/covariate preparation.
# Run after 01--05; no run of 07/08/09/0701/0901 required.
# Usage: Rscript 0902_early_vs_later_cox.R CLIF_DATA_DIR [COHORT_DIR] [TREATMENT_DIR] [NEW_OUTPUT_DIR] [DATA_CUTOFF_DATE]
# Default TREATMENT_GRACE_HOURS=24; original missing-data rules are preserved.
# Requires 06_generalized_five_methods.R and early_initiation_comparison_helpers.R.
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2);library(survival)})
a <- commandArgs(TRUE)
sf <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(normalizePath(gsub("~+~"," ",sf,fixed=TRUE)))
source(file.path(root,"06_generalized_five_methods.R"))
source(file.path(root,"early_initiation_comparison_helpers.R"))
arg <- function(i,default) if(length(a)>=i && nzchar(a[i])) a[i] else default
data_dir <- arg(1,Sys.getenv("CLIF_DATA_DIR","/projects/emulate_icu/rclif_2_1"))
cohort_dir <- arg(2,Sys.getenv("COHORT_OUTPUT_DIR",file.path(root,"pipeline_output/sensitive/cohort")))
treatment_dir <- arg(3,Sys.getenv("STEROID_TREATMENT_DIR",file.path(root,"pipeline_output/sensitive/steroid_treatment")))
out <- arg(4,file.path(root,"pipeline_output",paste0("0902_early_vs_later_",format(Sys.time(),"%Y%m%d_%H%M%S"))))
cutoff_text <- arg(5,Sys.getenv("CLIF_DATA_CUTOFF_DATE","2026-10-31"))
grace_hours <- suppressWarnings(as.numeric(Sys.getenv("TREATMENT_GRACE_HOURS","24")))
if(!is.finite(grace_hours) || grace_hours<6 || grace_hours>=672 || grace_hours%%6!=0) stop("Grace hours must be a multiple of 6 below 672")
if(normalizePath(out,mustWork=FALSE)%in%vapply(c(data_dir,cohort_dir,treatment_dir),normalizePath,character(1),mustWork=FALSE)) stop("Use a separate output folder")
agg <- shareable_dir <- file.path(out,"aggregate")
sensitive <- sensitive_dir <- file.path(out,"sensitive")
dir.create(agg,recursive=TRUE,showWarnings=FALSE);dir.create(sensitive,recursive=TRUE,showWarnings=FALSE)
write_agg <- write_shareable_csv <- function(x,f) {
 if(any(c("patient_id","hospitalization_id","hospital_block_id",".uid","subject")%in%names(x))) stop("Identifier in aggregate output")
 fwrite(x,file.path(agg,f),na="")
}
write_sensitive_csv <- function(x,f) fwrite(x,file.path(sensitive,f),na="")
state <- prepare_early_comparison_inputs(data_dir,cohort_dir,treatment_dir,sensitive,agg,grace_hours,cutoff_text,0L)
invisible(list2env(state,envir=environment()));rm(state)
site_summary <- cohort[,.(patients=.N,deaths_28d=sum(delta)),by=site]
cohort[, initiation_group := fcase(S<=tau,"Within grace period",S<=M,"Later initiation",default="No observed initiation through day 28")]
write_agg(cohort[,.(patients=.N,deaths_28d=sum(delta)),by=.(site,initiation_group)],"initiation_groups_by_site.csv")
write_agg(data.table(parameter=c("strategy_1","strategy_0","grace_hours","horizon_days","missing_data_rules"),
 value=c(design$strategy_1,design$strategy_0,as.character(grace_hours),as.character(horizon_days),"original primary analysis; not complete-case sensitivity")),"early_initiation_design.csv")
writeLines(c("Initiation within grace versus no initiation within grace; later initiation is allowed in controls.",
 "Control clones are censored only for initiation within grace. Later initiation does not censor them.",
 "Treatment/adherence weights freeze at grace end; natural-censoring weights, if applicable, can continue.",
 "All other preparation, six-hour grid/event boundary conventions and 1st/99th percentile weight truncation follow the original methods.",
 "Original missingness indicators and imputation are retained; this is not the 0701/0901 complete-case sensitivity.",
 "Landmark already classified early versus later/never, so its definition and results should remain unchanged.",
 "Observed early-treatment labels do not replace CCW cloning: deaths during grace remain handled by the inherited clone risk sets.",
 "No diagnosis, eligibility or pre-existing steroid-use definitions have been changed."
 ),file.path(out,"EARLY_INITIATION_RULES.txt"))
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

# Align each Cox row to its own decision interval, rather than the m=1 row.
# The clinical values were selected at/before the START of that interval by
# prepare_cox_inputs(), then processed using the existing imputation rules.
tv <- design$time_varying_covariates
candidate <- c(design$baseline_covariates,tv)
covariates <- as.data.table(d)[,c(".uid",".time",candidate),with=FALSE]
covariates[, .uid := as.character(.uid)]
setnames(covariates,".time","interval")
if(anyDuplicated(covariates[,.(.uid,interval)])) stop("Duplicate patient-interval covariates")
retained <- character(); audit <- list()
for(v in candidate) {
  duplicate <- retained[vapply(retained,function(u) identical(covariates[[u]],covariates[[v]]),logical(1))]
  reason <- if(uniqueN(covariates[[v]])<=1L) "omitted: constant" else
    if(length(duplicate)) paste0("omitted: identical to ",duplicate[1]) else "included"
  if(reason=="included") retained <- c(retained,v)
  audit[[v]] <- data.table(variable=v,status=reason,measurement_window=
    if(v %in% tv) "current interval; aligned to interval start" else "fixed baseline feature")
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
      interval=ds$.time[keep],strategy=arm,
      start=(ds$.time[keep]-1)*step/24,stop=ds$.time[keep]*step/24,
      event=as.integer(D[cell][keep]),weight=W[cell][keep])
  }))
}
ccw_frames <- landmark_frames <- list()
for(site in names(site_data)) {
  message("Reconstructing preserved weights and risk sets: ",site)
  ds <- site_data[[site]]
  x <- capture(g6_local_components(ds,cfg,common),"CCW",site)
  ccw_frames[[site]] <- make_frame(ds,x,x$uid)
  if(nrow(land[[site]])) {
    ls <- capture(g6_landmark_local_stats(ds,cfg,p_num),"Landmark",site)
    landmark_frames[[site]] <- make_frame(ds,ls$local,as.character(land[[site]]$.uid))
  }
}
frames <- list(CCW=rbindlist(ccw_frames),Landmark=rbindlist(landmark_frames))
rm(ccw_frames,landmark_frames)
fits <- results <- counts <- model_terms <- list()
for(method in names(frames)) {
  message("Fitting weighted Cox model with current time-varying covariates: ",method)
  z <- merge(frames[[method]],covariates,by=c(".uid","interval"),all.x=TRUE,sort=FALSE)
  if(nrow(z)!=nrow(frames[[method]]) || anyNA(z[,..candidate])) stop("Cox covariate alignment failed")
  if(any(!is.finite(z$weight)|z$weight<=0) || any(z$stop<=z$start)) stop("Invalid Cox rows")
  # .uid is the original patient/block key, shared by both clones.
  z[,subject:=.uid]
  active <- retained[vapply(retained,function(v) uniqueN(z[[v]])>1L,logical(1))]
  form <- reformulate(c("strategy",active,"strata(site)"),response="Surv(start, stop, event)")
  model_terms[[method]] <- data.table(method=method,term=c("strategy",active,"strata(site)"))
  fit <- capture(coxph(form,data=z,weights=weight,cluster=subject,robust=TRUE,
    ties="breslow",x=TRUE,y=TRUE,model=TRUE),method,"cox_fit")
  if(!is.finite(coef(fit)["strategy"])) stop("Strategy coefficient not estimable for ",method)
  sm <- summary(fit)$coefficients
  beta <- coef(fit)["strategy"]; se <- sqrt(vcov(fit)["strategy","strategy"])
  results[[method]] <- data.table(method=paste0(method," + weighted Cox, current TV adjusted"),
    patients=uniqueN(z$subject),grace_hours=tau*step,hazard_ratio=unname(exp(beta)),
    ci_lower=unname(exp(beta-1.96*se)),ci_upper=unname(exp(beta+1.96*se)),
    log_hr_robust_se=unname(se),p_value=unname(2*pnorm(-abs(beta/se))),
    outcome_covariates=length(active),adjustment="baseline + current time-varying covariates",
    ties="breslow",time_scale="days since original time zero; six-hour interval endpoints")
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
write_agg(rbindlist(model_terms),"cox_actual_model_terms.csv")
write_agg(unique(warnings),"model_warnings.csv")
saveRDS(list(fits=fits,design=design,adjustment="baseline + current time-varying covariates"),
        file.path(sensitive,"weighted_cox_models.rds"))
p <- ggplot(result,aes(hazard_ratio,method))+
  geom_vline(xintercept=1,linetype=2,color="grey50")+
  geom_segment(aes(x=ci_lower,xend=ci_upper,yend=method))+geom_point(size=3)+scale_x_log10()+
  labs(title="Weighted Cox strategy HRs: current time-varying adjustment",x="Hazard ratio (95% robust CI)",y=NULL)+theme_minimal(base_size=12)
ggsave(file.path(agg,"weighted_cox_hazard_ratios.png"),p,width=12,height=4.5,dpi=180)
writeLines(c(
  "Two methods: pooled weighted adjusted Cox following CCW, and pooled weighted adjusted landmark Cox.",
  "Prepared directly from raw CLIF tables and scripts 01/04 outputs; no execution or outputs of 07 are required.",
  "Control CCW risk sets allow later initiation; treatment/adherence weights freeze at grace end. Other main-pipeline rules are preserved.",
  "Events are placed at their six-hour interval endpoints; this is not an exact-event-time Cox analysis.",
  "CCW strategy is the assigned clone arm, not observed time-varying treatment. Clone rows share original-subject clustering.",
  "Landmark begins follow-up after the grace period among survivors; treatment group is fixed by initiation through the landmark.",
  "Both control definitions allow initiation after grace; landmark conditions on survival to the landmark whereas CCW uses the inherited clone risk sets from time zero.",
  "Both Cox models adjust for baseline features PLUS the CURRENT interval's clinical covariates and missingness indicators; hospital baseline hazards are stratified.",
  "Cox covariates are joined by patient AND interval and aligned to interval start. Existing lookback, plausibility, imputation, and scaling rules are retained.",
  "Current and lagged histories remain in the weights. Lagged measurements are not additionally included in the Cox outcome formula.",
  "The reported HR is conditional on baseline features, current measured history, and site. Treatment-affected covariates can change interpretation and introduce bias; this is not automatically a total or direct causal effect.",
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
