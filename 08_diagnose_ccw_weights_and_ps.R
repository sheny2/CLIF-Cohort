#!/usr/bin/env Rscript
# Post-07 diagnostics; existing scripts and analysis results are not modified.
# Usage (on the server, from CLIF-Cohort):
#   Rscript 08_diagnose_ccw_weights_and_ps.R \
#     /path/to/sensitive/ccw_analysis /path/to/new/ccw_diagnostics
# Inputs: pneumonia_ccw_long.parquet and analysis_design.rds saved by 07.
# No raw CLIF tables are needed. The saved design determines grace period.
# Output: aggregate/ (unsuppressed summaries/plots) and sensitive/ (row-level
# probabilities, weights, and extreme-weight trajectories; keep on the server).
# Optional environment variables:
#   RUN_PS_MODEL_COMPARISON=1 (default); PS_CV_FOLDS=3 (2--5); DIAGNOSTIC_SEED=32284
# Candidate models use the EXISTING covariates: original logistic, spline
# logistic, and ridge logistic if glmnet is installed. Cross-validation holds
# out entire patients within sites. It assesses initiation prediction, not
# causal identification, and does not automatically select a primary model.

suppressPackageStartupMessages({library(arrow); library(data.table); library(ggplot2)})
argv <- commandArgs(TRUE)
sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
script_dir <- dirname(normalizePath(gsub("~+~", " ", sf, fixed = TRUE)))
source(file.path(script_dir, "06_generalized_five_methods.R"))
input_dir <- if (length(argv)) argv[1L] else
  file.path(script_dir, "pipeline_output", "sensitive", "ccw_analysis")
output_dir <- if (length(argv) >= 2L) argv[2L] else
  file.path(dirname(input_dir), paste0("ccw_diagnostics_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (normalizePath(input_dir, mustWork = FALSE) == normalizePath(output_dir, mustWork = FALSE))
  stop("Use a new diagnostic output directory, separate from the script-07 input directory.")
aggregate_dir <- file.path(output_dir, "aggregate")
sensitive_dir <- file.path(output_dir, "sensitive")
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitive_dir, recursive = TRUE, showWarnings = FALSE)
write_agg <- function(x, name) {
  if (any(c("patient_id", "hospital_block_id", ".uid", ".id") %in% names(x)))
    stop("Identifier found in aggregate output: ", name)
  fwrite(x, file.path(aggregate_dir, name), na = "")
}
for (f in c("analysis_design.rds", "pneumonia_ccw_long.parquet"))
  if (!file.exists(file.path(input_dir, f))) stop("Missing script-07 output: ", f)
design <- readRDS(file.path(input_dir, "analysis_design.rds"))
long <- as.data.frame(read_parquet(file.path(input_dir, "pneumonia_ccw_long.parquet")))
M <- as.integer(design$horizon_intervals); tau <- as.integer(design$tau_intervals)
step <- as.numeric(design$step_hours)
if (!length(M) || !length(tau) || tau < 1 || tau >= M) stop("Invalid saved design.")
seed <- as.integer(Sys.getenv("DIAGNOSTIC_SEED", "32284"))
folds <- as.integer(Sys.getenv("PS_CV_FOLDS", "3"))
if (!is.finite(seed) || !is.finite(folds) || !folds %in% 2:5)
  stop("DIAGNOSTIC_SEED must be an integer and PS_CV_FOLDS must be 2--5.")
run_comparison <- Sys.getenv("RUN_PS_MODEL_COMPARISON", "1") == "1"
prepared <- g6_prepare_data(long, "site", "hospital_block_id", "m", "event_m",
                           "censor_m", "init_m", NULL, design$baseline_covariates,
                           design$time_varying_covariates, M, lag_initial = 0)
d <- prepared$data
rm(long); invisible(gc())
terms <- list(
  initiation_denominator = c(".time", design$baseline_covariates,
    design$time_varying_covariates, prepared$lag_names),
  initiation_numerator = ".time",
  censoring_denominator = c(".time", design$baseline_covariates, ".trt"),
  censoring_numerator = ".time"
)
learner <- g6_default_binary_learner()
cfg <- list(tau = tau, horizon = M, terms = terms, trunc = c(0, 1),
            learners = list(initiation = learner, censoring = learner))
common_num <- g6_common_init_numerator(d, M)

# Extend a LOCAL COPY's return value to expose already-computed fitted hazards
# and fit objects. The source file is untouched; assert its expected interface.
# Setting truncation to [0,1] retains the full range of the raw weights.
diagnostic_components <- g6_local_components
b <- as.list(body(diagnostic_components)); last <- b[[length(b)]]
if (!is.call(last) || !identical(last[[1L]], as.name("list")))
  stop("The script-06 component interface changed; review the diagnostic adapter.")
b[[length(b)]] <- as.call(c(as.list(last), list(
  Hden = quote(Hden), Hnum = quote(Hnum), initiation_fit = quote(afit),
  censor_fit = quote(cfit), censor_numerator_fit = quote(cnfit))))
body(diagnostic_components) <- as.call(b)

warnings_log <- data.table(site = character(), stage = character(), message = character())
capture <- function(expr, site, stage) withCallingHandlers(expr, warning = function(w) {
  warnings_log <<- rbind(warnings_log, data.table(site = site, stage = stage,
                                                message = conditionMessage(w)))
  invokeRestart("muffleWarning")
})
q <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) as.numeric(quantile(x, p, names = FALSE)) else rep(NA_real_, length(p))
}
weight_summary <- function(w) {
  finite <- is.finite(w) & w > 0
  z <- w[finite]; n <- length(z)
  qs <- q(z, c(0.01, 0.5, 0.95, 0.99, 0.999))
  ess <- if (n) sum(z)^2 / sum(z^2) else NA_real_
  list(n = length(w), n_invalid = sum(!finite), minimum = if(n) min(z) else NA_real_,
       q01 = qs[1], median = qs[2], q95 = qs[3], q99 = qs[4], q999 = qs[5],
       maximum = if(n) max(z) else NA_real_, mean = if(n) mean(z) else NA_real_,
       cv = if(n > 1) sd(z) / mean(z) else NA_real_, ess = ess,
       ess_fraction = ess / n, n_gt_5 = sum(z > 5), n_gt_10 = sum(z > 10),
       n_gt_20 = sum(z > 20), n_gt_50 = sum(z > 50),
       top_1pct_weight_share = if(n) sum(sort(z, decreasing = TRUE)[seq_len(max(1, ceiling(n*.01)))])/sum(z) else NA_real_)
}
fit_summary <- function(fit, dataset, response, site, component) {
  model <- fit$model
  cf <- if (inherits(model, "glm")) coef(model) else numeric()
  data.table(site = site, component = component, rows = nrow(dataset),
    events = sum(dataset[[response]]), patients = uniqueN(dataset$.uid),
    constant_response = isTRUE(fit$constant),
    coefficients = length(cf), rank = if(length(cf)) model$rank else 0L,
    aliased_coefficients = sum(is.na(cf)),
    converged = if(length(cf)) isTRUE(model$converged) else TRUE,
    largest_abs_coefficient = if(any(is.finite(cf))) max(abs(cf[is.finite(cf)])) else NA_real_)
}
components <- list(); weight_rows <- list(); probability_rows <- list(); fit_rows <- list()
coefficient_rows <- list()
sites <- sort(unique(d$.site))
for (site in sites) {
  message("Reconstructing original models and untruncated weights: ", site)
  ds <- d[d$.site == site, , drop = FALSE]
  x <- capture(diagnostic_components(ds, cfg, common_num), site, "original_model_refit")
  if (any(!is.finite(x$W1[x$obs])) || any(!is.finite(x$W0[x$obs])))
    stop("Non-finite reconstructed weights at ", site, "; inspect model fit.")
  ad <- ds[ds$.time <= pmin(ds$.S, M), , drop = FALSE]
  cd <- ds[ds$.event == 0L, , drop = FALSE]
  fit_rows[[site]] <- rbindlist(list(
    fit_summary(x$initiation_fit, ad, ".init", site, "initiation"),
    fit_summary(x$censor_fit, cd, ".censor", site, "natural_censoring")))
  if (inherits(x$initiation_fit$model, "glm")) {
    co <- capture(summary(x$initiation_fit$model)$coefficients, site, "coefficient_summary")
    coefficient_rows[[site]] <- data.table(site = site, term = rownames(co),
      estimate = co[,1], standard_error = co[,2], z_value = co[,3], p_value = co[,4])
  }
  acell <- cbind(match(as.character(ad$.uid), x$uid), ad$.time)
  pd <- data.table(site = site, .uid = as.character(ad$.uid),
                   interval = ad$.time, initiated = ad$.init,
                   probability = x$Hden[acell], numerator = x$Hnum[acell])
  pd[, probability_before_clamp := capture(
    as.numeric(x$initiation_fit$learner$predict(x$initiation_fit$model, ad)),
    site, "original_probability_prediction")]
  pd[, observed_action_probability := fifelse(initiated == 1L, probability, 1-probability)]
  pd[, stabilized_increment := fifelse(initiated == 1L,
                                      numerator/probability, (1-numerator)/(1-probability))]
  pd[, within_grace := interval <= tau]
  probability_rows[[site]] <- pd
  ids <- as.character(ds$.uid)
  cell <- cbind(match(ids, x$uid), ds$.time)
  for (arm in c("treatment", "no_treatment")) {
    W <- if(arm == "treatment") x$W1 else x$W0
    R <- if(arm == "treatment") x$R1 else x$R0
    # Exactly match script 06: clipping thresholds use all observed cells,
    # including cells outside the strategy-specific risk set.
    lim <- q(W[x$obs], c(.01, .99))
    w <- W[cell]; used <- R[cell]
    weight_rows[[paste(site, arm)]] <- data.table(
      site = site, .uid = ids, interval = ds$.time, strategy = arm,
      in_risk_set = used, raw_weight = w,
      truncated_weight = pmin(pmax(w, lim[1]), lim[2]),
      lower_limit = lim[1], upper_limit = lim[2])
  }
  x[c("Hden", "Hnum", "initiation_fit", "censor_fit", "censor_numerator_fit")] <- NULL
  components[[site]] <- x
}
weights_all <- rbindlist(weight_rows)
probabilities <- rbindlist(probability_rows)
write_agg(rbindlist(fit_rows), "original_model_fit_diagnostics.csv")
if (length(coefficient_rows)) write_agg(rbindlist(coefficient_rows), "original_initiation_model_coefficients.csv")
weights <- weights_all[in_risk_set == TRUE]
write_parquet(weights_all, file.path(sensitive_dir, "weights_by_patient_interval.parquet"))
write_parquet(probabilities, file.path(sensitive_dir, "initiation_probabilities.parquet"))
weights_long <- melt(weights, id.vars = c("site", ".uid", "interval", "strategy"),
                     measure.vars = c("raw_weight", "truncated_weight"),
                     variable.name = "weight_version", value.name = "weight")
by_site_time <- weights_long[, weight_summary(weight), by = .(site, strategy, weight_version, interval)]
pooled_time <- weights_long[, weight_summary(weight), by = .(strategy, weight_version, interval)]
pooled_time[, site := "All sites"]
ess <- rbindlist(list(by_site_time, pooled_time), use.names = TRUE)
ess[, interval_end_hours := interval * step]
write_agg(ess, "weight_and_patient_ess_by_interval.csv")
overall <- weights_long[, weight_summary(weight), by = .(site, strategy, weight_version)]
overall[, unit := "person-intervals; not patient ESS"]
write_agg(overall, "weight_summary_all_person_intervals.csv")
clipping <- weights[, .(risk_set_person_intervals = .N,
  percent_clipped_low = 100*mean(raw_weight < lower_limit),
  percent_clipped_high = 100*mean(raw_weight > upper_limit),
  lower_limit = lower_limit[1], upper_limit = upper_limit[1]), by = .(site, strategy)]
write_agg(clipping, "truncation_thresholds_and_fraction_clipped.csv")

prob_summary <- function(p, action) list(rows = length(p), minimum = min(p),
  q01 = q(p,.01), median = median(p), q99 = q(p,.99), maximum = max(p),
  percent_p_lt_001 = 100*mean(p < .01), percent_p_gt_099 = 100*mean(p > .99),
  percent_observed_action_probability_lt_001 = 100*mean(action < .01))
write_agg(probabilities[, prob_summary(probability, observed_action_probability),
  by = .(site, interval, initiated)], "initiation_probability_summary.csv")
write_agg(probabilities[, .(patients_at_risk = uniqueN(.uid), initiations = sum(initiated),
  initiation_percent = 100*mean(initiated),
  percent_predictions_clamped = 100*mean(probability_before_clamp < 1e-6 |
                                         probability_before_clamp > 1-1e-6)),
  by = .(site, interval)], "initiation_risk_set_counts.csv")
missing_flags <- grep("_missing$", design$time_varying_covariates, value = TRUE)
if (length(missing_flags)) {
  decision_data <- as.data.table(d[d$.time <= pmin(d$.S,tau), , drop = FALSE])
  miss <- melt(decision_data[, c(".site", ".time", missing_flags), with = FALSE],
               id.vars = c(".site", ".time"), variable.name = "variable", value.name = "missing")
  miss <- miss[, .(risk_set_rows = .N, missing_rows = sum(missing),
                   missing_percent = 100*mean(missing)), by = .(.site, .time, variable)]
  setnames(miss, c(".site", ".time"), c("site", "interval"))
  write_agg(miss, "untreated_decision_risk_set_missingness.csv")
}
write_agg(probabilities[, .(rows = .N, largest_increment = max(stabilized_increment),
  q99_increment = q(stabilized_increment,.99),
  percent_increment_gt_5 = 100*mean(stabilized_increment > 5)),
  by = .(site, within_grace, initiated)], "initiation_weight_increment_summary.csv")

# Identify extreme trajectories and descriptive covariate associations. These
# associations are diagnostic, not causal explanations of treatment effects.
top <- weights[, .(peak_raw_weight = max(raw_weight)), by = .(site, strategy, .uid)]
setorder(top, site, strategy, -peak_raw_weight)
top <- top[, head(.SD, 20L), by = .(site, strategy)]
trajectory <- merge(weights, top, by = c("site", "strategy", ".uid"))
trace_vars <- unique(c(design$baseline_covariates, design$time_varying_covariates, prepared$lag_names))
covdata <- as.data.table(d)[, c(".uid", ".time", trace_vars), with = FALSE]
covdata[, .uid := as.character(.uid)]; setnames(covdata, ".time", "interval")
trajectory <- merge(trajectory, covdata, by = c(".uid", "interval"), all.x = TRUE)
trajectory <- merge(trajectory, probabilities, by = c("site", ".uid", "interval"), all.x = TRUE)
write_parquet(trajectory, file.path(sensitive_dir, "top_weight_patient_trajectories.parquet"))
profiles <- rbindlist(lapply(sites, function(site) {
  selected_site <- site
  rbindlist(lapply(c("treatment", "no_treatment"), function(arm) {
    z <- merge(weights[site == selected_site & strategy == arm], covdata,
               by = c(".uid", "interval"))
    threshold <- q(z$raw_weight, .99)
    high <- z$raw_weight >= threshold
    rbindlist(lapply(trace_vars, function(v) {
      value <- z[[v]]; s <- sd(value)
      data.table(site = site, strategy = arm, variable = v, tail_threshold = threshold,
        tail_rows = sum(high), other_rows = sum(!high),
        mean_in_tail = mean(value[high]), mean_other = mean(value[!high]),
        standardized_difference = if(is.finite(s) && s > 0 && any(!high))
          (mean(value[high])-mean(value[!high]))/s else NA_real_)
    }))
  }))
}))
write_agg(profiles, "extreme_weight_covariate_profiles.csv")

# Reuse the same fitted probabilities to vary only truncation; this is a
# sensitivity analysis, not a search for a preferred treatment-effect sign.
settings <- list(none = c(0,1), p005_p995 = c(.005,.995),
                 original_p01_p99 = c(.01,.99), p025_p975 = c(.025,.975))
sensitivity <- rbindlist(lapply(names(settings), function(label) {
  stats <- lapply(components, function(x) {
    z <- x
    for (nm in c("W1", "W0")) {
      lim <- q(z[[nm]][z$obs], settings[[label]])
      z[[nm]] <- pmin(pmax(z[[nm]],lim[1]),lim[2])
    }
    g6_local_stats(z)
  })
  fit <- g6_central_ccw(stats)
  data.table(truncation = label, death_risk_treatment = fit$psi1,
    death_risk_no_treatment = fit$psi0, risk_difference = fit$RD,
    risk_ratio = fit$RR, rd_ci_lower_fixed_weights = fit$CI_RD[1],
    rd_ci_upper_fixed_weights = fit$CI_RD[2])
}))
write_agg(sensitivity, "ccw_truncation_sensitivity.csv")
original_rds <- file.path(input_dir, "five_method_results.rds")
if (file.exists(original_rds)) {
  saved <- readRDS(original_rds)$fits$fed_ccw
  refit <- sensitivity[truncation == "original_p01_p99"]
  if (!is.null(saved)) {
    check <- data.table(quantity = c("psi1","psi0","RD","RR"),
      original = c(saved$psi1,saved$psi0,saved$RD,saved$RR),
      reconstructed = c(refit$death_risk_treatment,refit$death_risk_no_treatment,
                        refit$risk_difference,refit$risk_ratio))
    check[, `:=`(absolute_difference = abs(original-reconstructed),
                  passed = abs(original-reconstructed) < 1e-8)]
    write_agg(check, "original_07_reproduction_check.csv")
    if (!all(check$passed)) warning("Original 07 estimates did not reproduce; inspect reproduction check before interpreting diagnostics.")
  }
}

# Patient-grouped, site-specific CV. No mortality outcome enters the model
# comparison. Scores use out-of-fold predictions; these are NOT used to refit
# the primary treatment-effect analysis or its weights.
group_folds <- function(a, k) {
  units <- unique(a[, .(.uid, ever_initiated = as.integer(is.finite(.S)))])
  units[, fold := sample(rep(seq_len(k), length.out = .N)), by = ever_initiated]
  units$fold[match(a$.uid, units$.uid)]
}
prediction_metrics <- function(y, p) {
  ok <- is.finite(p); y <- y[ok]; p <- p[ok]
  if (!length(y)) return(list(n = 0L, events = 0L, brier = NA_real_,
    log_loss = NA_real_, auc = NA_real_, calibration_intercept = NA_real_, calibration_slope = NA_real_))
  pc <- g6_clamp_prob(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  auc <- if(n1 && n0) (sum(rank(p)[y == 1])-n1*(n1+1)/2)/(n1*n0) else NA_real_
  cal <- if(n1 && n0 && sd(qlogis(pc)) > 0) tryCatch(
    suppressWarnings(coef(glm(y ~ qlogis(pc), family = binomial()))),
    error = function(e) c(NA_real_,NA_real_)) else c(NA_real_,NA_real_)
  list(n = length(y), events = sum(y), brier = mean((y-p)^2),
    log_loss = -mean(y*log(pc)+(1-y)*log1p(-pc)), auc = auc,
    calibration_intercept = unname(cal[1]), calibration_slope = unname(cal[2]))
}
oof <- list(); comparison_status <- list()
if (run_comparison) {
  set.seed(seed)
  for (site in sites) {
    message("Patient-grouped propensity-model cross-validation: ", site)
    ad <- as.data.table(d[d$.site == site & d$.time <= pmin(d$.S,M), , drop = FALSE])
    ad[, .uid := as.character(.uid)]
    k <- min(folds, uniqueN(ad$.uid))
    if (k < 2) next
    fold_id <- group_folds(ad, k)
    nonlinear_vars <- intersect(c(".time", "age_z", "temperature_model", "wbc_model",
                                  "glucose_model", "sf_ratio_model"), terms$initiation_denominator)
    nonlinear_vars <- nonlinear_vars[vapply(nonlinear_vars, function(v) uniqueN(ad[[v]]) >= 6L, logical(1))]
    spline_terms <- terms$initiation_denominator
    spline_terms[match(nonlinear_vars, spline_terms)] <- paste0("splines::ns(",nonlinear_vars,", df=3)")
    models <- c("original_logistic", "spline_logistic", "ridge_logistic")
    for (model_name in models) {
      if (model_name == "ridge_logistic" && !requireNamespace("glmnet", quietly = TRUE)) {
        comparison_status[[paste(site,model_name)]] <- data.table(site=site,model=model_name,fold=NA_integer_,status="skipped: glmnet not installed")
        next
      }
      for (fold in seq_len(k)) {
        train <- as.data.frame(ad[fold_id != fold]); test <- as.data.frame(ad[fold_id == fold])
        result <- tryCatch(capture({
          if (uniqueN(train$.init) < 2L) stop("Training fold has only one initiation class")
          rhs <- if(model_name == "spline_logistic") spline_terms else terms$initiation_denominator
          form <- g6_make_formula(".init", rhs)
          if (model_name == "ridge_logistic") {
            tr <- model.matrix(form,train)[,-1,drop=FALSE]
            te <- model.matrix(form,test)[,-1,drop=FALSE]
            inner_k <- min(5L,uniqueN(train$.uid))
            if (inner_k < 3L) stop("Too few patients for grouped inner ridge CV")
            inner <- group_folds(as.data.table(train),inner_k)
            fit <- glmnet::cv.glmnet(tr,train$.init,family="binomial",alpha=0,
                                     foldid=inner,type.measure="deviance")
            as.numeric(predict(fit,newx=te,s="lambda.min",type="response"))
          } else {
            fit <- glm(form,train,family=binomial())
            if (!isTRUE(fit$converged)) warning("Candidate logistic model did not converge")
            as.numeric(predict(fit,newdata=test,type="response"))
          }
        }, site, paste(model_name,"fold",fold)), error=function(e)e)
        success <- !inherits(result,"error") && length(result)==nrow(test) && all(is.finite(result))
        comparison_status[[paste(site,model_name,fold)]] <- data.table(
          site=site,model=model_name,fold=fold,status=if(success) "ok" else
            if(inherits(result,"error")) conditionMessage(result) else "non-finite or incomplete predictions")
        # Keep failed rows as NA so incomplete CV cannot silently look superior.
        oof[[paste(site,model_name,fold)]] <- data.table(site=site,model=model_name,
          fold=fold,.uid=test$.uid,interval=test$.time,initiated=test$.init,
          probability=if(success) result else rep(NA_real_,nrow(test)))
      }
    }
    write_agg(data.table(site=site,model="spline_logistic",term=spline_terms),
              paste0("spline_model_terms_",gsub("[^A-Za-z0-9]","_",site),".csv"))
  }
}
if (length(comparison_status)) write_agg(rbindlist(comparison_status), "ps_model_comparison_status.csv")
if (length(oof)) {
  cv <- rbindlist(oof); cv[, within_grace := interval <= tau]
  write_parquet(cv,file.path(sensitive_dir,"out_of_fold_initiation_probabilities.parquet"))
  metrics <- cv[, c(list(expected_rows=.N,missing_predictions=sum(!is.finite(probability))),
                    prediction_metrics(initiated,probability)), by=.(site,model,within_grace)]
  write_agg(metrics,"ps_model_comparison_patient_cv.csv")
  cal <- cv[is.finite(probability)]
  cal[, probability_bin := cut(probability,breaks=seq(0,1,.1),include.lowest=TRUE)]
  cal <- cal[,.(n=.N,mean_predicted=mean(probability),observed_initiation=mean(initiated)),
             by=.(site,model,within_grace,probability_bin)]
  write_agg(cal,"out_of_fold_calibration.csv")
  if (nrow(cal[within_grace == TRUE])) {
    p <- ggplot(cal[within_grace == TRUE],aes(mean_predicted,observed_initiation,color=model))+
    geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+
    geom_point(aes(size=n))+geom_line()+facet_wrap(~site)+
    labs(title="Out-of-fold initiation calibration during the grace period",
         x="Mean predicted interval initiation probability",y="Observed initiation proportion")+
    theme_minimal()
    ggsave(file.path(aggregate_dir,"plot_ps_model_calibration.png"),p,width=11,height=7,dpi=180)
  }
}

theme_set(theme_minimal(base_size=12))
probabilities[, initiation_group := fifelse(initiated == 1L,"Initiated in interval","Did not initiate in interval")]
p <- ggplot(probabilities[within_grace == TRUE],aes(probability,fill=initiation_group))+
  geom_histogram(aes(y=after_stat(density)),binwidth=.025,position="identity",alpha=.45)+
  facet_wrap(~site,scales="free_y")+labs(title="Fitted initiation probabilities during the grace period",
    subtitle="Original model, in-sample; untreated risk-set rows only; patients may contribute multiple intervals",
    x="Probability of starting steroids in this interval",y="Density",fill=NULL)
ggsave(file.path(aggregate_dir,"plot_initiation_probability_overlap.png"),p,width=11,height=7,dpi=180)
weights_long[, log10_weight := log10(weight)]
p <- ggplot(weights_long,aes(log10_weight,color=weight_version))+
  geom_freqpoly(bins=70)+facet_grid(site~strategy,scales="free_y")+
  labs(title="Raw and truncated weights in strategy-specific risk sets",
       x="log10(weight): 0 = 1; 1 = 10; 2 = 100",y="Person-interval count",color=NULL)
ggsave(file.path(aggregate_dir,"plot_raw_and_truncated_weights.png"),p,width=11,height=11,dpi=180)
p <- ggplot(ess,aes(interval_end_hours,ess_fraction,color=weight_version,linetype=strategy))+
  geom_line()+facet_wrap(~site)+geom_vline(xintercept=tau*step,linetype=3)+
  labs(title="Patient risk-set ESS as a fraction of patients at each interval",
       subtitle="ESS = sum(w)^2 / sum(w^2); strategies contain overlapping clones, not independent samples",
       x="Hours since time zero",y="ESS / risk-set patients",color=NULL,linetype=NULL)
ggsave(file.path(aggregate_dir,"plot_ess_over_time.png"),p,width=12,height=8,dpi=180)
p <- ggplot(ess[interval <= tau],aes(interval_end_hours,ess,color=weight_version,linetype=strategy))+
  geom_line()+geom_point()+facet_wrap(~site,scales="free_y")+
  labs(title="Patient risk-set ESS during the grace period",x="Interval end, hours",y="ESS",color=NULL,linetype=NULL)
ggsave(file.path(aggregate_dir,"plot_ess_during_grace.png"),p,width=12,height=8,dpi=180)
write_agg(unique(warnings_log),"model_warnings.csv")
write_agg(data.table(parameter=c("input_dir","grace_hours","horizon_days","seed","cv_folds",
  "probability_epsilon","original_truncation","cv_models_requested"),value=as.character(c(
    normalizePath(input_dir),tau*step,design$horizon_days,seed,folds,1e-6,"0.01;0.99",run_comparison))),"diagnostic_settings.csv")
writeLines(c(
  "Start with original_07_reproduction_check.csv (when original results RDS exists).",
  "Then inspect model_warnings, original_model_fit_diagnostics, raw weight tails and overlap plots.",
  "Use original_initiation_model_coefficients and extreme_weight_covariate_profiles to investigate predictors associated with instability; coefficient magnitude alone is not variable importance.",
  "PS here is P(initiation in this interval | not yet initiated, recorded history), not a baseline ever-treated score.",
  "Original hazards are fitted on ALL eligible intervals through day 28, matching script 07, not only the grace period.",
  "Weight summaries and plots use strategy-specific risk sets; clipping thresholds match script 06's all-observed-cell convention.",
  "ESS by interval is patient risk-set ESS; all-interval ESS is person-interval ESS. Never sum arm ESS or interpret repeated intervals as independent patients.",
  "Extreme-weight covariate profiles pool repeated intervals and describe associations, not causal mechanisms.",
  "Small observed-action probabilities, large weight increments, separated/singular fits and sparse site/interval cells can identify sources of instability.",
  "Compare candidate models using COMPLETE patient-grouped CV, calibration and log loss/Brier; high AUC alone does not establish useful causal weights.",
  "CV predictions do not replace the original weights. Before adopting a model, rerun CCW and assess pre-treatment balance, overlap, ESS and truncation sensitivity.",
  "Spline models add nonlinear time, age, temperature, WBC, glucose and S/F terms where there are enough distinct values; ridge shrinks the original terms.",
  "No additional confounders are invented here. New variables require clinically justified, pre-treatment definitions and extraction before fitting.",
  "Assess covariate availability in the actual untreated decision risk sets; global 28-day missingness can be driven by post-discharge intervals.",
  "Contemporaneous covariates after initiation can be treatment affected; poor cross-arm balance on these is not by itself proof of inadequate confounding control.",
  "Truncation sensitivity varies only truncation of the SAME reconstructed models. CIs treat fitted weights/cutoffs as fixed, as in 07.",
  "All aggregate plots/CSVs are unsuppressed; patient-level files are confined to sensitive/. No source files or original results were changed."
),file.path(output_dir,"START_HERE.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir,"SESSION_INFO.txt"))
writeLines(paste("Completed", format(Sys.time(), tz="UTC", usetz=TRUE)),
           file.path(output_dir,"COMPLETE.txt"))
message("Diagnostics complete: ", normalizePath(output_dir))
