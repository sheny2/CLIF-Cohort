#!/usr/bin/env Rscript
# Complete-case sensitivity after 01--05; no execution of 07/08/09 is required.
# Usage: Rscript 0701_complete_case_five_methods.R CLIF_DATA_DIR [COHORT_DIR] [TREATMENT_DIR] [NEW_OUTPUT_DIR] [DATA_CUTOFF_DATE]
# Baseline missingness: patient level. Time-varying missingness: all observed
# intervals through day 28, before deletion. Keep <30%; omit >=30%.
# Exclude entire blocks if ANY retained model variable is missing in ANY row.
# No covariate imputation or missingness indicators. Six-hour lookbacks remain.
# Requires 06_generalized_five_methods.R and complete_case_sensitivity_helpers.R
# alongside this script. Outputs called preselection_* describe the original
# candidates BEFORE complete-case restriction, not the final sensitivity cohort.
suppressPackageStartupMessages({library(arrow);library(data.table);library(dplyr);library(ggplot2);library(survival)})
a <- commandArgs(TRUE)
sf <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(normalizePath(gsub("~+~"," ",sf,fixed=TRUE)))
source(file.path(root,"06_generalized_five_methods.R"))
source(file.path(root,"complete_case_sensitivity_helpers.R"))
arg <- function(i,default) if(length(a)>=i && nzchar(a[i])) a[i] else default
data_dir <- arg(1,Sys.getenv("CLIF_DATA_DIR","/projects/emulate_icu/rclif_2_1"))
cohort_dir <- arg(2,Sys.getenv("COHORT_OUTPUT_DIR",file.path(root,"pipeline_output/sensitive/cohort")))
treatment_dir <- arg(3,Sys.getenv("STEROID_TREATMENT_DIR",file.path(root,"pipeline_output/sensitive/steroid_treatment")))
out <- arg(4,file.path(root,"pipeline_output",paste0("0701_complete_case_",format(Sys.time(),"%Y%m%d_%H%M%S"))))
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
state <- prepare_complete_case_inputs(data_dir,cohort_dir,treatment_dir,sensitive,agg,grace_hours,cutoff_text,0L)
invisible(list2env(state,envir=environment()));rm(state)
writeLines(c("Strict complete-case sensitivity, not the primary analysis.",
 "Retain confounders with missingness <30%; omit >=30%, assessed once before case deletion.",
 "Baseline: one row per patient; time-varying: full observed follow-up, before imputation (none performed).",
 "Exclude entire blocks with any missing retained confounder at any observed interval. No individual intervals are discarded.",
 "The SOFA proxy requires all five component scores; missing components are not scored zero.",
 "Unknown sex/location remain missing until the restriction. Absence of a qualifying comorbidity code remains zero as in the original phenotype.",
 "No missingness indicators enter any model. Initial lag=0 and existing bounded lookbacks are retained.",
 "Selection uses future measurement availability and can select on treatment, discharge and outcome; interpret as a selected-cohort sensitivity analysis.",
 "Inspect sensitivity_variable_selection.csv and sensitivity_attrition.csv before estimates.",
 "Files prefixed preselection_ summarize candidates before restriction. Sensitive outputs contain the final restricted prepared data."
),file.path(out,"SENSITIVITY_RULES.txt"))
site_summary <- cohort[,.(patients=.N,deaths_28d=sum(delta)),by=site]
cohort_mortality_summary <- data.table(measure=c("Sensitivity patients","28-day deaths"),n=c(nrow(cohort),sum(cohort$delta)))
write_shareable_csv(site_summary,"sensitivity_cohort_by_site.csv")
write_shareable_csv(cohort_mortality_summary,"cohort_mortality_summary.csv")
message("Running the five methods...")
analysis_warnings <- character()
analysis <- withCallingHandlers(
  run_pneumonia_five_methods(
    data = as.data.frame(long),
    id_col = "hospital_block_id",
    site_col = "site",
    time_col = "m",
    event_col = "event_m",
    censor_col = "censor_m",
    initiation_col = "init_m",
    baseline_covariates = baseline_covariates,
    time_varying_covariates = time_varying_covariates,
    grace_period = tau,
    horizon = M,
    weight_truncation = c(0.01, 0.99),
    include_lagged_tv = TRUE,
    lag_initial = 0,
    model_terms = list(
      initiation_denominator = c(
        ".time", baseline_covariates, time_varying_covariates,
        sensitivity_lags(time_varying_covariates)
      ),
      initiation_numerator = ".time",
      censoring_denominator = c(".time", baseline_covariates, ".trt"),
      censoring_numerator = ".time"
    ),
    verbose = TRUE
  ),
  warning = function(w) {
    analysis_warnings <<- c(analysis_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

# Federated diagnostics expose only aggregate weight distributions and the
# exact weighted risk/event-set totals transmitted to the coordinator.
weight_diagnostics <- rbindlist(lapply(
  names(analysis$site_components),
  function(site_name) {
    component <- analysis$site_components[[site_name]]
    rbindlist(lapply(c("treatment", "no_treatment"), function(strategy) {
      values <- if (strategy == "treatment") {
        component$W1[component$R1]
      } else component$W0[component$R0]
      values <- values[is.finite(values)]
      data.table(
        site = site_name,
        strategy = strategy,
        weighted_person_intervals = length(values),
        minimum = if (length(values)) min(values) else NA_real_,
        q01 = if (length(values)) quantile(values, 0.01, names = FALSE) else NA_real_,
        median = if (length(values)) median(values) else NA_real_,
        q99 = if (length(values)) quantile(values, 0.99, names = FALSE) else NA_real_,
        maximum = if (length(values)) max(values) else NA_real_,
        effective_sample_size = if (length(values) && sum(values^2) > 0) {
          sum(values)^2 / sum(values^2)
        } else NA_real_
      )
    }))
  }
))
write_shareable_csv(
  weight_diagnostics, "federated_ccw_weight_diagnostics.csv"
)

# ---- Covariate-balance diagnostics -----------------------------------------
#
# Weight ranges and effective sample sizes do not establish exchangeability.
# The diagnostics below calculate interval-specific standardized mean
# differences (SMDs) using the exact risk sets and weights used by each
# estimator.  The plotted value is the largest absolute SMD during the four
# treatment-decision intervals in the configured grace period.  The complete
# interval-specific results are also saved so balance later in follow-up can be
# inspected without exposing patient-level data.
#
# Federated and pooled site-stratified CCW use the same local nuisance models,
# weights, and weighted risk sets by design, so a single CCW Love plot applies
# to both estimators.  Per-protocol is unweighted; its before/after points are
# intentionally identical.

balance_variables <- unique(c(
  analysis$config$baseline,
  analysis$config$tv,
  analysis$config$lag_names
))

weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(keep) || sum(w[keep]) <= 0) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

unweighted_variance_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
}

balance_one_variable <- function(x1, x0, w1, w0) {
  keep1 <- is.finite(x1) & is.finite(w1) & w1 >= 0
  keep0 <- is.finite(x0) & is.finite(w0) & w0 >= 0
  x1 <- x1[keep1]
  x0 <- x0[keep0]
  w1 <- w1[keep1]
  w0 <- w0[keep0]
  if (length(x1) < 2L || length(x0) < 2L ||
      sum(w1) <= 0 || sum(w0) <= 0) {
    return(c(
      smd_before = NA_real_, smd_after = NA_real_,
      n1 = length(x1), n0 = length(x0),
      ess1 = NA_real_, ess0 = NA_real_
    ))
  }

  # Use the same unweighted pooled standardization factor before and after
  # weighting.  This ensures changes in the SMD reflect changes in means from
  # weighting rather than changes in the denominator.
  standardizer <- sqrt(
    (unweighted_variance_safe(x1) + unweighted_variance_safe(x0)) / 2
  )
  if (!is.finite(standardizer) || standardizer <= 0) {
    return(c(
      smd_before = NA_real_, smd_after = NA_real_,
      n1 = length(x1), n0 = length(x0),
      ess1 = sum(w1)^2 / sum(w1^2),
      ess0 = sum(w0)^2 / sum(w0^2)
    ))
  }

  c(
    smd_before = (mean(x1) - mean(x0)) / standardizer,
    smd_after = (
      weighted_mean_safe(x1, w1) - weighted_mean_safe(x0, w0)
    ) / standardizer,
    n1 = length(x1),
    n0 = length(x0),
    ess1 = sum(w1)^2 / sum(w1^2),
    ess0 = sum(w0)^2 / sum(w0^2)
  )
}

component_arm_data <- function(site_name, interval, method_key, arm) {
  component <- analysis$site_components[[site_name]]
  site_interval <- analysis$data[
    analysis$data$.site == site_name & analysis$data$.time == interval,
    ,
    drop = FALSE
  ]
  row_index <- match(component$uid, as.character(site_interval$.uid))
  observed <- !is.na(row_index)
  early_treatment <- component$S <= analysis$config$tau

  if (method_key == "ccw") {
    eligible <- if (arm == 1L) {
      component$R1[, interval]
    } else component$R0[, interval]
    weights <- if (arm == 1L) {
      component$W1[, interval]
    } else component$W0[, interval]
  } else if (method_key == "ipcw_no_clone") {
    eligible <- observed & if (arm == 1L) early_treatment else !early_treatment
    weights <- if (arm == 1L) {
      component$W1[, interval]
    } else component$W0[, interval]
  } else if (method_key == "perprotocol_naive") {
    eligible <- observed & if (arm == 1L) early_treatment else !early_treatment
    weights <- rep(1, component$n)
  } else {
    stop("Unknown component-based balance method: ", method_key)
  }

  keep <- which(eligible & observed & is.finite(weights))
  if (!length(keep)) return(NULL)
  answer <- as.data.table(
    as.data.frame(site_interval)[
      row_index[keep], balance_variables, drop = FALSE
    ]
  )
  answer[, `.balance_weight` := weights[keep]]
  answer
}

component_method_balance <- function(method_key, intervals = seq_len(M)) {
  rbindlist(lapply(intervals, function(interval) {
    arm1 <- rbindlist(lapply(
      names(analysis$site_components),
      component_arm_data,
      interval = interval,
      method_key = method_key,
      arm = 1L
    ), fill = TRUE)
    arm0 <- rbindlist(lapply(
      names(analysis$site_components),
      component_arm_data,
      interval = interval,
      method_key = method_key,
      arm = 0L
    ), fill = TRUE)
    if (!nrow(arm1) || !nrow(arm0)) return(NULL)

    rbindlist(lapply(balance_variables, function(variable) {
      values <- balance_one_variable(
        arm1[[variable]], arm0[[variable]],
        arm1$.balance_weight, arm0$.balance_weight
      )
      data.table(
        method_key = method_key,
        interval = interval,
        hours_from_time_zero = interval * step_hours,
        variable = variable,
        smd_before = unname(values[["smd_before"]]),
        smd_after = unname(values[["smd_after"]]),
        n_treatment = as.integer(values[["n1"]]),
        n_no_treatment = as.integer(values[["n0"]]),
        effective_n_treatment = unname(values[["ess1"]]),
        effective_n_no_treatment = unname(values[["ess0"]])
      )
    }))
  }), fill = TRUE)
}

landmark_balance <- function() {
  d <- analysis$data
  tau_balance <- analysis$config$tau
  landmark_by_site <- lapply(split(d, d$.site), function(site_data) {
    site_data[
      site_data$.time == tau_balance &
        site_data$.event == 0L & site_data$.censor == 0L,
      ,
      drop = FALSE
    ]
  })
  treated_total <- sum(vapply(
    landmark_by_site,
    function(x) sum(x$.S <= tau_balance),
    numeric(1)
  ))
  eligible_total <- sum(vapply(landmark_by_site, nrow, integer(1)))
  p_num <- g6_clamp_prob(treated_total / eligible_total)
  propensity_terms <- setdiff(
    analysis$config$terms$initiation_denominator,
    c(".time", ".time_f", ".trt", ".trt_lag")
  )

  landmark <- rbindlist(Map(function(site_data, site_name) {
    if (!nrow(site_data)) return(NULL)
    site_data$.A <- as.integer(site_data$.S <= tau_balance)
    if (length(unique(site_data$.A)) < 2L) {
      p_den <- rep(g6_clamp_prob(mean(site_data$.A)), nrow(site_data))
    } else {
      fit <- g6_fit_binary(
        site_data,
        ".A",
        propensity_terms,
        analysis$config$learners$initiation
      )
      p_den <- g6_predict_binary(fit, site_data)
    }
    weights <- ifelse(
      site_data$.A == 1L,
      p_num / p_den,
      (1 - p_num) / (1 - p_den)
    )
    limits <- quantile(
      weights[is.finite(weights)],
      probs = analysis$config$trunc,
      names = FALSE,
      na.rm = TRUE
    )
    site_data$.balance_weight <- pmin(pmax(weights, limits[1L]), limits[2L])
    site_data$.balance_site <- site_name
    as.data.table(site_data)
  }, landmark_by_site, names(landmark_by_site)), fill = TRUE)

  arm1 <- landmark[.A == 1L]
  arm0 <- landmark[.A == 0L]
  rbindlist(lapply(balance_variables, function(variable) {
    values <- balance_one_variable(
      arm1[[variable]], arm0[[variable]],
      arm1$.balance_weight, arm0$.balance_weight
    )
    data.table(
      method_key = "landmark_ipw",
      interval = tau_balance,
      hours_from_time_zero = tau_balance * step_hours,
      variable = variable,
      smd_before = unname(values[["smd_before"]]),
      smd_after = unname(values[["smd_after"]]),
      n_treatment = as.integer(values[["n1"]]),
      n_no_treatment = as.integer(values[["n0"]]),
      effective_n_treatment = unname(values[["ess1"]]),
      effective_n_no_treatment = unname(values[["ess0"]])
    )
  }))
}

message("Calculating covariate-balance diagnostics...")
balance_by_interval <- rbindlist(list(
  component_method_balance("ccw"),
  landmark_balance(),
  component_method_balance("ipcw_no_clone"),
  component_method_balance("perprotocol_naive")
), fill = TRUE)

balance_method_labels <- c(
  ccw = "CCW (federated and pooled site-stratified)",
  landmark_ipw = "Landmark IPW",
  ipcw_no_clone = "IPCW without cloning",
  perprotocol_naive = "Per-protocol, unweighted"
)
balance_by_interval[, method := unname(balance_method_labels[method_key])]
balance_by_interval[, within_grace_period := interval <= tau]
write_shareable_csv(
  balance_by_interval,
  "covariate_balance_by_method_and_interval.csv"
)

# The Love plots use the worst absolute imbalance observed in the grace period.
# This avoids concealing a poorly balanced treatment-decision interval by
# averaging SMDs over time.  Landmark has one assessment at the end of the grace period.
balance_love_data <- balance_by_interval[
  within_grace_period == TRUE &
    is.finite(smd_before) & is.finite(smd_after),
  .(
    abs_smd_before = max(abs(smd_before), na.rm = TRUE),
    abs_smd_after = max(abs(smd_after), na.rm = TRUE),
    interval_of_worst_before = interval[which.max(abs(smd_before))],
    interval_of_worst_after = interval[which.max(abs(smd_after))]
  ),
  by = .(method_key, method, variable)
]
write_shareable_csv(balance_love_data, "covariate_balance_love_plot_data.csv")

balance_summary <- balance_love_data[, .(
  variables_assessed = .N,
  median_abs_smd_before = median(abs_smd_before),
  median_abs_smd_after = median(abs_smd_after),
  maximum_abs_smd_before = max(abs_smd_before),
  maximum_abs_smd_after = max(abs_smd_after),
  variables_over_0_10_before = sum(abs_smd_before > 0.10),
  variables_over_0_10_after = sum(abs_smd_after > 0.10),
  variables_over_0_05_before = sum(abs_smd_before > 0.05),
  variables_over_0_05_after = sum(abs_smd_after > 0.05)
), by = .(method_key, method)]
write_shareable_csv(balance_summary, "covariate_balance_summary.csv")

pretty_balance_name <- function(x) {
  x <- sub("^lag_", "Prior interval: ", x)
  x <- gsub("_model$", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(x)
}

for (method_key_value in names(balance_method_labels)) {
  plot_data <- copy(balance_love_data[method_key == method_key_value])
  if (!nrow(plot_data)) next
  plot_data[, variable_label := pretty_balance_name(variable)]
  plot_data[, variable_label := factor(
    variable_label,
    levels = variable_label[order(pmax(abs_smd_before, abs_smd_after))]
  )]
  long_plot <- melt(
    plot_data,
    id.vars = c("variable", "variable_label"),
    measure.vars = c("abs_smd_before", "abs_smd_after"),
    variable.name = "weighting",
    value.name = "absolute_smd"
  )
  long_plot[, weighting := factor(
    weighting,
    levels = c("abs_smd_before", "abs_smd_after"),
    labels = c("Before weighting", "After weighting")
  )]

  love_plot <- ggplot(
    long_plot,
    aes(x = absolute_smd, y = variable_label, color = weighting)
  ) +
    geom_vline(xintercept = 0.10, linetype = "dashed", color = "#B22222") +
    geom_vline(xintercept = 0.05, linetype = "dotted", color = "grey50") +
    geom_line(
      aes(group = variable_label),
      color = "grey75",
      linewidth = 0.45
    ) +
    geom_point(size = 2.4) +
    scale_color_manual(values = c(
      "Before weighting" = "#D55E00",
      "After weighting" = "#0072B2"
    )) +
    labs(
      title = paste0(balance_method_labels[[method_key_value]],
                     ": covariate balance"),
      subtitle = if (method_key_value == "landmark_ipw") {
        "Absolute SMD at the end-of-grace-period landmark"
      } else if (method_key_value == "perprotocol_naive") {
        paste(
          "Worst absolute SMD across treatment-decision intervals in the grace period;",
          "before and after are identical because this method is unweighted"
        )
      } else {
        paste(
          "Worst absolute SMD across treatment-decision intervals in the grace period;",
          "dashed line = 0.10"
        )
      },
      x = "Absolute standardized mean difference",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  ggsave(
    file.path(
      shareable_dir,
      paste0("love_plot_", method_key_value, ".png")
    ),
    love_plot,
    width = 9,
    height = max(7, 0.27 * nrow(plot_data) + 2.5),
    dpi = 300,
    limitsize = FALSE
  )
}

local_stats <- lapply(analysis$site_components, g6_local_stats)
sum_stat <- function(name) Reduce(`+`, lapply(local_stats, `[[`, name))
risk_event_counts <- data.table(
  interval = seq_len(M),
  day_end = seq_len(M) * step_hours / 24,
  weighted_risk_treatment = sum_stat("rw1"),
  weighted_events_treatment = sum_stat("dw1"),
  weighted_risk_no_treatment = sum_stat("rw0"),
  weighted_events_no_treatment = sum_stat("dw0")
)
write_sensitive_csv(
  risk_event_counts, "federated_ccw_weighted_risk_event_counts.csv"
)

method_labels <- c(
  fed_ccw = "Federated CCW",
  pooled_ccw_site_stratified = "Pooled CCW: site-stratified",
  landmark_ipw = "Federated landmark IPW",
  ipcw_no_clone = "IPCW without cloning",
  perprotocol_naive = "Per-protocol, unweighted"
)

summary <- rbindlist(lapply(names(analysis$fits), function(method) {
  fit <- analysis$fits[[method]]
  data.table(
    method_key = method,
    method = unname(method_labels[method]),
    death_risk_treatment = as.numeric(fit$psi1),
    death_risk_no_treatment = as.numeric(fit$psi0),
    risk_difference = as.numeric(fit$RD),
    risk_difference_se = fit$SE_RD %||% NA_real_,
    risk_difference_lower = if (!is.null(fit$CI_RD)) fit$CI_RD[1L] else NA_real_,
    risk_difference_upper = if (!is.null(fit$CI_RD)) fit$CI_RD[2L] else NA_real_,
    risk_ratio = as.numeric(fit$RR),
    log_risk_ratio_se = fit$SE_logRR %||% NA_real_,
    risk_ratio_lower = if (!is.null(fit$CI_RR)) fit$CI_RR[1L] else NA_real_,
    risk_ratio_upper = if (!is.null(fit$CI_RR)) fit$CI_RR[2L] else NA_real_,
    odds_ratio = as.numeric(fit$OR),
    log_odds_ratio_se = fit$SE_logOR %||% NA_real_,
    odds_ratio_lower = if (!is.null(fit$CI_OR)) fit$CI_OR[1L] else NA_real_,
    odds_ratio_upper = if (!is.null(fit$CI_OR)) fit$CI_OR[2L] else NA_real_,
    rmst_treatment_days = as.numeric(fit$RMST1) * step_hours / 24,
    rmst_treatment_se_days = (fit$SE_RMST1 %||% NA_real_) * step_hours / 24,
    rmst_no_treatment_days = as.numeric(fit$RMST0) * step_hours / 24,
    rmst_no_treatment_se_days = (fit$SE_RMST0 %||% NA_real_) * step_hours / 24,
    rmst_difference_days = as.numeric(fit$RMST_diff) * step_hours / 24,
    rmst_difference_se_days = (fit$SE_RMSTdiff %||% NA_real_) * step_hours / 24,
    rmst_difference_lower_days = if (!is.null(fit$CI_RMSTdiff)) {
      fit$CI_RMSTdiff[1L] * step_hours / 24
    } else NA_real_,
    rmst_difference_upper_days = if (!is.null(fit$CI_RMSTdiff)) {
      fit$CI_RMSTdiff[2L] * step_hours / 24
    } else NA_real_,
    landmark_eligible = fit$landmark_eligible %||% NA_integer_,
    variance_type = fit$variance_type %||% NA_character_
  )
}))
write_shareable_csv(summary, "five_method_results.csv")
write_shareable_csv(as.data.table(analysis$estimates), "five_method_tidy_estimates.csv")
write_shareable_csv(
  as.data.table(analysis$identity_check),
  "fed_vs_site_stratified_identity_check.csv"
)

curves <- rbindlist(lapply(names(analysis$fits), function(method) {
  fit <- analysis$fits[[method]]
  data.table(
    method_key = method,
    method = unname(method_labels[method]),
    strategy = rep(
      c("Glucocorticoid strategy", "No-glucocorticoid strategy"),
      each = M + 1L
    ),
    day = rep(seq.int(0, M) * step_hours / 24, times = 2L),
    survival = c(fit$S1, fit$S0)
  )
}))
curves[, method := factor(method, levels = unname(method_labels))]
write_shareable_csv(curves, "five_method_survival_curves.csv")

survival_plot <- ggplot(
  curves,
  aes(x = day, y = survival, color = strategy, linetype = strategy)
) +
  geom_step(linewidth = 0.75) +
  facet_wrap(~ method, ncol = 2) +
  scale_color_manual(values = c(
    "Glucocorticoid strategy" = "#C44E52",
    "No-glucocorticoid strategy" = "#4C72B0"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Severe CAP: five-method survival comparison",
    subtitle = paste(
      "Qualifying glucocorticoid initiation within the grace period vs",
      "no initiation through day 28; five Hopkins hospitals"
    ),
        x = "Days since respiratory-failure start",
    y = "Estimated survival probability",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(
  file.path(shareable_dir, "five_method_survival_comparison.png"),
  survival_plot,
  width = 11,
  height = 10,
  dpi = 250
)

effect_data <- rbindlist(list(
  summary[, .(
    method,
    estimand = "28-day risk difference",
    estimate = risk_difference,
    lower = risk_difference_lower,
    upper = risk_difference_upper
  )],
  summary[, .(
    method,
    estimand = "28-day RMST difference (days)",
    estimate = rmst_difference_days,
    lower = rmst_difference_lower_days,
    upper = rmst_difference_upper_days
  )]
))
effect_data[, method := factor(method, levels = rev(unname(method_labels)))]
effect_plot <- ggplot(effect_data, aes(x = estimate, y = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0.18,
    orientation = "y",
    na.rm = TRUE
  ) +
  geom_point(size = 2.5, color = "#7A3E9D") +
  facet_wrap(~ estimand, scales = "free_x", ncol = 1) +
  labs(
    title = "Estimated glucocorticoid strategy contrasts",
    subtitle = "Positive risk difference indicates higher mortality; positive RMST difference indicates longer survival",
    x = "Estimate",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(
  file.path(shareable_dir, "five_method_effect_comparison.png"),
  effect_plot,
  width = 9,
  height = 7,
  dpi = 250
)

writeLines(unique(analysis_warnings), file.path(shareable_dir, "method_warnings.txt"))

compact_fits <- lapply(analysis$fits, function(fit) {
  fit[intersect(names(fit), c(
    "S1", "S0", "psi1", "psi0", "RD", "RR", "OR",
    "RMST1", "RMST0", "RMST_diff", "SE_RD", "SE_logRR",
    "SE_logOR", "SE_RMSTdiff", "SE_RMST1", "SE_RMST0", "CI_RD",
    "CI_RR", "CI_OR", "CI_RMSTdiff", "variance_type",
    "landmark_eligible", "landmark_treated"
  ))]
})
saveRDS(
  list(summary = summary, curves = curves, fits = compact_fits, design = design),
  file.path(sensitive_dir, "five_method_results.rds")
)

cat("\nObserved cohort mortality summary:\n")
print(cohort_mortality_summary)
cat("\nHospital summary:\n")
print(site_summary)
cat("\nFive-method summary:\n")
print(summary)
cat("\nSensitive results saved to:\n", normalizePath(sensitive_dir), "\n")
cat("Shareable results saved to:\n", normalizePath(shareable_dir), "\n")

writeLines("Complete",file.path(out,"COMPLETE.txt"))
