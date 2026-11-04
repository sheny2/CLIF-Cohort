#!/usr/bin/env Rscript
# =============================================================================
# Generalized grace-period estimators used by the Hopkins CLIF analysis.
#
# PNEUMONIA FIVE-METHOD WRAPPER: run_pneumonia_five_methods()
#   1. Federated clone-censor-weight (CCW)
#   2. Pooled CCW with site-stratified nuisance models (identity control)
#   3. Federated landmark IPW
#   4. Federated IPW without cloning
#   5. Federated unweighted per-protocol analysis
#
# The file retains reusable g-computation functions from the parent generalized
# natural-censoring implementation, but the Hopkins wrapper below runs exactly
# the five estimators listed here.
#
# INPUT FORMAT
#   Long person-period data with one row per subject and interval, beginning
#   at interval 1 and continuing without gaps until event, natural censoring,
#   or the administrative horizon. A natural-censoring row is retained: L is
#   observed during that interval and censoring occurs at its end. Event and
#   censoring indicators must therefore never both equal 1 on the same row.
#
# INTERVENTIONS
#   g = 1: initiate treatment by `grace_period`; retain the observed start
#          time if it is within the grace period and otherwise set it to tau.
#   g = 0: never initiate treatment through the analysis horizon.
#
# BASIC USE
#   ans <- run_pneumonia_five_methods(
#     data = long_data,
#     site_col = "site", id_col = "id", time_col = "month",
#     event_col = "event", censor_col = "lost_to_followup",
#     initiation_col = "initiated_now",
#     baseline_covariates = c("age", "sex"),
#     time_varying_covariates = c("severity", "biomarker"),
#     grace_period = 3, horizon = 12
#   )
#   ans$estimates
#
# TREATMENT INPUT
#   Supply either:
#     * initiation_col: 1 only in the interval treatment starts; or
#     * treatment_col: an absorbing 0/1 treatment-status indicator.
#
# CUSTOMIZATION
#   `model_terms` can replace the default RHS term vectors. Its recognized
#   names are initiation_denominator, initiation_numerator,
#   censoring_denominator, censoring_numerator, outcome, and transition.
#   Internally generated names available to these vectors are:
#     .time_f, .trt, .trt_lag, and lag_<time-varying-covariate>.
#
#   `learners` can replace any binary GLM with a custom learner. Each element
#   (initiation, censoring, outcome) must be a list containing:
#     fit(formula, data) and predict(model, newdata)
#   where predict returns P(Y=1). See g6_default_binary_learner().
#
# ASSUMPTIONS
#   This is analysis code, not an automatic identification check. Users remain
#   responsible for consistency, positivity, exchangeability, correct temporal
#   ordering, and defensible nuisance/transition models.
#
# INFERENCE AND DEPLOYMENT
#   All five methods are represented as weighted discrete-time Kaplan-Meier
#   estimators and return influence-function standard errors and confidence
#   intervals. As in ArraySim5-6 Uncertainty, fitted nuisance weights, the
#   estimated common numerator, and adaptive truncation cutoffs are treated as
#   fixed. Bootstrap inference can be layered around the wrapper when needed.
#
#   This file is a single-process reference implementation: "federated" means
#   nuisance models are fit separately after splitting by site and only
#   interval summaries are pooled for point estimation. A physically
#   distributed deployment should expose g6_local_components()/the local
#   influence sums as site-side calls and g6_central_ccw() as the coordinator.
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

g6_clamp_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

g6_default_binary_learner <- function() {
  list(
    fit = function(formula, data) {
      glm(formula, data = data, family = binomial())
    },
    predict = function(model, newdata) {
      predict(model, newdata = newdata, type = "response")
    }
  )
}

g6_validate_learner <- function(x, name) {
  if (!is.list(x) || !is.function(x$fit) || !is.function(x$predict))
    stop("Learner '", name, "' must contain fit() and predict() functions.")
  x
}

g6_make_formula <- function(response, terms) {
  terms <- unique(terms[nzchar(terms)])
  reformulate(if (length(terms)) terms else "1", response = response)
}

g6_first_by_unit <- function(d) {
  d[!duplicated(d$.uid), , drop = FALSE]
}

g6_prepare_data <- function(data, site_col, id_col, time_col, event_col,
                            censor_col, initiation_col, treatment_col,
                            baseline_covariates, time_varying_covariates,
                            horizon, lag_initial) {
  if (!is.data.frame(data)) data <- as.data.frame(data)

  user_cols <- c(id_col, time_col, event_col, censor_col, initiation_col,
                 treatment_col, site_col, baseline_covariates,
                 time_varying_covariates)
  user_cols <- unique(user_cols[!is.na(user_cols) & nzchar(user_cols)])
  missing_cols <- setdiff(user_cols, names(data))
  if (length(missing_cols))
    stop("Missing input columns: ", paste(missing_cols, collapse = ", "))

  covars <- unique(c(baseline_covariates, time_varying_covariates))
  if (length(intersect(baseline_covariates, time_varying_covariates)))
    stop("A covariate cannot be both baseline and time-varying.")
  if (any(make.names(covars) != covars))
    stop("Covariate names must be syntactically valid R names.")

  reserved <- c(".site", ".id", ".uid", ".time", ".time_f", ".event",
                ".censor", ".init", ".trt", ".trt_lag", ".S")
  if (length(intersect(names(data), reserved)))
    stop("Input uses reserved internal columns: ",
         paste(intersect(names(data), reserved), collapse = ", "))

  d <- data
  d$.site <- if (is.null(site_col)) "site_1" else as.character(d[[site_col]])
  d$.id <- as.character(d[[id_col]])
  d$.uid <- interaction(d$.site, d$.id, drop = TRUE, lex.order = TRUE)
  d$.time <- as.integer(d[[time_col]])
  d$.event <- as.integer(d[[event_col]])
  d$.censor <- if (is.null(censor_col)) 0L else as.integer(d[[censor_col]])

  if (anyNA(d[, c(".site", ".id", ".time", ".event", ".censor", covars),
                    drop = FALSE]))
    stop("Missing values are not allowed in required columns or covariates.")
  if (any(!d$.event %in% 0:1) || any(!d$.censor %in% 0:1))
    stop("Event and censoring indicators must be binary.")
  if (any(d$.event + d$.censor > 1))
    stop("Event and natural censoring cannot occur in the same interval.")

  d <- d[order(d$.site, d$.id, d$.time), , drop = FALSE]
  rownames(d) <- NULL
  split_rows <- split(seq_len(nrow(d)), d$.uid)
  bad_sequence <- vapply(split_rows, function(ii) {
    !identical(d$.time[ii], seq_len(length(ii)))
  }, logical(1))
  if (any(bad_sequence))
    stop("Every subject must have consecutive rows starting at interval 1.")
  if (any(d$.time < 1L | d$.time > horizon))
    stop("Observed interval values must lie in 1:horizon.")

  terminal_bad <- vapply(split_rows, function(ii) {
    term <- which(d$.event[ii] == 1L | d$.censor[ii] == 1L)
    length(term) > 1L || (length(term) == 1L && term != length(ii)) ||
      (length(term) == 0L && d$.time[ii[length(ii)]] != horizon)
  }, logical(1))
  if (any(terminal_bad))
    stop("Each subject must end at event/censoring or at the horizon.")

  if (!is.null(initiation_col)) {
    d$.init <- as.integer(d[[initiation_col]])
    if (any(!d$.init %in% 0:1))
      stop("initiation_col must be binary.")
    n_init <- vapply(split(seq_len(nrow(d)), d$.uid),
                     function(ii) sum(d$.init[ii]), numeric(1))
    if (any(n_init > 1L))
      stop("initiation_col may equal 1 at most once per subject.")
  } else if (!is.null(treatment_col)) {
    d$.trt <- as.integer(d[[treatment_col]])
    if (any(!d$.trt %in% 0:1))
      stop("treatment_col must be binary.")
    nonabsorbing <- vapply(split(seq_len(nrow(d)), d$.uid), function(ii) {
      any(diff(d$.trt[ii]) < 0)
    }, logical(1))
    if (any(nonabsorbing))
      stop("treatment_col must be absorbing (it cannot return from 1 to 0).")
    d$.init <- unlist(lapply(split(d$.trt, d$.uid),
                             function(x) as.integer(x == 1L &
                               c(0L, head(x, -1L)) == 0L)),
                      use.names = FALSE)
  } else {
    stop("Supply initiation_col or treatment_col.")
  }

  S_map <- vapply(split(seq_len(nrow(d)), d$.uid), function(ii) {
    z <- d$.time[ii][d$.init[ii] == 1L]
    if (length(z)) z[1] else Inf
  }, numeric(1))
  d$.S <- unname(S_map[as.character(d$.uid)])
  d$.trt <- as.integer(d$.time >= d$.S)
  d$.trt_lag <- unlist(lapply(split(d$.trt, d$.uid),
                              function(x) c(0L, head(x, -1L))),
                       use.names = FALSE)
  d$.time_f <- factor(d$.time, levels = seq_len(horizon))

  baseline_bad <- vapply(baseline_covariates, function(v) {
    any(vapply(split(d[[v]], d$.uid),
               function(x) length(unique(x)) != 1L, logical(1)))
  }, logical(1))
  if (any(baseline_bad))
    stop("Baseline covariates vary within subject: ",
         paste(baseline_covariates[baseline_bad], collapse = ", "))

  lag_names <- paste0("lag_", time_varying_covariates)
  for (j in seq_along(time_varying_covariates)) {
    v <- time_varying_covariates[j]
    d[[lag_names[j]]] <- unlist(lapply(split(d[[v]], d$.uid),
                                       function(x) c(lag_initial,
                                                     head(x, -1L))),
                                use.names = FALSE)
  }

  list(data = d, lag_names = lag_names)
}

g6_default_terms <- function(baseline, tv, lags, include_lags,
                             transition_include_treatment,
                             transition_include_time) {
  history <- c(baseline, tv, if (include_lags) lags)
  trans <- c(baseline, lags,
             if (transition_include_treatment) ".trt_lag",
             if (transition_include_time) ".time_f")
  list(
    initiation_denominator = c(".time_f", history),
    initiation_numerator = ".time_f",
    censoring_denominator = c(".time_f", history, ".trt"),
    censoring_numerator = ".time_f",
    outcome = c(".time_f", history, ".trt"),
    transition = trans
  )
}

g6_merge_terms <- function(defaults, user) {
  if (is.null(user)) return(defaults)
  bad <- setdiff(names(user), names(defaults))
  if (length(bad))
    stop("Unknown model_terms entries: ", paste(bad, collapse = ", "))
  for (nm in names(user)) defaults[[nm]] <- user[[nm]]
  defaults
}

g6_fit_binary <- function(data, response, terms, learner) {
  form <- g6_make_formula(response, terms)
  observed_response <- data[[response]]
  if (length(unique(observed_response)) < 2L) {
    probability <- g6_clamp_prob(mean(observed_response))
    constant_learner <- list(
      fit = function(formula, data) list(probability = probability),
      predict = function(model, newdata) {
        rep(model$probability, nrow(newdata))
      }
    )
    return(list(
      model = constant_learner$fit(form, data),
      learner = constant_learner,
      formula = form,
      constant = TRUE
    ))
  }
  list(model = learner$fit(form, data), learner = learner, formula = form)
}

g6_predict_binary <- function(fit, newdata) {
  g6_clamp_prob(fit$learner$predict(fit$model, newdata))
}

g6_common_init_numerator <- function(d, horizon) {
  ad <- d[d$.time <= pmin(d$.S, horizon), , drop = FALSE]
  out <- numeric(horizon)
  for (m in seq_len(horizon)) {
    rows <- ad$.time == m
    out[m] <- if (any(rows)) mean(ad$.init[rows]) else 0
  }
  g6_clamp_prob(out)
}

g6_local_components <- function(d, cfg, common_init_num = NULL) {
  M <- cfg$horizon
  units <- g6_first_by_unit(d)
  uid <- as.character(units$.uid)
  n <- length(uid)
  ii <- match(as.character(d$.uid), uid)
  mm <- d$.time
  cell <- cbind(ii, mm)

  obs <- event <- censor <- matrix(FALSE, n, M)
  obs[cell] <- TRUE
  event[cell] <- d$.event == 1L
  censor[cell] <- d$.censor == 1L

  init_rows <- d$.time <= pmin(d$.S, M)
  ad <- d[init_rows, , drop = FALSE]
  afit <- g6_fit_binary(ad, ".init", cfg$terms$initiation_denominator,
                        cfg$learners$initiation)
  Hden <- Hnum <- matrix(NA_real_, n, M)
  acell <- cbind(match(as.character(ad$.uid), uid), ad$.time)
  Hden[acell] <- g6_predict_binary(afit, ad)
  if (is.null(common_init_num)) {
    nfit <- g6_fit_binary(ad, ".init", cfg$terms$initiation_numerator,
                          cfg$learners$initiation)
    Hnum[acell] <- g6_predict_binary(nfit, ad)
  } else {
    if (length(common_init_num) != M || any(!is.finite(common_init_num))) {
      stop("common_init_num must contain one finite probability per interval.")
    }
    Hnum[acell] <- g6_clamp_prob(common_init_num[ad$.time])
  }

  cd <- d[d$.event == 0L, , drop = FALSE]
  cfit <- g6_fit_binary(cd, ".censor", cfg$terms$censoring_denominator,
                        cfg$learners$censoring)
  cnfit <- g6_fit_binary(cd, ".censor", cfg$terms$censoring_numerator,
                         cfg$learners$censoring)
  cd1 <- cd
  Son1 <- ifelse(cd$.S <= cfg$tau, cd$.S, cfg$tau)
  cd1$.trt <- as.integer(cd1$.time >= Son1)
  cd0 <- cd
  cd0$.trt <- 0L
  Cden1 <- Cden0 <- Cnum <- matrix(NA_real_, n, M)
  ccell <- cbind(match(as.character(cd$.uid), uid), cd$.time)
  Cden1[ccell] <- g6_predict_binary(cfit, cd1)
  Cden0[ccell] <- g6_predict_binary(cfit, cd0)
  Cnum[ccell] <- g6_predict_binary(cnfit, cd)

  S <- units$.S
  C1 <- ifelse(S <= cfg$tau, M + 1, cfg$tau - 1L)
  C0 <- ifelse(is.finite(S) & S <= M, S - 1L, M + 1)
  W1 <- W0 <- matrix(1, n, M)

  for (i in seq_len(n)) {
    aw1 <- aw0 <- nw1 <- nw0 <- 1
    for (m in seq_len(M)) {
      hd <- Hden[i, m]; hn <- Hnum[i, m]
      if (is.finite(hd) && is.finite(hn)) {
        if (m <= cfg$tau) {
          if (S[i] == m) aw1 <- aw1 * hn / hd
          if (S[i] > m) aw1 <- aw1 * (1 - hn) / (1 - hd)
        }
        if (S[i] > m) aw0 <- aw0 * (1 - hn) / (1 - hd)
      }

      W1[i, m] <- aw1 * nw1
      W0[i, m] <- aw0 * nw0

      if (obs[i, m] && !event[i, m] && !censor[i, m] &&
          is.finite(Cden1[i, m]) && is.finite(Cden0[i, m]) &&
          is.finite(Cnum[i, m])) {
        nw1 <- nw1 * (1 - Cnum[i, m]) / (1 - Cden1[i, m])
        nw0 <- nw0 * (1 - Cnum[i, m]) / (1 - Cden0[i, m])
      }
    }
  }

  clip <- function(W) {
    vals <- W[obs & is.finite(W)]
    q <- quantile(vals, probs = cfg$trunc, names = FALSE, na.rm = TRUE)
    pmin(pmax(W, q[1]), q[2])
  }
  W1 <- clip(W1); W0 <- clip(W0)

  R1 <- R0 <- D1 <- D0 <- matrix(FALSE, n, M)
  for (m in seq_len(M)) {
    R1[, m] <- obs[, m] & m <= C1
    R0[, m] <- obs[, m] & m <= C0
    D1[, m] <- R1[, m] & event[, m]
    D0[, m] <- R0[, m] & event[, m]
  }

  list(n = n, M = M, uid = uid, S = S, obs = obs, event = event,
       censor = censor, W1 = W1, W0 = W0, R1 = R1, R0 = R0,
       D1 = D1, D0 = D0)
}

g6_local_stats <- function(x) {
  list(
    M = x$M,
    dw1 = colSums(x$W1 * x$D1), rw1 = colSums(x$W1 * x$R1),
    dw0 = colSums(x$W0 * x$D0), rw0 = colSums(x$W0 * x$R0),
    local = x
  )
}

g6_contrasts <- function(S1, S0, M) {
  psi1 <- 1 - S1[M + 1L]
  psi0 <- 1 - S0[M + 1L]
  RMST1 <- sum(S1[seq_len(M)])
  RMST0 <- sum(S0[seq_len(M)])
  list(S1 = S1, S0 = S0, psi1 = psi1, psi0 = psi0,
       RD = psi1 - psi0, RR = psi1 / psi0,
       OR = (psi1 / (1 - psi1)) / (psi0 / (1 - psi0)),
       RMST1 = RMST1, RMST0 = RMST0, RMST_diff = RMST1 - RMST0)
}

g6_local_influence <- function(x, lambda1, lambda0, S1, S0, RW1, RW0,
                               est) {
  n <- x$n; M <- x$M
  cumU1 <- cumU0 <- numeric(n)
  IF_RMST1 <- IF_RMST0 <- numeric(n)
  for (m in seq_len(M)) {
    IF_RMST1 <- IF_RMST1 - S1[m] * cumU1
    IF_RMST0 <- IF_RMST0 - S0[m] * cumU0
    if (RW1[m] > 0 && lambda1[m] < 1)
      cumU1 <- cumU1 + (x$W1[, m] / RW1[m]) *
        (x$D1[, m] - lambda1[m] * x$R1[, m]) / (1 - lambda1[m])
    if (RW0[m] > 0 && lambda0[m] < 1)
      cumU0 <- cumU0 + (x$W0[, m] / RW0[m]) *
        (x$D0[, m] - lambda0[m] * x$R0[, m]) / (1 - lambda0[m])
  }
  IF1 <- S1[M + 1L] * cumU1
  IF0 <- S0[M + 1L] * cumU0
  p1 <- min(max(est$psi1, 1e-8), 1 - 1e-8)
  p0 <- min(max(est$psi0, 1e-8), 1 - 1e-8)
  list(
    sum_RD2 = sum((IF1 - IF0)^2),
    sum_logRR2 = sum((IF1 / p1 - IF0 / p0)^2),
    sum_logOR2 = sum((IF1 / (p1 * (1 - p1)) -
                       IF0 / (p0 * (1 - p0)))^2),
    sum_RMSTdiff2 = sum((IF_RMST1 - IF_RMST0)^2),
    sum_RMST1_2 = sum(IF_RMST1^2),
    sum_RMST0_2 = sum(IF_RMST0^2)
  )
}

g6_central_ccw <- function(stats) {
  M <- stats[[1]]$M
  add <- function(nm) Reduce(`+`, lapply(stats, `[[`, nm))
  DW1 <- add("dw1"); RW1 <- add("rw1")
  DW0 <- add("dw0"); RW0 <- add("rw0")
  lambda1 <- ifelse(RW1 > 0, DW1 / RW1, 0)
  lambda0 <- ifelse(RW0 > 0, DW0 / RW0, 0)
  S1 <- c(1, cumprod(1 - lambda1))
  S0 <- c(1, cumprod(1 - lambda0))
  est <- g6_contrasts(S1, S0, M)
  infl <- lapply(stats, function(s)
    g6_local_influence(s$local, lambda1, lambda0, S1, S0, RW1, RW0, est))
  agg <- function(nm) sqrt(sum(vapply(infl, `[[`, numeric(1), nm)))
  z <- qnorm(0.975)
  se_rd <- agg("sum_RD2"); se_rr <- agg("sum_logRR2")
  se_or <- agg("sum_logOR2"); se_rm <- agg("sum_RMSTdiff2")
  c(est, list(
    lambda1 = lambda1, lambda0 = lambda0,
    SE_RD = se_rd, SE_logRR = se_rr, SE_logOR = se_or,
    SE_RMSTdiff = se_rm,
    SE_RMST1 = agg("sum_RMST1_2"), SE_RMST0 = agg("sum_RMST0_2"),
    CI_RD = est$RD + c(-1, 1) * z * se_rd,
    CI_RR = exp(log(est$RR) + c(-1, 1) * z * se_rr),
    CI_OR = exp(log(est$OR) + c(-1, 1) * z * se_or),
    CI_RMSTdiff = est$RMST_diff + c(-1, 1) * z * se_rm,
    variance_type = paste(
      "federated influence function; fitted weights, estimated common",
      "numerator, and adaptive truncation cutoffs treated as fixed"
    )
  ))
}

g6_run_fed_ccw <- function(d, cfg) {
  site_data <- split(d, d$.site)
  stats <- lapply(site_data, function(x)
    g6_local_stats(g6_local_components(x, cfg)))
  g6_central_ccw(stats)
}

g6_run_pooled_ccw <- function(d, cfg) {
  g6_central_ccw(list(g6_local_stats(g6_local_components(d, cfg))))
}

g6_run_ipcw_no_clone <- function(d, cfg) {
  components <- lapply(split(d, d$.site), g6_local_components, cfg = cfg)
  M <- cfg$horizon
  DW1 <- RW1 <- DW0 <- RW0 <- numeric(M)
  for (x in components) {
    grp1 <- x$S <= cfg$tau
    for (m in seq_len(M)) {
      DW1[m] <- DW1[m] + sum(x$W1[grp1, m] *
                               x$event[grp1, m] * x$obs[grp1, m])
      RW1[m] <- RW1[m] + sum(x$W1[grp1, m] * x$obs[grp1, m])
      DW0[m] <- DW0[m] + sum(x$W0[!grp1, m] *
                               x$event[!grp1, m] * x$obs[!grp1, m])
      RW0[m] <- RW0[m] + sum(x$W0[!grp1, m] * x$obs[!grp1, m])
    }
  }
  l1 <- ifelse(RW1 > 0, DW1 / RW1, 0)
  l0 <- ifelse(RW0 > 0, DW0 / RW0, 0)
  g6_contrasts(c(1, cumprod(1 - l1)), c(1, cumprod(1 - l0)), M)
}

g6_run_perprotocol <- function(d, cfg) {
  M <- cfg$horizon
  units <- g6_first_by_unit(d)
  Smap <- setNames(units$.S, as.character(units$.uid))
  d$.group1 <- Smap[as.character(d$.uid)] <= cfg$tau
  curve <- function(g) {
    z <- d[d$.group1 == g, , drop = FALSE]
    lam <- numeric(M)
    for (m in seq_len(M)) {
      r <- z$.time == m
      lam[m] <- if (any(r)) sum(z$.event[r]) / sum(r) else 0
    }
    c(1, cumprod(1 - lam))
  }
  g6_contrasts(curve(TRUE), curve(FALSE), M)
}

g6_make_plain_frame <- function(d) d

g6_make_cloned_frame <- function(d, cfg) {
  M <- cfg$horizon
  S <- d$.S
  C1 <- ifelse(S <= cfg$tau, M + 1, cfg$tau - 1L)
  C0 <- ifelse(is.finite(S) & S <= M, S - 1L, M + 1)
  d1 <- d[d$.time <= C1, , drop = FALSE]
  Son1 <- ifelse(d1$.S <= cfg$tau, d1$.S, cfg$tau)
  d1$.trt <- as.integer(d1$.time >= Son1)
  d1$.trt_lag <- as.integer((d1$.time - 1L) >= Son1)
  d1$.event <- as.integer(d1$.event == 1L & d1$.time <=
                            ifelse(d1$.S <= cfg$tau, M + 1,
                                   cfg$tau - 1L))
  d1$.arm <- 1L
  d0 <- d[d$.time <= C0, , drop = FALSE]
  d0$.trt <- 0L
  d0$.trt_lag <- 0L
  d0$.event <- as.integer(d0$.event == 1L & d0$.time <=
                            ifelse(is.finite(d0$.S) & d0$.S <= M,
                                   d0$.S - 1L, M + 1))
  d0$.arm <- 0L
  rbind(d1, d0)
}

g6_fit_transitions <- function(frame, cfg) {
  if (!length(cfg$tv)) return(list())
  tr <- frame[frame$.time >= 2L, , drop = FALSE]
  out <- vector("list", length(cfg$tv)); names(out) <- cfg$tv
  for (v in cfg$tv) {
    fam <- cfg$tv_families[[v]]
    form <- g6_make_formula(v, cfg$terms$transition)
    if (fam == "binomial") {
      fit <- glm(form, data = tr, family = binomial())
      out[[v]] <- list(family = fam, model = fit, sd = NA_real_)
    } else {
      fit <- lm(form, data = tr)
      out[[v]] <- list(family = fam, model = fit,
                       sd = sd(residuals(fit), na.rm = TRUE))
    }
  }
  out
}

g6_fit_gcomp <- function(frame, cfg) {
  list(
    outcome = g6_fit_binary(frame, ".event", cfg$terms$outcome,
                            cfg$learners$outcome),
    transitions = g6_fit_transitions(frame, cfg)
  )
}

g6_gcomp_standardize <- function(d, fit, cfg, seed) {
  set.seed(seed)
  M <- cfg$horizon; mc <- cfg$mc
  units <- d[d$.time == 1L, , drop = FALSE]
  units <- units[!duplicated(units$.uid), , drop = FALSE]
  n <- nrow(units)
  draw_id <- rep(seq_len(n), each = mc)
  nd_base <- units[draw_id, cfg$baseline, drop = FALSE]
  S1 <- ifelse(units$.S <= cfg$tau, units$.S, cfg$tau)
  S0 <- rep(Inf, n)

  simulate_arm <- function(Sarm, arm_seed) {
    set.seed(arm_seed)
    Son <- Sarm[draw_id]
    surv <- rep(1, length(draw_id))
    curve <- numeric(M + 1L); curve[1] <- 1
    prev <- matrix(cfg$lag_initial, nrow = length(draw_id),
                   ncol = length(cfg$tv),
                   dimnames = list(NULL, cfg$lag_names))

    for (m in seq_len(M)) {
      trt <- as.integer(m >= Son)
      if (m == 1L) {
        current <- as.data.frame(units[draw_id, cfg$tv, drop = FALSE])
      } else {
        trans_nd <- nd_base
        for (j in seq_along(cfg$lag_names))
          trans_nd[[cfg$lag_names[j]]] <- prev[, j]
        trans_nd$.trt_lag <- as.integer((m - 1L) >= Son)
        trans_nd$.time_f <- factor(rep(m, nrow(trans_nd)),
                                   levels = seq_len(M))
        current <- vector("list", length(cfg$tv)); names(current) <- cfg$tv
        for (v in cfg$tv) {
          tf <- fit$transitions[[v]]
          mu <- if (tf$family == "binomial")
            predict(tf$model, newdata = trans_nd, type = "response") else
            predict(tf$model, newdata = trans_nd)
          current[[v]] <- if (tf$family == "binomial")
            rbinom(length(mu), 1L, g6_clamp_prob(mu)) else
            rnorm(length(mu), mu, tf$sd)
        }
        current <- as.data.frame(current)
      }

      out_nd <- nd_base
      for (v in cfg$tv) out_nd[[v]] <- current[[v]]
      for (j in seq_along(cfg$lag_names))
        out_nd[[cfg$lag_names[j]]] <- prev[, j]
      out_nd$.trt <- trt
      out_nd$.time_f <- factor(rep(m, nrow(out_nd)), levels = seq_len(M))
      haz <- g6_predict_binary(fit$outcome, out_nd)
      surv <- surv * (1 - haz)
      curve[m + 1L] <- mean(surv)
      if (length(cfg$tv))
        prev[,] <- as.matrix(current[, cfg$tv, drop = FALSE])
    }
    curve
  }

  curve1 <- simulate_arm(S1, seed)
  curve0 <- simulate_arm(S0, seed + 1L)
  g6_contrasts(curve1, curve0, M)
}

g6_run_cc_gcomp <- function(d, cfg, seed) {
  frame <- g6_make_cloned_frame(d, cfg)
  fit <- g6_fit_gcomp(frame, cfg)
  g6_gcomp_standardize(d, fit, cfg, seed)
}

g6_run_direct_gcomp <- function(d, cfg, seed) {
  frame <- g6_make_plain_frame(d)
  fit <- g6_fit_gcomp(frame, cfg)
  g6_gcomp_standardize(d, fit, cfg, seed)
}

g6_result_rows <- function(fit, method) {
  estimands <- c("psi1", "psi0", "RD", "RR", "OR",
                 "RMST1", "RMST0", "RMST_diff")
  se_map <- c(RD = "SE_RD", RR = "SE_logRR", OR = "SE_logOR",
              RMST1 = "SE_RMST1", RMST0 = "SE_RMST0",
              RMST_diff = "SE_RMSTdiff")
  ci_map <- c(RD = "CI_RD", RR = "CI_RR", OR = "CI_OR",
              RMST_diff = "CI_RMSTdiff")
  do.call(rbind, lapply(estimands, function(e) {
    se_name <- unname(se_map[e])
    ci_name <- unname(ci_map[e])
    se <- if (is.na(se_name)) NA_real_ else fit[[se_name]] %||% NA_real_
    ci <- if (is.na(ci_name)) c(NA_real_, NA_real_) else
      fit[[ci_name]] %||% c(NA_real_, NA_real_)
    data.frame(method = method, estimand = e, estimate = fit[[e]],
               se = se, ci_lo = ci[1], ci_hi = ci[2],
               stringsAsFactors = FALSE)
  }))
}

#' Run six grace-period methods with natural censoring
#'
#' @return A list with `estimates`, method-specific `fits`, processed `data`,
#'   and the resolved `config`.
run_six_methods_natural_censoring <- function(
    data,
    id_col,
    time_col,
    event_col,
    grace_period,
    baseline_covariates = character(),
    time_varying_covariates = character(),
    site_col = NULL,
    censor_col = NULL,
    initiation_col = NULL,
    treatment_col = NULL,
    horizon = NULL,
    weight_truncation = c(0.01, 0.99),
    mc_gcomp = 20,
    seed = 2026,
    tv_families = NULL,
    include_lagged_tv = TRUE,
    transition_include_treatment = FALSE,
    transition_include_time = FALSE,
    lag_initial = 0,
    model_terms = NULL,
    learners = NULL,
    verbose = TRUE) {

  horizon <- horizon %||% max(as.integer(data[[time_col]]), na.rm = TRUE)
  horizon <- as.integer(horizon)
  tau <- as.integer(grace_period)
  if (length(tau) != 1L || is.na(tau) || tau < 1L || tau > horizon)
    stop("grace_period must be one integer in 1:horizon.")
  if (length(weight_truncation) != 2L ||
      any(weight_truncation < 0 | weight_truncation > 1) ||
      weight_truncation[1] > weight_truncation[2])
    stop("weight_truncation must be two ordered probabilities.")
  if (mc_gcomp < 1) stop("mc_gcomp must be positive.")

  prepared <- g6_prepare_data(
    data, site_col, id_col, time_col, event_col, censor_col,
    initiation_col, treatment_col, baseline_covariates,
    time_varying_covariates, horizon, lag_initial
  )
  d <- prepared$data
  lag_names <- prepared$lag_names

  if (is.null(tv_families))
    tv_families <- setNames(rep("gaussian", length(time_varying_covariates)),
                            time_varying_covariates)
  tv_families <- as.list(tv_families)
  missing_fam <- setdiff(time_varying_covariates, names(tv_families))
  if (length(missing_fam))
    stop("tv_families is missing: ", paste(missing_fam, collapse = ", "))
  bad_fam <- !vapply(tv_families[time_varying_covariates],
                     function(x) x %in% c("gaussian", "binomial"), logical(1))
  if (any(bad_fam))
    stop("tv_families values must be 'gaussian' or 'binomial'.")
  for (v in time_varying_covariates) {
    if (!is.numeric(d[[v]]))
      stop("Time-varying covariate '", v, "' must be numeric.")
    if (tv_families[[v]] == "binomial" && any(!d[[v]] %in% 0:1))
      stop("Binomial time-varying covariate '", v, "' must contain only 0/1.")
  }

  default_learner <- g6_default_binary_learner()
  learners <- learners %||% list()
  resolved_learners <- list(
    initiation = g6_validate_learner(learners$initiation %||%
                                       default_learner, "initiation"),
    censoring = g6_validate_learner(learners$censoring %||%
                                      default_learner, "censoring"),
    outcome = g6_validate_learner(learners$outcome %||%
                                    default_learner, "outcome")
  )
  terms <- g6_default_terms(
    baseline_covariates, time_varying_covariates, lag_names,
    include_lagged_tv, transition_include_treatment,
    transition_include_time
  )
  terms <- g6_merge_terms(terms, model_terms)
  unknown_terms <- setdiff(unique(unlist(terms)),
                           c(names(d), ".time_f", ".trt", ".trt_lag",
                             lag_names))
  if (length(unknown_terms))
    stop("Unknown variables in model_terms: ",
         paste(unknown_terms, collapse = ", "))

  cfg <- list(
    tau = tau, horizon = horizon, baseline = baseline_covariates,
    tv = time_varying_covariates, lag_names = lag_names,
    tv_families = tv_families, terms = terms,
    learners = resolved_learners, trunc = weight_truncation,
    mc = as.integer(mc_gcomp), lag_initial = lag_initial
  )

  if (verbose) {
    message(sprintf("Running six methods: N=%d, sites=%d, tau=%d, horizon=%d",
                    length(unique(d$.uid)), length(unique(d$.site)),
                    tau, horizon))
  }
  fits <- list(
    fed_ccw = g6_run_fed_ccw(d, cfg),
    pooled_ccw = g6_run_pooled_ccw(d, cfg),
    direct_gcomp = g6_run_direct_gcomp(d, cfg, seed + 100L),
    cc_gcomp = g6_run_cc_gcomp(d, cfg, seed + 100L),
    ipcw_no_clone = g6_run_ipcw_no_clone(d, cfg),
    perprotocol_naive = g6_run_perprotocol(d, cfg)
  )
  estimates <- do.call(rbind, Map(g6_result_rows, fits, names(fits)))
  rownames(estimates) <- NULL
  invisible(list(estimates = estimates, fits = fits, data = d, config = cfg))
}

# =============================================================================
# Hopkins pneumonia study's five-method set
#
# This matches ArraySim5-6 Uncertainty:
#   1. federated CCW;
#   2. pooled CCW with site-stratified nuisance models (identity control);
#   3. federated landmark IPW;
#   4. federated IPW without cloning;
#   5. federated unweighted per-protocol analysis.
#
# The landmark analysis is conditional on being alive and observed at tau and
# therefore has a different target population from the other five methods.
# =============================================================================

g6_grouped_stats_from_component <- function(x, tau, weighted = TRUE) {
  group1 <- x$S <= tau
  arm1 <- matrix(group1, nrow = x$n, ncol = x$M)
  arm0 <- !arm1
  local <- list(
    n = x$n,
    M = x$M,
    W1 = if (weighted) x$W1 else matrix(1, x$n, x$M),
    W0 = if (weighted) x$W0 else matrix(1, x$n, x$M),
    R1 = x$obs & arm1,
    R0 = x$obs & arm0,
    D1 = x$obs & x$event & arm1,
    D0 = x$obs & x$event & arm0
  )
  g6_local_stats(local)
}

g6_ipcw_no_clone_from_components <- function(components, cfg) {
  stats <- lapply(
    components,
    g6_grouped_stats_from_component,
    tau = cfg$tau,
    weighted = TRUE
  )
  g6_central_ccw(stats)
}

g6_perprotocol_from_components <- function(components, cfg) {
  stats <- lapply(
    components,
    g6_grouped_stats_from_component,
    tau = cfg$tau,
    weighted = FALSE
  )
  g6_central_ccw(stats)
}

g6_landmark_local_stats <- function(site_data, cfg, p_num) {
  tau <- cfg$tau
  M <- cfg$horizon
  landmark <- site_data[
    site_data$.time == tau & site_data$.event == 0L & site_data$.censor == 0L,
    ,
    drop = FALSE
  ]
  if (!nrow(landmark)) {
    return(NULL)
  }

  landmark$.A <- as.integer(landmark$.S <= tau)
  propensity_terms <- setdiff(
    cfg$terms$initiation_denominator,
    c(".time_f", ".trt", ".trt_lag")
  )
  if (length(unique(landmark$.A)) < 2L) {
    p_den <- rep(g6_clamp_prob(mean(landmark$.A)), nrow(landmark))
  } else {
    fit <- g6_fit_binary(
      landmark,
      ".A",
      propensity_terms,
      cfg$learners$initiation
    )
    p_den <- g6_predict_binary(fit, landmark)
    if (any(!is.finite(p_den))) {
      p_den <- rep(g6_clamp_prob(mean(landmark$.A)), nrow(landmark))
    }
  }

  weights <- ifelse(
    landmark$.A == 1L,
    p_num / p_den,
    (1 - p_num) / (1 - p_den)
  )
  limits <- quantile(
    weights[is.finite(weights)],
    probs = cfg$trunc,
    names = FALSE,
    na.rm = TRUE
  )
  weights <- pmin(pmax(weights, limits[1L]), limits[2L])
  uid <- as.character(landmark$.uid)
  names(weights) <- uid
  n <- length(uid)
  obs <- event <- matrix(FALSE, n, M)
  rows <- site_data$.uid %in% landmark$.uid & site_data$.time > tau
  post <- site_data[rows, , drop = FALSE]
  if (nrow(post)) {
    cell <- cbind(match(as.character(post$.uid), uid), post$.time)
    obs[cell] <- TRUE
    event[cell] <- post$.event == 1L
  }
  arm1 <- matrix(landmark$.A == 1L, nrow = n, ncol = M)
  post_time <- matrix(rep(seq_len(M) > tau, each = n), nrow = n)
  local <- list(
    n = n,
    M = M,
    W1 = matrix(weights, nrow = n, ncol = M),
    W0 = matrix(weights, nrow = n, ncol = M),
    R1 = obs & arm1 & post_time,
    R0 = obs & !arm1 & post_time,
    D1 = obs & event & arm1 & post_time,
    D0 = obs & event & !arm1 & post_time
  )
  out <- g6_local_stats(local)
  out$eligible <- n
  out$treated <- sum(landmark$.A)
  out
}

g6_run_landmark_ipw_general <- function(d, cfg) {
  site_data <- split(d, d$.site)
  counts <- vapply(site_data, function(x) {
    landmark <- x[
      x$.time == cfg$tau & x$.event == 0L & x$.censor == 0L,
      ,
      drop = FALSE
    ]
    c(
      treated = sum(landmark$.S <= cfg$tau),
      eligible = nrow(landmark)
    )
  }, numeric(2))
  totals <- rowSums(counts)
  if (totals[["eligible"]] <= 0) {
    stop("No participants are alive and observed at the landmark.")
  }
  p_num <- g6_clamp_prob(totals[["treated"]] / totals[["eligible"]])
  stats <- Filter(
    Negate(is.null),
    lapply(site_data, g6_landmark_local_stats, cfg = cfg, p_num = p_num)
  )
  fit <- g6_central_ccw(stats)
  fit$landmark_eligible <- as.integer(totals[["eligible"]])
  fit$landmark_treated <- as.integer(totals[["treated"]])
  fit
}

#' Run the five estimators used for the Hopkins pneumonia comparison
#'
#' Input is long person-interval data. The treatment strategy for arm 1 is
#' initiation by `grace_period`; arm 0 is no initiation through the horizon.
#' Event/censor rows must be terminal. Time-varying covariates must be numeric.
run_pneumonia_five_methods <- function(
    data,
    id_col,
    time_col,
    event_col,
    grace_period,
    baseline_covariates = character(),
    time_varying_covariates = character(),
    site_col = NULL,
    censor_col = NULL,
    initiation_col = NULL,
    treatment_col = NULL,
    horizon = NULL,
    weight_truncation = c(0.01, 0.99),
    include_lagged_tv = TRUE,
    lag_initial = 0,
    model_terms = NULL,
    learners = NULL,
    verbose = TRUE) {

  horizon <- as.integer(horizon %||% max(data[[time_col]], na.rm = TRUE))
  tau <- as.integer(grace_period)
  if (length(tau) != 1L || is.na(tau) || tau < 1L || tau >= horizon) {
    stop("grace_period must be one integer in 1:(horizon - 1).")
  }

  prepared <- g6_prepare_data(
    data = data,
    site_col = site_col,
    id_col = id_col,
    time_col = time_col,
    event_col = event_col,
    censor_col = censor_col,
    initiation_col = initiation_col,
    treatment_col = treatment_col,
    baseline_covariates = baseline_covariates,
    time_varying_covariates = time_varying_covariates,
    horizon = horizon,
    lag_initial = lag_initial
  )
  d <- prepared$data
  lag_names <- prepared$lag_names

  default_learner <- g6_default_binary_learner()
  learners <- learners %||% list()
  resolved_learners <- list(
    initiation = g6_validate_learner(
      learners$initiation %||% default_learner,
      "initiation"
    ),
    censoring = g6_validate_learner(
      learners$censoring %||% default_learner,
      "censoring"
    ),
    outcome = g6_validate_learner(
      learners$outcome %||% default_learner,
      "outcome"
    )
  )
  terms <- g6_default_terms(
    baseline_covariates,
    time_varying_covariates,
    lag_names,
    include_lagged_tv,
    transition_include_treatment = FALSE,
    transition_include_time = FALSE
  )
  terms <- g6_merge_terms(terms, model_terms)
  cfg <- list(
    tau = tau,
    horizon = horizon,
    baseline = baseline_covariates,
    tv = time_varying_covariates,
    lag_names = lag_names,
    terms = terms,
    learners = resolved_learners,
    trunc = weight_truncation,
    lag_initial = lag_initial
  )

  if (verbose) {
    message(sprintf(
      "Running Hopkins pneumonia five-method set: N=%d, sites=%d, tau=%d, horizon=%d",
      length(unique(d$.uid)),
      length(unique(d$.site)),
      tau,
      horizon
    ))
  }

  common_init_num <- g6_common_init_numerator(d, horizon)
  site_components <- lapply(
    split(d, d$.site),
    g6_local_components,
    cfg = cfg,
    common_init_num = common_init_num
  )
  site_stats <- lapply(site_components, g6_local_stats)
  fits <- list(
    fed_ccw = g6_central_ccw(site_stats),
    pooled_ccw_site_stratified = g6_central_ccw(site_stats),
    landmark_ipw = g6_run_landmark_ipw_general(d, cfg),
    ipcw_no_clone = g6_ipcw_no_clone_from_components(site_components, cfg),
    perprotocol_naive = g6_perprotocol_from_components(site_components, cfg)
  )
  identity_fields <- c(
    "psi1", "psi0", "RD", "RR", "OR", "RMST1", "RMST0", "RMST_diff",
    "SE_RD", "SE_logRR", "SE_logOR", "SE_RMSTdiff", "SE_RMST1", "SE_RMST0"
  )
  identity_differences <- vapply(identity_fields, function(field) {
    abs(
      as.numeric(fits$fed_ccw[[field]]) -
        as.numeric(fits$pooled_ccw_site_stratified[[field]])
    )
  }, numeric(1))
  if (any(identity_differences > 1e-12, na.rm = TRUE)) {
    stop("Federated CCW failed the site-stratified pooled identity check.")
  }
  estimates <- do.call(rbind, Map(g6_result_rows, fits, names(fits)))
  rownames(estimates) <- NULL

  invisible(list(
    estimates = estimates,
    fits = fits,
    data = d,
    config = cfg,
    site_components = site_components,
    common_init_numerator = common_init_num,
    identity_check = data.frame(
      quantity = identity_fields,
      absolute_difference = identity_differences,
      tolerance = 1e-12,
      passed = identity_differences <= 1e-12,
      stringsAsFactors = FALSE
    )
  ))
}
