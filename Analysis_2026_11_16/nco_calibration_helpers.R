# =============================================================
# Summary-level NN and robust mixture Normal-Normal calibration.
# Adapted from ArraySim6-1 NCO/NCO_Calibration.R; exploratory small-panel use.
# Inputs are historical log prevalence ratios, not prospective outcome logRRs.
# =============================================================

.component_mle <- function(y, se, w, start_mu, start_sigma) {
  if (sum(w) < 1e-6)
    return(c(mu = start_mu, sigma = max(start_sigma, 1e-4)))
  objective <- function(par) {
    mu <- par[1]; sigma2 <- exp(par[2])
    v <- se^2 + sigma2
    0.5 * sum(w * (log(v) + (y - mu)^2 / v))
  }
  fit <- optim(c(start_mu, log(max(start_sigma^2, 1e-8))), objective,
               method = "BFGS", control = list(maxit = 200))
  if (fit$convergence != 0L || !is.finite(fit$value)) stop("Component optimizer did not converge")
  c(mu = fit$par[1], sigma = sqrt(exp(fit$par[2])))
}

fit_nn_calibration <- function(y, se) {
  keep <- is.finite(y) & is.finite(se) & se > 0
  y <- y[keep]; se <- se[keep]
  if (length(y) < 3L) stop("At least three finite NCOs are required.")
  start_sigma <- sqrt(max(var(y) - mean(se^2), 1e-6))
  p <- .component_mle(y, se, rep(1, length(y)), mean(y), start_sigma)
  v <- se^2 + p[["sigma"]]^2
  ll <- sum(dnorm(y, p[["mu"]], sqrt(v), log = TRUE))
  list(mu_valid = p[["mu"]], sigma_valid = p[["sigma"]],
       pi_valid = 1, logLik = ll, converged = TRUE, n_used = length(y),
       model = "NN")
}

.one_mnn_start <- function(y, se, pi0, mu_v0, mu_i0, sigma0,
                           tol = 1e-8, max_iter = 500L) {
  pi_v <- pi0; mu_v <- mu_v0; mu_i <- mu_i0
  sig_v <- sig_i <- max(sigma0, 1e-3)
  ll_old <- -Inf

  for (iter in seq_len(max_iter)) {
    lv <- log(pi_v) + dnorm(y, mu_v, sqrt(se^2 + sig_v^2), log = TRUE)
    li <- log1p(-pi_v) + dnorm(y, mu_i, sqrt(se^2 + sig_i^2), log = TRUE)
    mx <- pmax(lv, li)
    r <- exp(lv - mx) / (exp(lv - mx) + exp(li - mx))

    # The valid component is defined as the strict majority component.
    if (mean(r) < 0.5) r <- 1 - r
    pi_v <- min(max(mean(r), 0.5001), 0.9999)
    pv <- .component_mle(y, se, r, mu_v, sig_v)
    pi <- .component_mle(y, se, 1 - r, mu_i, sig_i)
    mu_v <- pv[["mu"]]; sig_v <- pv[["sigma"]]
    mu_i <- pi[["mu"]]; sig_i <- pi[["sigma"]]

    lv <- log(pi_v) + dnorm(y, mu_v, sqrt(se^2 + sig_v^2), log = TRUE)
    li <- log1p(-pi_v) + dnorm(y, mu_i, sqrt(se^2 + sig_i^2), log = TRUE)
    mx <- pmax(lv, li)
    ll <- sum(mx + log(exp(lv - mx) + exp(li - mx)))
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * (1 + abs(ll_old))) break
    ll_old <- ll
  }
  list(mu_valid = mu_v, sigma_valid = sig_v, pi_valid = pi_v,
       mu_invalid = mu_i, sigma_invalid = sig_i, logLik = ll,
       converged = iter < max_iter, iterations = iter)
}

fit_mnn_calibration <- function(y, se) {
  keep <- is.finite(y) & is.finite(se) & se > 0
  y <- y[keep]; se <- se[keep]
  # Four is a computational guard, not a scientifically adequate panel size.
  if (length(y) < 4L) stop("At least four distinct finite NCOs are required for exploratory MNN.")
  spread <- max(sd(y), 0.05)
  med <- median(y)
  starts <- expand.grid(
    pi0 = c(0.55, 0.70, 0.85, 0.95),
    direction = c(-1, 1), KEEP.OUT.ATTRS = FALSE
  )
  fits <- lapply(seq_len(nrow(starts)), function(k) {
    d <- starts$direction[k]
    try(.one_mnn_start(
      y, se, pi0 = starts$pi0[k],
      mu_v0 = med, mu_i0 = med + d * spread,
      sigma0 = max(spread / 2, 0.01)
    ), silent = TRUE)
  })
  ok <- !vapply(fits, inherits, logical(1), "try-error") &
    vapply(fits, function(x) if (inherits(x, "try-error")) FALSE else
      is.finite(x$logLik) && isTRUE(x$converged), logical(1))
  if (!any(ok)) stop("All MNN starts failed.")
  fits <- fits[ok]
  best <- fits[[which.max(vapply(fits, `[[`, numeric(1), "logLik"))]]
  best$n_used <- length(y); best$model <- "MNN"
  best
}

calibrate_primary <- function(primary_estimate, primary_se, fit,
                              method_label) {
  estimate <- primary_estimate - fit$mu_valid
  se <- sqrt(primary_se^2 + fit$sigma_valid^2)
  z <- qnorm(0.975)
  data.frame(
    method = method_label, estimate_logRR = estimate, se_logRR = se,
    ci_low = estimate - z * se, ci_high = estimate + z * se,
    p_value = 2 * pnorm(-abs(estimate / se)),
    mu_valid_hat = fit$mu_valid, sigma_valid_hat = fit$sigma_valid,
    pi_valid_hat = fit$pi_valid, calibration_converged = fit$converged,
    stringsAsFactors = FALSE
  )
}
