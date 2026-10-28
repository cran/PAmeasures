# Event-observation IPCW from PA measure new/pam.survreg.update.R.
# Working outcomes may be restricted, but censoring is estimated and evaluated
# at the ORIGINAL follow-up times. Only observed events of any cause contribute.
# These helpers are separate from the ordinary right-censoring estimators.

#' @noRd
pa_cr_event_ipcw <- function(ftime, fstatus, X = NULL,
                             censoring = c("cox", "km")) {
  censoring <- match.arg(censoring)
  n <- length(ftime)
  if (!is.numeric(ftime) || !n || length(fstatus) != n ||
      anyNA(ftime) || any(!is.finite(ftime)) || any(ftime < 0) ||
      !is.numeric(fstatus) || anyNA(fstatus) ||
      any(!is.finite(fstatus)) || any(fstatus < 0) ||
      any(fstatus != floor(fstatus))) {
    stop("Original follow-up must be finite and non-negative, with matching non-negative integer status codes.",
         call. = FALSE)
  }
  delta <- as.numeric(fstatus > 0)
  if (!is.null(X)) {
    if (is.atomic(X) && is.null(dim(X))) X <- matrix(X, ncol = 1L)
    X <- as.data.frame(X)
    if (nrow(X) != n || !ncol(X) || anyNA(X) ||
        any(!vapply(X, is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(X)))) {
      stop("'censor_covariates' must contain finite numeric columns with one row per subject.",
           call. = FALSE)
    }
    if (censoring == "km") {
      stop("'censor_covariates' require censoring = 'cox'.", call. = FALSE)
    }
    names(X) <- paste0("x", seq_len(ncol(X)))
  }

  if (all(delta == 1) || all(delta == 0)) {
    G <- rep(1, n)
  } else if (censoring == "km") {
    fit <- survival::survfit(survival::Surv(ftime, 1 - delta) ~ 1)
    # A common left limit at a tied time must not depend on subject row order.
    idx <- findInterval(ftime, fit$time, left.open = TRUE)
    G <- c(1, fit$surv)[idx + 1L]
  } else {
    dat <- data.frame(original_time = ftime, censored = 1 - delta)
    if (!is.null(X)) dat <- cbind(dat, X)
    fit <- survival::coxph(
      survival::Surv(original_time, censored) ~ ., data = dat,
      ties = "breslow")
    # Center the baseline and linear predictors together. The source script's
    # uncentered baseline / centered predictor combination is not reproduced.
    bh <- survival::basehaz(fit, centered = TRUE)
    risk <- exp(stats::predict(fit, type = "lp", reference = "sample"))
    idx <- findInterval(ftime, bh$time, left.open = TRUE)
    H0 <- c(0, bh$hazard)[idx + 1L]
    G <- exp(-H0 * risk)
  }

  # Retain the source strategy's explicit numerical probability floor.
  G <- pmax(G, 1e-4)
  if (any(!is.finite(G))) {
    stop("The censoring model produced non-finite survival probabilities.", call. = FALSE)
  }
  q <- delta / G
  weight <- if (sum(q) > 0) q / sum(q) else rep(0, n)
  list(probability = G, observed_event = delta, weight = weight)
}

#' @noRd
pa_cr_weighted_r2 <- function(outcome, prediction, weight) {
  missing_result <- list(R.squared = NA_real_, L.squared = NA_real_)
  if (sum(weight > 0) < 2L) {
    warning("Competing-risks R-squared requires at least two observed events.",
            call. = FALSE)
    return(missing_result)
  }
  fit <- stats::lm(outcome ~ prediction, weights = weight)
  corrected <- as.numeric(stats::predict(fit))
  outcome_mean <- sum(weight * outcome)
  total <- sum(weight * (outcome - outcome_mean)^2)
  original_error <- sum(weight * (outcome - prediction)^2)
  list(
    R.squared = if (total > 0)
      round(sum(weight * (corrected - outcome_mean)^2) / total, 4) else NA_real_,
    L.squared = if (original_error > 0)
      round(sum(weight * (outcome - corrected)^2) / original_error, 4) else NA_real_
  )
}
