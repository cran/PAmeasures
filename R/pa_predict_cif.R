
#' Predict Subject-Specific Cumulative Incidence Function (CIF)
#'
#' This function calculates cumulative incidence function (CIF) predictions for individual subjects 
#' across multiple types of competing risks models. It supports cause-specific Cox models, 
#' accelerated failure time (AFT) models (Weibull, Log-normal, Log-logistic), Fine-Gray models, 
#' and competing risks random survival forests.
#'
#' @param model1 The fitted model for the target cause (either a `coxph` or `survreg` object).
#' @param model2 The fitted model for competing causes (either a `coxph` or `survreg` object).
#' @param fg_model Optional Fine-Gray model object (`crr` from `cmprsk` package).
#' @param cr_model Optional competing risks random forest model object.
#' @param newdata A data frame containing observed event times (`time`), event indicators (`status`).
#' @param event.type The event type of interest for which CIF is to be predicted.
#' @param tau The maximum time horizon for prediction.
#' @param covs A character vector specifying the names of the covariates used in the model.
#' @details For cause-specific AFT models, CIFs are evaluated by adaptive
#' quadrature on the target cause's probability scale. This handles Weibull
#' densities that are unbounded at time zero, including fixed scales above 1.
#' The time summary integrates one minus the CIF on the returned time grid.
#' @return A list with the following components:
#' \describe{
#'   \item{times}{Numeric vector of observed times from \code{newdata$time}.}
#'   \item{status}{Integer vector of event indicators from \code{newdata$status}.}
#'   \item{cif_pred}{A matrix where the first column is time, and subsequent columns are predicted CIFs for each subject.}
#'   \item{pred}{Numeric vector of the integral of one minus the CIF from \code{0} to \code{tau}.}
#'   \item{linear.pred}{Numeric vector of subject-specific linear predictors.}
#' }
#' 
#' @examples
#' data(pbc, package = "survival")
#' pbc <- pbc[!is.na(pbc$trt), ]
#' pbc$log_albumin <- log(pbc$albumin)
#' pbc$log_bili    <- log(pbc$bili)
#' pbc$log_protime <- log(pbc$protime)
#' covariates <- c("age", "log_albumin", "log_bili", "log_protime", "edema")
#' dat <- pbc[, c("time", "status", covariates)]
#' dat <- dat[stats::complete.cases(dat), ]
#'
#' # Cause-specific Cox models: event of interest (1) and competing event (2)
#' m1 <- survival::coxph(survival::Surv(time, status == 1) ~ age + log_albumin +
#'                         log_bili + log_protime + edema,
#'                       data = dat, x = TRUE, y = TRUE)
#' m2 <- survival::coxph(survival::Surv(time, status == 2) ~ age + log_albumin +
#'                         log_bili + log_protime + edema,
#'                       data = dat, x = TRUE, y = TRUE)
#'
#' pred <- pa_predict_cif(model1 = m1, model2 = m2, newdata = dat,
#'                        covs = covariates, event.type = 1,
#'                        tau = max(dat$time))
#' str(pred$pred)
#'
#' # Fine-Gray backend, available when cmprsk is installed
#' if (requireNamespace("cmprsk", quietly = TRUE)) {
#'   fg <- cmprsk::crr(ftime = dat$time, fstatus = dat$status,
#'                     cov1 = as.matrix(dat[, covariates]),
#'                     failcode = 1, cencode = 0)
#'   pred_fg <- pa_predict_cif(fg_model = fg, newdata = dat,
#'                             covs = covariates, event.type = 1,
#'                             tau = max(dat$time))
#'   str(pred_fg$pred)
#' }
#'
#' @export


pa_predict_cif <- function(model1 = NULL, model2 = NULL,
                                    fg_model = NULL, cr_model = NULL,
                                    tau = NULL, newdata, 
                                    event.type = 1, covs) {
  
  #quantile_list <- seq(0, 1, length.out = num.grid)
  #LB_tau <- min(newdata$time[newdata$status == event.type])
  #UB_tau <- max(newdata$time[newdata$status == event.type])
  #tau_list <- qunif(quantile_list, min = LB_tau, max = min(UB_tau, tau))
  
  if(is.null(tau)) tau <- max(newdata$time)
  ###  Handle Fine-Gray model (crr model fg_model)
  if (!is.null(fg_model)) {
    if (!requireNamespace("cmprsk", quietly = TRUE)) {
      stop("The 'fg_model' argument requires the 'cmprsk' package. ",
           "Install it with install.packages(\"cmprsk\").", call. = FALSE)
    }
    m_pred <- predict(fg_model, as.matrix(newdata[, covs, drop = FALSE]))
    if (is.null(dim(m_pred[, -1]))) return(NULL)
    pred <- apply(m_pred[, -1], 2, m_cif, time.cif = m_pred[, 1], tau = tau)
    rs <- as.numeric(as.matrix(newdata[, covs, drop = FALSE]) %*% 
                       fg_model$coef)
    return(list(times = newdata$time, status = newdata$status,
                cif_pred = m_pred, pred = pred, linear.pred = rs))
  }
  
  if (!is.null(cr_model)) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      stop("The 'cr_model' argument requires the 'randomForestSRC' package. ",
           "Install it with install.packages(\"randomForestSRC\").",
           call. = FALSE)
    }
    cr_pred <- predict(cr_model, newdata = newdata)
    cif_event <- cr_pred$cif[,,event.type]
    m_pred <- cbind(cr_pred$time.interest, t(cif_event))
    pred <- apply(m_pred[, -1], 2, m_cif, time.cif = m_pred[, 1], tau = tau)
    return(list(times = newdata$time, status = newdata$status,
                cif_pred = m_pred, pred = pred,
                linear.pred = cr_pred$predicted[, event.type]))
  }
  
  if (inherits(model1, "coxph")) {
    model.info1 <- get_CIF(model1, list(model2), newdata)
    m_pred <- cbind(model.info1$y, model.info1$CIF)
    pred <- apply(m_pred[, -1], 2, m_cif, time.cif = m_pred[, 1], tau = tau)
    return(list(times = newdata$time, status = newdata$status,
                cif_pred = m_pred, pred = pred,
                linear.pred = predict(model1, newdata, type = "lp")))
  }
  
  if (inherits(model1, "survreg")) {
    model.info1 <- get_CIF_aft(newdata$time, newdata$status, model1, list(model2),
                               newX = newdata[, covs, drop = FALSE])
    m_pred <- cbind(model.info1$y, model.info1$CIF)
    
    pred <- apply(m_pred[, -1], 2, m_cif, time.cif = m_pred[, 1], tau = tau)
    return(list(times = newdata$time, status = newdata$status,
                cif_pred = m_pred, pred = pred,
                linear.pred = predict(model1, newdata, type = "lp")))
  }
  
  stop("Unknown model type")
}

handle_error <- function(e) {
  warning(conditionMessage(e), call. = FALSE)
  return()
}


#######################################################################
# Internal helper function: calculate CIF based on Cox cause-specific hazard models
#######################################################################


get_CIF <- function(cox.focus, cox.other, X){
  
  fit <- tryCatch({survfit(cox.focus, newdata = X, se.fit = FALSE)}, 
                  error = handle_error)
  Ht <- fit$cumhaz[fit$n.event>=1,]
  if(is.null(dim(Ht))) return()
  hi <- apply(Ht, 2, function(x) diff(c(0, x)))
  
  for(cox.model in cox.other){
    fit_other <- tryCatch({survfit(cox.model, newdata = X, se.fit = FALSE)}, 
                          error = handle_error)
    Ht_other <- fit_other$cumhaz[fit$n.event>=1,]
    Ht <- Ht + Ht_other
    if(is.null(dim(Ht))) return()
  }
  
  CIF <- apply(exp(-Ht) * hi, 2, cumsum)
  return(list(y = fit$time[fit$n.event>=1], CIF = CIF))
}


#######################################################################
# Internal helper function: calculate CIF based on AFT models
#######################################################################


get_CIF_aft <- function(ftime, fstatus, aft.focus, aft.other, newX = NULL) {
  grid <- sort(unique(ftime))
  if (any(!is.finite(grid)) || any(grid < 0))
    stop("AFT CIF prediction requires finite, non-negative follow-up times.")
  fits <- c(list(aft.focus), aft.other)
  params <- lapply(fits, function(fit) {
    if (!inherits(fit, "survreg") ||
        !fit$dist %in% c("weibull", "exponential", "lognormal", "loglogistic"))
      stop("AFT CIF prediction requires Weibull, exponential, lognormal, or loglogistic survreg models.")
    if (length(fit$scale) != 1L || !is.finite(fit$scale) || fit$scale <= 0)
      stop("AFT CIF prediction requires a single positive scale per model.")
    lp <- if (is.null(newX)) stats::predict(fit, type = "lp") else
      stats::predict(fit, newdata = as.data.frame(newX), type = "lp")
    if (any(!is.finite(lp))) stop("AFT CIF linear predictors must be finite.")
    list(dist = fit$dist, scale = fit$scale, lp = as.numeric(lp))
  })
  n <- length(params[[1]]$lp)
  if (any(vapply(params, function(p) length(p$lp), integer(1)) != n))
    stop("AFT CIF models must predict the same subjects.")

  # F1(t) transforms f1(t) * product(S_other(t)) dt to a bounded
  # integrand on [0, F1(t)]. In particular, no zero-density assumption
  # is imposed at t = 0 when the Weibull shape is below 1.
  cdf <- function(z, dist) switch(dist,
    weibull = -expm1(-exp(z)), exponential = -expm1(-exp(z)),
    lognormal = stats::pnorm(z), loglogistic = stats::plogis(z))
  quantile_log <- function(u, p, i) {
    z <- switch(p$dist,
      weibull = log(-log1p(-u)), exponential = log(-log1p(-u)),
      lognormal = stats::qnorm(u), loglogistic = stats::qlogis(u))
    p$lp[i] + p$scale * z
  }
  log_survival <- function(log_t, p, i) {
    z <- (log_t - p$lp[i]) / p$scale
    switch(p$dist,
      weibull = -exp(z), exponential = -exp(z),
      lognormal = stats::pnorm(z, lower.tail = FALSE, log.p = TRUE),
      loglogistic = stats::plogis(z, lower.tail = FALSE, log.p = TRUE))
  }
  target <- params[[1]]
  cif <- matrix(0, nrow = length(grid), ncol = n)
  for (i in seq_len(n)) {
    u <- cdf((log(grid) - target$lp[i]) / target$scale, target$dist)
    integrand <- function(prob) {
      log_t <- quantile_log(prob, target, i)
      log_s <- rep(0, length(prob))
      for (p in params[-1]) log_s <- log_s + log_survival(log_t, p, i)
      exp(log_s)
    }
    # Integrate successive probability intervals to preserve monotonicity
    # and avoid repeatedly integrating the entire interval from zero.
    lower <- 0
    total <- 0
    for (j in seq_along(u)) {
      if (u[j] > lower) {
        total <- total + stats::integrate(integrand, lower, u[j],
          rel.tol = 1e-8, abs.tol = 1e-10)$value
      }
      cif[j, i] <- total
      lower <- u[j]
    }
  }
  list(y = grid, CIF = cif)
}
