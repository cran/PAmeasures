#' @title Prediction Accuracy Measures for Cox Proportional Hazards Model
#'
#' @description This function calculates a set of measures, $R^2$ and $L^2$, $Pseudo-R^2$, for a Cox proportional hazards model. When predict = TRUE, 
#' it returns the predicted survival times alongside observed values and event status. 
#' R-squared quantifies the proportion of variability in the observed survival times explained by the model's predictions,
#' while L-squared measures the proportion of prediction error explained by a corrected prediction function. Together, these metrics provide a comprehensive evaluation of the predictive power of the Cox model.
#'
#' @param model An object of class `coxph` representing a fitted Cox proportional hazards regression model. The model must be fitted with `x = TRUE` and `y = TRUE` to include the design matrix and response vector.
#' @param covs A character vector specifying the names of the covariates used in the model.
#' @param tau (Optional) A restriction time for computing restricted survival time. If provided, the function evaluates predictions up to this time point. Default is 10e10.
#' @param new_data Optional evaluation data. Its outcomes are returned for metric
#'   computation and define the evaluation time grid; they never estimate the
#'   baseline hazard. All prediction parameters come from the training fit.
#' @param predict A logical value indicating whether to return individual predictions. Default is TRUE.
#' @details Probabilities use the fitted training baseline, including its case
#'   weights, strata, tie method and offsets. Offset centering follows the
#'   weighted training mean used by [survival::survfit()].
#'   With `predict = FALSE`, R-squared metrics return missing values with a
#'   warning if censoring weights cannot be normalized, for example when the
#'   training censoring survival is zero at an observed evaluation event.
#' @return
#' A **named list** of model-specific prediction objects.  
#' Each element of the returned list corresponds to one fitted model and contains:
#' \itemize{
#'   \item \code{model} -- the fitted survival model object.
#'   \item \code{times} -- numeric vector of observed follow-up times for the evaluation dataset.
#'   \item \code{status} -- event indicators (1 = event, 0 = censored).
#'   \item \code{surv_prob} -- an \eqn{n \times K} matrix of predicted
#'         subject-specific survival probabilities on a common time grid.
#'   \item \code{pred} -- numerical integral of the training-model survival curve
#'         on the evaluation grid, up to its maximum time or tau, whichever is
#'         smaller. This preserves the legacy integration rule; it is not an
#'         identified unrestricted mean. Fix the grid and tau for comparisons.
#'   \item \code{covs} -- character vector of covariate names used for prediction.
#'   \item \code{new_data} -- dataset on which the predictions were computed.
#' }
#' 
#' @references
#' Li, G., & Wang, X. (2016). Prediction Accuracy Measures for a Nonlinear Model and for Right-Censored Time-to-Event Data.
#' arXiv preprint arXiv:1611.03063. Available at \url{https://arxiv.org/abs/1611.03063}.
#'
#' @examples
#' data(pbc, package = "survival")
#' pbc <- pbc[!is.na(pbc$trt), ]
#' pbc$log_albumin <- log(pbc$albumin)
#' pbc$log_bili    <- log(pbc$bili)
#' pbc$log_protime <- log(pbc$protime)
#' pbc$status <- ifelse(pbc$status == 2, 1, 0)
#' covs <- c("age", "log_albumin", "log_bili", "log_protime", "edema")
#' df <- pbc[, c("time", "status", covs)]
#' df <- df[stats::complete.cases(df), ]
#'
#' set.seed(2025)
#' idx <- sample(seq_len(nrow(df)), round(2 / 3 * nrow(df)))
#' train_data <- df[idx, ]
#' test_data  <- df[-idx, ]
#'
#' m.cox <- survival::coxph(survival::Surv(time, status) ~ .,
#'                          data = train_data, x = TRUE, y = TRUE)
#' cox_pred <- pa_predict_coxph(model = m.cox, covs = covs,
#'                                  new_data = test_data)
#' str(cox_pred$pred)
#'
#' @export
pa_predict_coxph <- function(model, covs, tau = 10e10, new_data = NULL, predict = TRUE)
{
  if (!inherits(model, "coxph") || inherits(model, "coxphms")) {
    stop("model must be a right-censored coxph fit.")
  }
  if (is.null(model$y) || ncol(model$y) != 2L ||
      attr(model$y, "type") != "right") {
    stop("model must include y = TRUE and a right-censored response.")
  }
  if (!is.null(tau) && (length(tau) != 1L || !is.finite(tau) || tau <= 0)) {
    stop("tau must be NULL or a positive finite number.")
  }
  y.unsorted <- model$y[, 1]
  censor.unsorted <- model$y[, 2]
  if (is.null(new_data)) {
    prediction_data <- stats::model.frame(model)
    y.order.new <- NULL
  } else {
    if (!all(c("time", "status") %in% names(new_data))) {
      stop("new_data require time and status columns")
    }
    if (!all(covs %in% names(new_data))) {
      stop("All covariates must be present in new_data.")
    }
    if (!nrow(new_data) || any(!is.finite(new_data$time)) ||
        any(new_data$time < 0) || anyNA(new_data$status) ||
        any(!new_data$status %in% c(0, 1))) {
      stop("new_data require finite non-negative time and binary status.")
    }
    prediction_data <- new_data
    y.unsorted.new <- new_data$time
    censor.unsorted.new <- new_data$status
    y.order.new <- order(y.unsorted.new)
  }
  y.order <- order(y.unsorted)
  y <- y.unsorted[y.order]
  delta <- censor.unsorted[y.order]
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  if (!is.null(y.order.new)){
    y.order <- y.order.new
    y <- y.unsorted.new[y.order]
    delta <- censor.unsorted.new[y.order]
  } 
  km.censoring.minus <- pa_censoring_left_limit(fit.km.censoring, y)
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  prediction_order <- if (is.null(y.order.new)) order(y.unsorted) else y.order.new
  prediction_data <- prediction_data[prediction_order, , drop = FALSE]
  temp_rs <- if (is.null(new_data)) {
    as.numeric(stats::predict(model, type = "lp", reference = "sample"))[prediction_order]
  } else {
    as.numeric(stats::predict(model, newdata = prediction_data,
                              type = "lp", reference = "sample"))
  }
  # The baseline, including training weights and the fitted tie convention,
  # comes exclusively from the training fit. Validation outcomes only supply
  # the evaluation grid and the responses passed to downstream metrics.
  S.hat.x <- t(pa_cox_survival(model, prediction_data, temp_rs, y))

  # Truncate data if tau is provided
  delta.origin <- delta
  y.origin <- y
  if (!is.null(tau)) {
    delta <- ifelse(y <= tau, delta, 1)
    y <- pmin(y, tau)
    y.input <- ifelse(tau > y, y, tau)
  } else {
    y.input <- NULL
  }
  
  
  
  if(is.null(y.input)){ 
    t1 <- y
  }else{
    t1 <- y.input #[y.order]
  }
  
  t2 <- c(0, t1[1:length(t1) - 1])
  delta.t <- t1 - t2
  t.predicted <- colSums(delta.t * S.hat.x)
  t.predicted2 <- integrate_survival(t(S.hat.x), t1, delta, tau) # restricted mean survival time
  if(predict){return(list(pred = t.predicted2, #t.predicted, 
                          times = y.origin,
                          status = delta.origin,
                          surv_prob = t(S.hat.x),
                          linear.pred = temp_rs,
                          model = model,
                          new_data = new_data,
                          covs = covs
                                     ))}

  if (any(!is.finite(weight.km)) || !any(weight.km > 0)) {
    warning("R-squared metrics require positive censoring survival at observed ",
            "evaluation events and at least one observed event.", call. = FALSE)
    return(list(R.squared = NA_character_, L.squared = NA_character_,
                Psuedo.R = NA_character_))
  }
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * 
                                                        y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- round(num.rho2/denom.rho2, digits = 4)
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - t.predicted)^2)
  L2 <- round(num.L2/denom.L2, digits = 4)
  SR <- round(R2 * L2, digits = 4)
  return(list(R.squared = format(R2, nsmall = 4), 
              L.squared = format(L2, nsmall = 4), 
              Psuedo.R = format(SR, nsmall = 4)))
}


# Training-baseline survival on an arbitrary sorted evaluation grid. Chunking
# avoids allocating a training-times x validation-subjects survfit object.
pa_cox_survival <- function(model, new_data, linear_predictor, times) {
  # The predictions below are subject-specific, so survival's warning about the
  # default curve of an interaction model does not apply.
  bh <- tryCatch(
    withCallingHandlers(
      survival::basehaz(model, centered = TRUE),
      warning = function(w) {
        if (grepl("contains interactions", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }),
    error = function(e) {
      if (grepl("strata by covariate interaction", conditionMessage(e))) NULL
      else stop(e)
    })
  if (is.null(bh)) return(pa_cox_survival_by_row(model, new_data, times))
  if (!nrow(bh)) stop("The training fit has no baseline time grid.")
  if (length(attr(stats::terms(model), "offset"))) {
    # survfit/basehaz centers offsets at their weighted training mean. predict's
    # lp centering differs between in-sample and newdata calls, so construct the
    # relative risk on the baseline's scale from the fitted design and offsets.
    training_frame <- stats::model.frame(model)
    training_offset <- stats::model.offset(training_frame)
    training_weight <- stats::model.weights(training_frame)
    offset_center <- if (is.null(training_weight)) mean(training_offset) else
      stats::weighted.mean(training_offset, training_weight)
    prediction_frame <- if (!is.null(attr(new_data, "terms"))) new_data else
      stats::model.frame(stats::delete.response(stats::terms(model)), new_data,
                        xlev = model$xlevels, na.action = stats::na.pass)
    design <- stats::model.matrix(model, data = prediction_frame)
    beta <- stats::coef(model)
    beta[is.na(beta)] <- 0
    linear_predictor <- as.numeric(sweep(design, 2, model$means, "-") %*% beta) +
      stats::model.offset(prediction_frame) - offset_center
  }
  if (is.null(bh$strata)) {
    bh$strata <- "all"
    subject_strata <- rep("all", nrow(new_data))
  } else {
    tt <- stats::delete.response(stats::terms(model))
    specials <- attr(tt, "specials")$strata
    term_names <- rownames(attr(tt, "factors"))[specials]
    mf <- if (all(term_names %in% names(new_data))) new_data else
      stats::model.frame(tt, new_data, xlev = model$xlevels, na.action = stats::na.pass)
    # delete.response updates the special indices to the response-free frame.
    sf <- mf[, term_names, drop = FALSE]
    subject_strata <- if (ncol(sf) == 1L) as.character(sf[[1L]]) else
      as.character(survival::strata(sf, shortlabel = TRUE))
    if (any(!subject_strata %in% as.character(bh$strata))) {
      stop("New data contain a stratum absent from the training fit.")
    }
  }
  out <- matrix(NA_real_, nrow(new_data), length(times))
  for (s in unique(subject_strata)) {
    baseline <- bh[as.character(bh$strata) == s, , drop = FALSE]
    at <- findInterval(times, baseline$time)
    h0 <- c(0, baseline$hazard)[at + 1L]
    rows <- which(subject_strata == s)
    for (start in seq.int(1L, length(rows), by = 256L)) {
      block <- rows[start:min(start + 255L, length(rows))]
      out[block, ] <- exp(-outer(exp(linear_predictor[block]), h0))
    }
  }
  out
}

# Strata-by-covariate interactions have no shared stratum baseline, so each
# subject's curve comes from survfit() on the raw covariates.
#' @noRd
pa_cox_survival_by_row <- function(model, new_data, times) {
  if (!is.null(attr(new_data, "terms"))) {
    stop("Cox models with strata-by-covariate interactions need raw covariates; ",
         "pass the evaluation data as new_data.", call. = FALSE)
  }
  out <- matrix(NA_real_, nrow(new_data), length(times))
  for (i in seq_len(nrow(new_data))) {
    fit <- survival::survfit(model, newdata = new_data[i, , drop = FALSE],
                             se.fit = FALSE)
    out[i, ] <- c(1, fit$surv)[findInterval(times, fit$time) + 1L]
  }
  out
}
