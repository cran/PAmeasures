#' Calibration of Predicted Event Probabilities at One Time Point
#'
#' Assess calibration from externally supplied probabilities and observed
#' outcomes, without requiring the original prediction model.
#'
#' @param p Numeric vector of predicted event probabilities at \code{t_star},
#'   in the same subject order as \code{event_time} and \code{status}. For
#'   ordinary survival supply \code{1 - S(t_star)}; for competing risks supply
#'   the CIF of \code{event_type} at \code{t_star}.
#' @param event_time Numeric vector of finite, non-negative observed times.
#' @param status Event indicators: 0 for censoring and 1 for an event in
#'   \code{pa_calibration()}; non-negative integer event codes in
#'   \code{pa_calibration_cr()} (0 for censoring).
#' @param t_star One positive finite evaluation time, not exceeding the maximum
#'   observed follow-up. Use the same time units as \code{event_time}.
#' @param nk Number of restricted cubic spline knots (at least 3; default 5).
#' @param case_weights Optional positive finite inverse-inclusion sampling
#'   weights, aligned with the subjects. For \code{pa_calibration()} only.
#'   Used in both the calibration model and the metric summaries.
#' @param event_type Positive integer code for the cause of interest.
#' @param bandwidth Nearest-neighborhood smoothing bandwidth passed to
#'   \code{prodlim::prodlim()}. \code{NULL} uses its automatic bandwidth.
#'
#' @details
#' \code{pa_calibration()} uses \code{pmcalibration::pmcalibration()} with a
#' restricted cubic spline and no confidence interval computation. The knot
#' count is capped at the number of distinct fitting probabilities; with fewer
#' than five distinct probabilities or fewer than six subjects it uses the
#' backend's single-predictor Cox calibration instead. Constant predictions
#' are compared with the marginal Kaplan--Meier risk.
#'
#' The right-censored backend requires probabilities strictly inside (0, 1).
#' Probabilities within \code{sqrt(.Machine$double.eps)} of 0 or 1 are clipped
#' to that distance from the boundary for fitting, with a warning (except for
#' the constant-prediction Kaplan--Meier branch). Metrics compare calibrated values
#' with the original probabilities and retain every subject.
#'
#' \code{pa_calibration_cr()} uses local Aalen--Johansen estimates from
#' \code{prodlim}, conditional on the predicted CIF. All nonzero causes other
#' than \code{event_type} are pooled as competing events, not censored.
#' Constant or few-valued predictions use the backend's discrete risk groups.
#' This is a different curve estimator from Fine--Gray spline calibration.
#'
#' For unweighted evaluations, both functions use
#' \code{pmcalibration::cal_metrics()} to summarize absolute
#' differences between predicted and curve-implied probabilities at the original
#' subjects. Only \code{ici} (the backend's \code{Eavg}), \code{e50}, and
#' \code{e90} are returned. They are probability-scale errors, with 0 ideal.
#' ICI averages over subjects, not over follow-up time or a uniform risk grid.
#'
#' With \code{case_weights}, \code{pa_calibration()} fits a sampling-weighted
#' Cox model directly, using the same complementary log-log transformation
#' and restricted cubic spline knot-placement rule as the ordinary backend.
#' Both coefficients and baseline cumulative hazard come from the weighted
#' fit. Constant predictions use weighted Kaplan--Meier risk; with fewer than
#' five distinct fitting probabilities or six subjects a single linear
#' predictor is used. Knot locations are recorded in \code{knots}.
#' ICI is the normalized sampling-weighted mean absolute error. E50 and E90
#' use weighted order-statistic interpolation, with weights normalized to the
#' sample size (\code{Hmisc::wtd.quantile(..., normwt = TRUE)}).
#' Multiplying every sampling weight by the same constant has no effect.
#'
#' Supply subject inclusion weights for the actual case--cohort/NCC design,
#' rather than weights already combined with censoring adjustments. Only
#' sampled subjects with positive inclusion probabilities should be supplied.
#' The calibration Cox model assumes conditionally independent censoring given
#' its predictor and an appropriate proportional-hazards spline relationship.
#' The competing-risks function remains an unweighted cohort evaluation.
#' No optimism correction or design-based confidence interval is applied;
#' preferably supply held-out or out-of-fold predictions.
#'
#' @return A list with \code{metrics}, a named numeric vector
#'   \code{c(ici, e50, e90)}; \code{curve}, a data frame with \code{predicted}
#'   and \code{calibrated} probabilities in input subject order;
#'   \code{t_star}; \code{event_type}; and \code{method}, the curve estimator.
#'   Weighted results also include the original \code{case_weights} and spline
#'   \code{knots} on the complementary log-log scale (empty for simpler fits).
#'
#' @references
#' Austin PC, Steyerberg EW (2019). The Integrated Calibration Index (ICI) and
#' related metrics for quantifying the calibration of logistic regression
#' models. Statistics in Medicine 38, 4051--4065. \doi{10.1002/sim.8281}.
#' Austin PC, Harrell FE Jr, van Klaveren D (2020). Graphical calibration curves
#' and the integrated calibration index (ICI) for survival models. Statistics
#' in Medicine 39, 2714--2742. \doi{10.1002/sim.8570}.
#' Breslow NE, Wellner JA (2007). Weighted likelihood for semiparametric models
#' and two-phase stratified samples, with application to Cox regression.
#' Scandinavian Journal of Statistics 34, 86--102.
#' \doi{10.1111/j.1467-9469.2006.00523.x}.
#' Samuelsen SO (1997). A pseudolikelihood approach to analysis of nested
#' case-control studies. Biometrika 84, 379--394.
#' \doi{10.1093/biomet/84.2.379}.
#'
#' @examples
#' set.seed(19)
#' x <- stats::rnorm(100)
#' event <- stats::rexp(100, rate = exp(x))
#' censor <- stats::rexp(100, rate = 0.3)
#' time <- pmin(event, censor)
#' status <- as.integer(event <= censor)
#' p <- 1 - exp(-exp(x) * 0.5)
#' pa_calibration(p, time, status, t_star = 0.5)$metrics
#'
#' # Externally supplied CIFs; no fitted model object is needed.
#' competing <- stats::rexp(100, rate = 0.4)
#' cr_time <- pmin(event, competing, censor)
#' cr_status <- ifelse(censor < pmin(event, competing), 0L,
#'                     ifelse(event <= competing, 1L, 2L))
#' cif <- exp(x) / (exp(x) + 0.4) * (1 - exp(-(exp(x) + 0.4) * 0.5))
#' pa_calibration_cr(cif, cr_time, cr_status, t_star = 0.5)$metrics
#' @export
pa_calibration <- function(p, event_time, status, t_star, nk = 5,
                           case_weights = NULL) {
  pa_validate_calibration(p, event_time, status, t_star, competing = FALSE)
  # A named quantile (e.g. the default "50%") must not become row names in
  # the backend's survival prediction frames.
  t_star <- unname(t_star)
  if (!is.numeric(nk) || length(nk) != 1L || !is.finite(nk) ||
      nk < 3 || nk != as.integer(nk)) {
    stop("'nk' must be an integer of at least 3.", call. = FALSE)
  }
  if (!is.null(case_weights)) {
    return(pa_calibration_weighted(p, event_time, status, t_star, nk,
                                    case_weights))
  }
  if (length(unique(p)) == 1L) {
    km <- survival::survfit(survival::Surv(event_time, status) ~ 1)
    calibrated <- rep(1 - summary(km, times = t_star, extend = TRUE)$surv,
                      length(p))
    method <- "kaplan_meier_constant"
  } else {
    fitting_p <- pmin(1 - sqrt(.Machine$double.eps),
                      pmax(sqrt(.Machine$double.eps), p))
    if (any(fitting_p != p)) {
      warning("Boundary probabilities were moved inside (0, 1) for calibration ",
              "fitting only; all subjects and original probabilities are ",
              "retained in the metrics.", call. = FALSE)
    }
    distinct <- length(unique(fitting_p))
    # Hmisc needs at least three interior knot locations; 3/4-valued risks
    # cannot provide these. Match the weighted path's linear Cox fallback.
    smooth <- if (distinct < 5L || length(p) < 6L) "none" else "rcs"
    cal <- pmcalibration::pmcalibration(
      y = survival::Surv(event_time, status), p = fitting_p, time = t_star,
      smooth = smooth, nk = min(nk, distinct), ci = "none", plot = FALSE)
    calibrated <- cal$p_c$p_c
    method <- paste0("pmcalibration_", smooth)
  }
  pa_calibration_result(p, calibrated, t_star, 1L, method)
}

# The weighted path is implemented here rather than passing unsupported weights
# through pmcalibration's smoothing arguments. Hmisc supplies the same spline
# basis as that backend; survival handles the weighted Cox estimating equations.
pa_calibration_weighted <- function(p, event_time, status, t_star, nk,
                                     case_weights) {
  pa_validate_calibration_weights(case_weights, length(p))
  # Relative weights keep the Cox fit invariant to weight units and avoid
  # numerical overflow when inclusion probabilities are small.
  weight <- case_weights / max(case_weights)
  weight <- weight / mean(weight)
  knots <- numeric()
  if (length(unique(p)) == 1L) {
    fit <- survival::survfit(survival::Surv(event_time, status) ~ 1,
                            weights = weight)
    calibrated <- rep(1 - summary(fit, times = t_star, extend = TRUE)$surv,
                      length(p))
    method <- "weighted_kaplan_meier_constant"
  } else {
    fitting_p <- pmin(1 - sqrt(.Machine$double.eps),
                      pmax(sqrt(.Machine$double.eps), p))
    if (any(fitting_p != p)) {
      warning("Boundary probabilities were moved inside (0, 1) for calibration ",
              "fitting only; all subjects and original probabilities are ",
              "retained in the metrics.", call. = FALSE)
    }
    z <- log(-log1p(-fitting_p))
    distinct <- length(unique(z))
    if (distinct < 5L || length(p) < 6L) {
      X <- matrix(z, ncol = 1L)
      method <- "weighted_cox_linear"
    } else {
      X <- Hmisc::rcspline.eval(z, nk = min(nk, distinct), inclx = TRUE)
      knots <- attr(X, "knots")
      method <- "weighted_cox_rcs"
    }
    colnames(X) <- paste0("spline", seq_len(ncol(X)))
    d <- data.frame(time = event_time, event = status, X, weight = weight)
    predictors <- colnames(X)
    formula <- stats::reformulate(predictors,
                                  response = "survival::Surv(time, event)")
    fit <- survival::coxph(formula, data = d, weights = weight,
                           x = TRUE, model = TRUE, robust = TRUE)
    # Expected cumulative hazard must use the weighted fit's baseline too.
    newdata <- data.frame(time = t_star, event = 1, X)
    calibrated <- -expm1(-as.numeric(stats::predict(
      fit, newdata = newdata, type = "expected")))
  }
  out <- pa_calibration_result(p, calibrated, t_star, 1L, method, case_weights)
  out$case_weights <- case_weights
  out$knots <- knots
  out
}

pa_validate_calibration_weights <- function(case_weights, n) {
  if (!is.numeric(case_weights) || !is.null(dim(case_weights)) ||
      length(case_weights) != n || any(!is.finite(case_weights)) ||
      any(case_weights <= 0)) {
    stop("'case_weights' must be a positive finite sampling-weight vector, ",
         "one per subject in prediction order.", call. = FALSE)
  }
}

#' @rdname pa_calibration
#' @export
pa_calibration_cr <- function(p, event_time, status, t_star, event_type = 1,
                              bandwidth = NULL) {
  pa_validate_calibration(p, event_time, status, t_star, competing = TRUE,
                         event_type = event_type)
  t_star <- unname(t_star)
  if (!is.null(bandwidth) &&
      !(is.character(bandwidth) && identical(bandwidth, "smooth")) &&
      (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
       !is.finite(bandwidth) || bandwidth <= 0)) {
    stop("'bandwidth' must be NULL, 'smooth', or a positive finite number.",
         call. = FALSE)
  }
  # Explicit recoding keeps arbitrary original cause labels from changing
  # the meaning of predict(..., cause = 1), including when a cause is absent.
  code <- ifelse(status == 0, 0L, ifelse(status == event_type, 1L, 2L))
  d <- data.frame(time = event_time, event = code, p = p)
  fit <- prodlim::prodlim(prodlim::Hist(time, event) ~ p, data = d,
                         bandwidth = bandwidth)
  calibrated <- as.numeric(stats::predict(
    fit, newdata = data.frame(p = p), times = t_star, cause = 1,
    type = "cuminc", mode = "matrix"))
  pa_calibration_result(p, calibrated, t_star, event_type, "prodlim_nne")
}

pa_validate_calibration <- function(p, event_time, status, t_star,
                                    competing, event_type = 1) {
  n <- length(event_time)
  if (!is.numeric(event_time) || !is.null(dim(event_time)) || n < 2L ||
      any(!is.finite(event_time)) || any(event_time < 0)) {
    stop("'event_time' must contain at least two finite non-negative times.",
         call. = FALSE)
  }
  if (!is.numeric(p) || !is.null(dim(p)) || length(p) != n ||
      any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("'p' must be a finite probability vector in [0, 1], one per subject.",
         call. = FALSE)
  }
  if (!(is.numeric(status) || is.logical(status)) || !is.null(dim(status)) ||
      length(status) != n || any(!is.finite(status)) ||
      any(status < 0 | status != floor(status)) ||
      (!competing && any(!status %in% c(0, 1)))) {
    stop(if (competing) "'status' must contain non-negative integer event codes."
         else "'status' must contain only 0 (censored) and 1 (event).",
         call. = FALSE)
  }
  if (!is.numeric(event_type) || length(event_type) != 1L ||
      !is.finite(event_type) || event_type <= 0 || event_type != floor(event_type)) {
    stop("'event_type' must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(t_star) || length(t_star) != 1L || !is.finite(t_star) ||
      t_star <= 0 || t_star > max(event_time)) {
    stop("'t_star' must be positive and within observed follow-up.", call. = FALSE)
  }
  if (!any(status == event_type)) {
    stop("Calibration requires at least one observed target event.", call. = FALSE)
  }
}

pa_calibration_result <- function(p, calibrated, t_star, event_type, method,
                                   case_weights = NULL) {
  if (length(calibrated) != length(p) || any(!is.finite(calibrated)) ||
      any(calibrated < 0 | calibrated > 1)) {
    stop("Calibration curve is not estimable for all subjects at 't_star'; ",
         "check follow-up support and smoothing settings.", call. = FALSE)
  }
  if (is.null(case_weights)) {
    values <- pmcalibration::cal_metrics(p = p, p_c = calibrated)
    metrics <- stats::setNames(unname(values[c("Eavg", "E50", "E90")]),
                               c("ici", "e50", "e90"))
  } else {
    error <- abs(p - calibrated)
    weight <- case_weights / max(case_weights)
    quantiles <- Hmisc::wtd.quantile(error, weights = weight,
      probs = c(.5, .9), normwt = TRUE)
    metrics <- c(ici = sum(weight * error) / sum(weight),
                  e50 = unname(quantiles[1]), e90 = unname(quantiles[2]))
  }
  list(metrics = metrics,
       curve = data.frame(predicted = p, calibrated = calibrated),
       t_star = t_star, event_type = event_type, method = method)
}

# Predictions on supplied grids are interpreted as right-continuous steps.
# Never use a future grid value to assess calibration at t_star.
pa_calibration_grid_risk <- function(probability, grid, t_star, n,
                                      survival = FALSE) {
  probability <- as.matrix(probability)
  if (!is.numeric(grid) || !length(grid) || any(!is.finite(grid)) ||
      any(grid < 0) || is.unsorted(grid) ||
      !is.numeric(probability) || nrow(probability) != n ||
      ncol(probability) != length(grid) ||
      any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
    stop("Calibration requires a finite probability matrix (subjects by time), ",
         "and a matching non-decreasing non-negative time grid.", call. = FALSE)
  }
  if (!is.numeric(t_star) || length(t_star) != 1L || !is.finite(t_star) ||
      t_star <= 0 || t_star > max(grid)) {
    stop("Calibration 't_star' must be positive and within the prediction grid.",
         call. = FALSE)
  }
  i <- findInterval(t_star, grid)
  if (i == 0L) return(rep(0, n))
  p <- probability[, i]
  if (survival) 1 - p else p
}

# A non-estimable calibration curve must not suppress other requested metrics.
pa_calibration_values <- function(fun, args) {
  tryCatch(do.call(fun, args)$metrics, error = function(e) {
    warning("Calibration unavailable: ", conditionMessage(e), call. = FALSE)
    stats::setNames(rep(NA_real_, 3L), c("ici", "e50", "e90"))
  })
}
