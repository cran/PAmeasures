#' Prediction Accuracy for Continuous Outcomes
#'
#' Evaluate observed and predicted continuous responses using the R-squared
#' and L-squared measures of Li and Wang (2019). Predictions may come from a
#' linear or nonlinear model; the original fitted model is not required.
#'
#' @param y A finite numeric vector of observed continuous responses.
#' @param y.predict A finite numeric vector of predictions, paired with
#'   \code{y} in the same order and of the same length. At least two pairs
#'   are required. Missing and infinite values are rejected.
#'
#' @details An ordinary least-squares regression of \code{y} on
#'   \code{y.predict}, with an intercept, gives the corrected predictions
#'   \eqn{\hat y_c}. The measures are
#'   \deqn{R^2 = \frac{\sum_i (\hat y_{c,i} - \bar y)^2}
#'                       {\sum_i (y_i - \bar y)^2}}
#'   and
#'   \deqn{L^2 = \frac{\sum_i (y_i - \hat y_{c,i})^2}
#'                       {\sum_i (y_i - y.predict_i)^2}.}
#'   R-squared describes explained response variation after affine correction;
#'   L-squared compares the corrected and original prediction errors. These
#'   are the original PAmeasures definitions, not \eqn{1 - L^2}.
#'
#'   Constant predictions use the intercept-only correction \code{mean(y)}.
#'   A constant observed response makes R-squared undefined. Predictions
#'   exactly equal to the observations make L-squared undefined because its
#'   denominator is zero. Undefined measures are \code{NA_real_}, with a
#'   warning; any other defined measure is retained. Inputs whose squared
#'   deviations overflow or whose nonzero deviations underflow to zero are
#'   rejected with an error and should be rescaled to common units.
#'
#'   \code{pam.nlm()} is the retained PAmeasures 0.1.0 compatibility entry
#'   point. It uses the same calculation, then rounds and formats each result
#'   to four decimal places as in the historical package. Its ordinary,
#'   defined results retain the historical field names and character format;
#'   input validation and undefined-result handling follow \code{pa_nlm()}.
#'   No deprecation warning is issued for this supported compatibility entry.
#'
#' @return A list with \code{R.squared} and \code{L.squared}.
#'   \code{pa_nlm()} returns unrounded numeric scalars.
#'   \code{pam.nlm()} returns the historical formatted character scalars.
#'
#' @references Li, G. and Wang, X. (2019). Prediction Accuracy Measures for a
#'   Nonlinear Model and for Right-Censored Time-to-Event Data.
#'   \emph{Journal of the American Statistical Association}, 114, 1815--1825.
#'   \doi{10.1080/01621459.2018.1515079}.
#'
#' @examples
#' y <- c(1, 2, 4, 5)
#' prediction <- c(1, 2, 3, 4)
#' pa_nlm(y, prediction)
#' pam.nlm(y, prediction)
#'
#' data(moore)
#' fit <- stats::lm(log2(count) ~ time, data = moore)
#' pa_nlm(moore$count, 2^stats::predict(fit))
#' @export
pa_nlm <- function(y, y.predict) {
  valid <- function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x))
  if (!valid(y) || !valid(y.predict))
    stop("`y` and `y.predict` must be real numeric vectors.", call. = FALSE)
  if (length(y) != length(y.predict) || length(y) < 2L)
    stop("`y` and `y.predict` must have equal lengths of at least two.", call. = FALSE)
  if (any(!is.finite(y)) || any(!is.finite(y.predict)))
    stop("`y` and `y.predict` must contain only finite values.", call. = FALSE)
  # Integer inputs must use double arithmetic before paired subtraction.
  y <- as.numeric(y)
  y.predict <- as.numeric(y.predict)

  total <- sum((y - mean(y))^2)
  original_error <- sum((y - y.predict)^2)
  if (any(!is.finite(c(total, original_error))))
    stop("Squared deviations overflow; rescale both inputs to common units.",
         call. = FALSE)
  constant_response <- all(y == y[1L])
  if ((total == 0 && !constant_response) ||
      (original_error == 0 && any(y != y.predict)))
    stop("Squared deviations underflow; rescale both inputs to common units.",
         call. = FALSE)
  if (constant_response || all(y.predict == y.predict[1L])) {
    corrected <- rep(mean(y), length(y))
  } else {
    corrected <- stats::fitted(stats::lm(y ~ y.predict))
  }
  explained <- sum((corrected - mean(y))^2)
  corrected_error <- sum((y - corrected)^2)
  if (any(!is.finite(c(total, explained, corrected_error, original_error))))
    stop("Squared deviations overflow; rescale both inputs to common units.",
         call. = FALSE)

  undefined <- character()
  r_squared <- if (total == 0) {
    undefined <- c(undefined, "R.squared (constant observed response)")
    NA_real_
  } else explained / total
  l_squared <- if (original_error == 0) {
    undefined <- c(undefined, "L.squared (zero original prediction error)")
    NA_real_
  } else corrected_error / original_error
  if (length(undefined))
    warning("Undefined prediction accuracy measure: ",
            paste(undefined, collapse = "; "), ".", call. = FALSE)
  list(R.squared = r_squared, L.squared = l_squared)
}

#' @rdname pa_nlm
#' @export
pam.nlm <- function(y, y.predict) {
  lapply(pa_nlm(y, y.predict), function(value)
    format(round(value, digits = 4), nsmall = 4))
}
