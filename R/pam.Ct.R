# Reverse Kaplan-Meier estimate of the CENSORING distribution G(t) = P(C > t).
#
# Two conventions matter and both were wrong before:
#
#  1. The censoring "event" is an observation that was censored (status == 0).
#     The previous coding, ifelse(status == min(status), 1, 0), inverted this
#     whenever a dataset contained no censoring at all: min(status) was then 1,
#     so every observation counted as a censoring event and G decayed from 1
#     instead of staying at 1.
#
#  2. At a tied time, censorings are taken to occur AFTER events, so the risk
#     set for the censoring hazard is Y(t) - d(t): the number at risk less the
#     events at that same time. A plain survfit(Surv(time, status == 0)) treats
#     the two symmetrically and yields a systematically larger G. This is the
#     reverse-Kaplan-Meier convention used by prodlim(reverse = TRUE), and
#     therefore by pec::ipcw(), which the tests check against.
#
# Returns a step function as a list of jump times and the value of G at each.
#' @keywords internal
#' @noRd
gt_censoring_fit <- function(object) {
  time   <- object[, 1]
  status <- object[, 2]

  utime <- sort(unique(time))
  n     <- length(time)

  n_risk  <- vapply(utime, function(u) sum(time >= u), numeric(1))
  n_event <- vapply(utime, function(u) sum(time == u & status != 0), numeric(1))
  n_cens  <- vapply(utime, function(u) sum(time == u & status == 0), numeric(1))

  # risk set for the censoring hazard: events at the same time are removed first
  at_risk_for_cens <- n_risk - n_event
  haz <- ifelse(at_risk_for_cens > 0, n_cens / at_risk_for_cens, 0)

  list(time = utime, surv = cumprod(1 - haz), n = n)
}

Gt <- function(object, timepoint, left = FALSE) {
  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }

  if (missing(object)) {
    stop("The survival object is missing")
  }

  if (missing(timepoint)) {
    stop("The time is missing with no default")
  }

  if (any(timepoint <= 0)) {
    stop("The timepoint must be positive")
  }

  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (length(timepoint) != 1) {
    stop("Gt can only be calculated at a single time point")
  }

  if (is.na(timepoint)) {
    stop("Cannot calculate Gt at NA")
  }

  fit <- gt_censoring_fit(object)

  # Beyond the last observed time the censoring distribution is not identified.
  # pec::ipcw() returns NA there; so do we.
  if (timepoint > max(fit$time)) {
    return(NA_real_)
  }

  # Step-function evaluation. G(t) is right-continuous, so it takes the value at
  # the last jump time at or before t; G(t-) takes the value strictly before t.
  # Before the first jump both are 1.
  idx <- if (left) which(fit$time < timepoint) else which(fit$time <= timepoint)
  Gvalue <- if (length(idx) == 0L) 1 else fit$surv[max(idx)]

  if (!is.na(Gvalue) && Gvalue <= 0) {
    stop("The censoring distribution has reached zero at or before time ",
         format(timepoint), ", so the inverse-probability-of-censoring weight ",
         "is undefined. Evaluate at an earlier time.", call. = FALSE)
  }

  Gvalue
}
