# -------------------------------------------------------------------
# PUBLIC: NCC weights (handles matched if `strata` given)
# -------------------------------------------------------------------
#' Nested case-control (NCC) sampling weights
#'
#' Computes inverse-selection weights for NCC designs. If `strata` is
#' provided, matched NCC weights are computed.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric \{0,1\}, 1 = event, 0 = censored.
#' @param strata Optional factor/character defining matching strata in the full
#'   phase-one cohort, not IDs of sampled case-control sets.
#' @param m Positive integer, number of controls per case (required).
#' @details
#' Supply full-cohort follow-up and event information, not only sampled rows.
#' Controls are sampled without replacement within each case's risk set and
#' independently across cases. At tied failure times, all simultaneous cases
#' are excluded from the eligible controls. Events precede censoring at ties.
#' Matching restricts both cases and controls to the same stratum. If fewer
#' than \code{m} controls are eligible, all eligible controls are sampled.
#' A non-case with no sampling opportunity has infinite inverse probability;
#' such a subject cannot belong to the sampled analysis dataset. Subset the
#' returned weights to the subjects actually sampled and align their order
#' with the prediction objects before evaluation.
#' @return Numeric vector of weights, in full-cohort input order.
#' @export
#' @examples
#' # Unmatched nested case-control design with m = 2 controls per case.
#' d <- pa_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' set.seed(2025)
#' w <- pa_nested_case_control_weights(time = d$time, status = d$status, m = 2)
#' summary(w)
pa_nested_case_control_weights <- function(time, status, strata = NULL, m = NULL) {
  if (is.null(m)) stop("`m` must be provided for NCC weights.")
  design <- if (is.null(strata)) "ncc" else "matched_ncc"
  weighted_param(time = time, status = status, design = design,
                 subcohort = NULL, strata = strata, m = m)
}

