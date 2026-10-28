# -------------------------------------------------------------------
# PUBLIC: Case-cohort weights (handles stratified if `strata` given)
# -------------------------------------------------------------------
#' Case-cohort sampling weights
#'
#' Computes Prentice-style weights for case-cohort designs. If `strata`
#' is provided, stratified case-cohort weights are computed.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric \{0,1\}, 1 = event, 0 = censored.
#' @param subcohort Logical vector; TRUE if in subcohort.
#' @param strata Optional factor/character for strata (stratified cc).
#' @return Numeric vector of weights.
#' @export
#' @examples
#' # Compute weights using the full cohort, then retain the random subcohort
#' # and all cases for analysis. Sampled non-cases represent cohort non-cases.
#' d <- pa_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' set.seed(2025)
#' w <- pa_case_cohort_weights(time = d$time, status = d$status,
#'                             subcohort = stats::rbinom(nrow(d), 1, 0.4))
#' summary(w)
pa_case_cohort_weights <- function(time, status, subcohort = NULL, strata = NULL) {
  design <- if (is.null(strata)) "casecohort" else "strat_casecohort"
  weighted_param(time = time, status = status, design = design,
                 subcohort = subcohort, strata = strata, m = NULL)
}
