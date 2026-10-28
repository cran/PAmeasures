#' Calculates sample weights for different types of survival study designs
#'
#' Calculates sampling weights for different types of survival study designs, including full cohort, case-cohort, and nested case-control (NCC) studies.
#' These weights can be used with standard survival analysis functions (like coxph from the survival package) to obtain valid parameter estimates for a Cox proportional hazards model in a subsampled cohort.
#'
#' @param time A numeric vector representing the observed follow-up time for each individual.
#' @param status A numeric vector representing the event status, where 1 indicates an event (e.g., death, disease onset) and 0 indicates censoring.
#' @param design A character string specifying the study design. It can be one of the following:
#'   - `full`: Full cohort analysis (no subsampling). Weights will all be 1.
#'   - `casecohort`:  Unstratified case-cohort design.
#'   - `strat_casecohort`: Stratified case-cohort design.
#'   - `ncc`: Unmatched nested case-control design.
#'   - `matched_ncc`: Matched nested case-control design.
#' @param subcohort A  logical vector (TRUE/FALSE) indicating whether an individual is part of the subcohort. This is required for case-cohort designs.
#' @param strata A factor or character vector specifying the stratum or matching set for each individual. This is required for stratified case-cohort and matched NCC designs.
#' @param  m An integer specifying the number of controls selected for each case. This is required for both NCC designs.
#'
#' @return A numeric vector of sampling weights, one for each individual in the cohort.
#'
# -------------------------------------------------------------------
# INTERNAL: keep but DO NOT export
# -------------------------------------------------------------------
#' @keywords internal
#' @noRd
weighted_param <- function(time, status,
                           design = c("full",
                                      "casecohort",
                                      "strat_casecohort", # Stratified case-cohort
                                      "ncc", # Unmatched nested case-control
                                      "matched_ncc"),
                           subcohort      = NULL,   # TRUE/FALSE for sampled controls
                           strata         = NULL,   # factor for stratum / matching set
                           m              = NULL   # number of controls per case (NCC)
) {
  design <- match.arg(design)
  n      <- length(time)
  stopifnot(length(status) == n,
            all(status %in% c(0, 1)))      # 1 = event, 0 = censored

  ## -----------------------------------------------------------------------
  ##  build the weights  w_k
  ## -----------------------------------------------------------------------
  w <- rep(1, n)                           # cases already = 1

  if (design == "casecohort") {
    if (is.null(subcohort))
      stop("`subcohort` logical vector is required for (unstratified) case-cohort")

    o0 <- sum(status == 0)                 # non-failures in full cohort
    g0 <- sum(subcohort & status == 0)     # sampled non-failures

    w[status == 0] <- o0 / g0              #   w_k = o^0 / g^0

  } else if (design == "strat_casecohort") {
    if (is.null(subcohort) || is.null(strata))
      stop("`strata` *and* `subcohort` must be supplied for stratified case-cohort")

    #w <- rep(NA_real_, n)
    for (s in levels(as.factor(strata))) {
      idx  <- strata == s & status == 0     # non-failures in stratum s
      o0s  <- sum(idx)
      g0s  <- sum(idx & subcohort)
      w[idx] <- o0s / g0s                   #   w_k = o_s^0 / g_s^0
    }
    #w[status == 1] <- 1

  } else if (design %in% c("ncc", "matched_ncc")) {

    if (is.null(m))
      stop("Number of controls per case `m` must be specified for NCC designs")
    if (length(m) != 1L || !is.numeric(m) || !is.finite(m) || m < 1 || m != floor(m))
      stop("`m` must be a positive integer number of controls per case.")

    if (design == "matched_ncc" &&
        (is.null(strata) || length(strata) != n || anyNA(strata)))
      stop("`strata` must have one non-missing matching stratum per cohort member.")

    evt_times <- sort(unique(time[status == 1]))              # all event times
    ## riskset_sizes(t_j)  --------------------------------------------------
    riskset_sizes <- sapply(evt_times, function(tj) sum(time >= tj))

    ## probability a *control* is *never* sampled up to its own t_k --
    temp_prob <- function(k) {
      if (status[k] == 1) return(0)            # cases: probability = 0
      prod_term <- 1
      for (j in seq_along(evt_times)) {
        tj <- evt_times[j]
        if (time[k] >= tj) {                   # events precede censoring at ties
          if (design == "matched_ncc") {
            if (!any(time == tj & status == 1 & strata == strata[k])) next
            Rj <- sum(time >= tj & strata == strata[k])
            failures <- sum(time == tj & status == 1 & strata == strata[k])
          } else {
            Rj <- riskset_sizes[j]
            failures <- sum(time == tj & status == 1)
          }
          # m is the number of CONTROLS, not the matched-set size.
          # For tied failures, each case draws independently from eligible
          # non-cases, without replacement within its control sample.
          eligible <- Rj - failures
          selection_prob <- min(1, m / eligible)
          prod_term <- prod_term * (1 - selection_prob)^failures
        }
      }
      prod_term                                # product over event times
    }

    prob_ever <- sapply(seq_len(n), function(k) 1 - temp_prob(k))
    w[status == 0] <- 1 / prob_ever[status == 0]

  }
  return(w)
}
