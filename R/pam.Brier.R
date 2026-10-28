#' @title Brier Score Calculation
#'
#' @description
#' The Brier Score, proposed by Glenn W. Brier in 1950, is a proper score function used to measure the accuracy of probabilistic predictions. It is commonly applied to assess model fits for survival data. The Brier Score can be calculated at any time point, regardless of whether it coincides with the event time.
#'
#' The Brier Score represents the mean squared difference between true classes and predicted probabilities, effectively serving as a cost function. A lower Brier Score indicates better-calibrated predictions. Its values range from zero to one, as it reflects the maximum possible squared difference between predicted probabilities and actual outcomes. 
#'
#' In the context of censored samples, where the exact time of an event (e.g., death) is unknown, direct calculation of residuals is not feasible. Thus, the Brier Score is widely utilized in survival analysis.
#'
#' The Brier Score is considered a strictly proper score (Gneiting and Raftery, 2007), meaning it achieves its minimum value only when the predicted probabilities align with empirical probabilities. Empirical evidence suggests that predictions of survival duration can be inaccurate; however, incorporating patient-specific survival probabilities along with the Brier Score improves the ability to differentiate between future survivors and failures.
#'
#' @importFrom survival Surv
#' @importFrom pec predictSurvProb
#' @importFrom stats median
#' @param object An object of class \code{Surv}, created by the \code{Surv} function, or a fitted survival model, such as those produced by \code{coxph}, \code{survreg}, or \code{rfsrc}.
#' @param pre_sp 
#' If \code{object} is a fitted survival model, this parameter should be a dataset on which you want to calculate the Brier Score. If \code{object} is a survival object, this parameter should be a vector of predicted survival probabilities for each observation at time \code{t_star}.
#' @param t_star A specified time point for calculating the Brier Score. This is necessary when \code{object} is a fitted survival model, as it indicates when the survival probability is predicted. If \code{object} is a survival object, this parameter can be ignored and does not affect the function's outcome.
#'
#' @return The Brier Score at time \code{t_star}, representing the difference between true classes and predicted probabilities.
#'
#' @references
#' Graf, E., Schmoor, C., Sauerbrei, W., & et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. *Statistical Medicine*, 18(17-18), 2529-2545.
#'
#' Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. *Monthly Weather Review*, 78.
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.
#' 
#' Zhou, H., Cheng, X., Wang, S., Zou, Y., & Wang, H. (2022). SurvMetrics: Predictive Evaluation Metrics in Survival Analysis. 
#' R package version 0.5.0. Available at \url{https://github.com/skyee1/SurvMetrics}.
#' @examples
#' library(survival)
#'
#'# Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#'
#'# Fit an exponential model with bilirubin
#' fit.coxph.full <- coxph(Surv(time, status) ~ age + log_albumin + 
#'                          log_bili + log_protime + edema, 
#'                       data = pbc,x=TRUE,y=TRUE)
#' taulist <- seq(0, max(pbc$time), 300)
#' # Choose a specific time point (t_star) for calculating the Brier Score 
#' # For simplicity, we use the median survival time as t_star.
#' t_star <- median(pbc$time)
#' 
#' pam.Brier(fit.coxph.full, pbc, t_star)
#' 
#' @keywords internal
#' @noRd

pam.Brier <- function(object, pre_sp, t_star = -1) {
  # case1, coxph AND testing set
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- pre_sp
    t_star0 <- t_star
    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))
    
    if ( t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point
    
    vec_coxph <-
      predictSurvProb(obj, test_data, t_star0) # get the survival probability vector
    object_coxph <- Surv(test_data$time, test_data$status)
    
    object <- object_coxph
    pre_sp <- vec_coxph
    t_star <- t_star0
  }
  
  
  # case2, RSF AND testing set
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- pre_sp
    t_star0 <- t_star
    
    # the interesting times of training set
    distime <- obj$time.interest
    
    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point
    
    med_index <- order(abs(distime - t_star0))[1]
    
    mat_rsf <-
      predict(obj, test_data)$survival # get the survival probability matrix
    vec_rsf <-
      mat_rsf[, med_index] # get the survival probability vector
    object_rsf <- Surv(test_data$time, test_data$status)
    
    object <- object_rsf
    pre_sp <- vec_rsf
    t_star <- t_star0
  }
  
  # case3 survreg AND testing set
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- pre_sp
    t_star0 <- t_star
    
    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))
    
    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point
    
    pre_sp <- predictSurvProb2survreg(obj, test_data, t_star0)
    object <- Surv(test_data$time, test_data$status)
    t_star <- t_star0
  }
  
  if (is.na(t_star)) {
    stop("Cannot calculate Brier Score at NA")
  }
  
  # default time point for Brier(t_star)
  if (t_star <= 0) {
    t_star <- median(object[, 1][object[, 2] == 1])
  }
  
  if (length(t_star) != 1) {
    stop("Brier Score can only be calculated at a single time point")
  }
  
  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }
  
  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }
  
  if (any(is.na(pre_sp))) {
    stop("The input probability vector cannot have NA")
  }
  
  if (length(object) != length(pre_sp)) {
    stop("The prediction survival probability and the survival object have different lengths")
  }
  
  time <- object[, 1]
  status <- object[, 2]
  
  t_order <- order(time)
  time <- time[t_order]
  status <- status[t_order]
  pre_sp <- pre_sp[t_order]
  
  # Initialization
  sum_before_t <- 0
  sum_after_t <- 0
  
  # Standard IPCW weighting (Graf et al. 1999): subjects with an event before
  # t_star are weighted by the LEFT limit G(T_i-), subjects still at risk at
  # t_star by G(t_star) itself. Both are step-function evaluations of the
  # censoring Kaplan-Meier; neither interpolates between jump times.
  Gtstar <- Gt(object, t_star)
  for (i in c(1:length(time))) {
    # survival time is less than t_star and sample died
    if (time[i] < t_star & (status[i] == 1)) {
      Gti <- Gt(object, time[i], left = TRUE)
      if (is.na(Gti)) {
        next
      }
      sum_before_t <- sum_before_t + 1 / Gti * (pre_sp[i]) ^ 2 # IPCW
      next
    }
    # survival time is greater than t_star
    if (time[i] >= t_star) {
      if (is.na(Gtstar)) {
        next
      }
      sum_after_t <- sum_after_t + 1 / Gtstar * (1 - pre_sp[i]) ^ 2
    } # IPCW
  }
  
  BSvalue <- (sum_before_t + sum_after_t) / length(time)
  names(BSvalue) <- "Brier Score"
  
  return(round(BSvalue, 6))
}