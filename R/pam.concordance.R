#' @title Compute the Concordance Statistic
#'
#' @description This function computes the concordance statistic, which measures the agreement
#' between an observed response and a predictor. It is closely related to Kendall's
#' tau-a and tau-b, Goodman's gamma, and Somers' d, all of which can also be calculated
#' from the results of this function. This function handles different model types and
#' can be used with new data to evaluate but not refit the models.
#'
#' @param object A fitted model or a formula of the form y ~ x or y ~ x + strata(z) with
#' a single numeric or survival response and a single predictor. Counts of concordant,
#' discordant and tied pairs are computed separately per stratum and then added.
#' @param data A data.frame used to interpret variables named in the formula. Only applicable
#' if object is a formula.
#' @param newdata Optional new data frame for evaluating models.
#' @param weights Optional vector of case weights, applicable only if object is a formula.
#' @param subset An expression indicating which subset of the rows of data should be used
#' in the fit. Only applicable if object is a formula.
#' @param na.action A missing-data filter function, applied to the model.frame after any
#' subset argument has been used. Default is options()$na.action. Only applicable if object
#' is a formula.
#' @param cluster Optional grouping vector for calculating robust variance.
#' @param ymin Optional lower bound for y values in the calculation.
#' @param ymax Optional upper bound for y values in the calculation.
#' @param timewt Weighting to be applied. The overall statistic is a weighted mean over
#' event times, with options like "n", "S", "S/G", "n/G2", "I".
#' @param influence Level of influence output to return, where 1=return the dfbeta vector,
#' 2=return the full influence matrix, 3=return both.
#' @param ranks If TRUE, returns a data frame with ranks.
#' @param reverse If TRUE, assume larger x values predict smaller response values y.
#' @param timefix If TRUE, correct for possible rounding errors related to tied times.
#' @param keepstrata Either TRUE, FALSE, or an integer value. Computations are always done
#' within stratum, then added. If the total number of strata is greater than keepstrata, or
#' keepstrata=FALSE, those subtotals are not kept in the output.
#' @return An object of class 'concordance' containing the estimated concordance value or values,
#' counts of concordant, discordant, and tied pairs, number of observations, estimated variance
#' of the concordance, and optionally, the data frame of ranks, dfbeta vector, and influence matrix.
#' 
#' @references
#' F Harrell, R Califf, D Pryor, K Lee and R Rosati, Evaluating the yield of medical tests, J Am Medical Assoc, 1982.
#' 
#' R Peto and J Peto, Asymptotically efficient rank invariant test procedures (with discussion), J Royal Stat Soc A, 1972.
#' 
#' M Schemper, Cox analysis of survival data with non-proportional hazard functions, The Statistician, 1992.
#' 
#' H Uno, T Cai, M Pencina, R D'Agnostino and Lj Wei, On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data, Statistics in Medicine, 2011.
#' 
#' Therneau, T. M., Lumley, T., Atkinson, E., Crowson, C. (2024). survival: Survival Analysis. 
#' R package version 3.7-0. DOI: \doi{10.32614/CRAN.package.survival}. Available at \url{https://CRAN.R-project.org/package=survival}.
#' 
#' @examples
#' library(survival)
#' library(dplyr)
#'
#' # Load the pbc dataset
#' data(pbc)
#'
#' # Data preparation
#' pbc <- pbc %>%
#'   filter(!is.na(trt)) %>%
#'   mutate(
#'     log_albumin = log(albumin),
#'     log_bili = log(bili),
#'     log_protime = log(protime),
#'     status = ifelse(status == 2, 1, 0)
#'   )
#'
#' # Fit a Cox Proportional Hazards model
#' cox_model <- cph(
#'   Surv(time, status) ~ age + log_albumin + log_bili + log_protime + edema,
#'   data = pbc,
#'   x = TRUE,
#'   y = TRUE
#' )
#'
#' # Compute the concordance statistic
#' concordance_result <- pam.concordance(cox_model)
#' print(concordance_result)                   
#' @keywords internal
#' @noRd

pam.concordance <- function(object, ...) {
  survival::concordance(object, ...)
}
