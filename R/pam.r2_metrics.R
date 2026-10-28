#' Compute R^2, L^2, and Pseudo-R^2 for Survival Data
#'
#' This function calculates R^2, L^2, and pseudo-R^2
#' for survival data using observed and predicted survival times.
#'
#' @param predicted_data A numeric vector of predicted survival times.
#' @param survival_time A numeric vector of observed survival times.
#' @param status A numeric vector indicating the event occurrence (1 for event, 0 for censoring).
#' @param tau An optional numeric value for restricted time horizon. Default is NULL (no restriction).
#' @param case_weight An optional numeric value for case specific weight. Default is NULL.
#'
#' @return A list with R^2 , L^2, and pseudoR^2.
#'
#' @examples
#' 
#' predicted_data <- c(5, 4, 8, 2)
#' survival_time <- c(6, 5, 10, 3)
#' status <- c(1, 1, 0, 1)
#' pam.r2_metrics(predicted_data, survival_time, status)
#' @keywords internal
#' @noRd
pam.r2_metrics <- function(predicted_data, survival_time, status, 
                           tau = NULL, case_weight = NULL) {
  
  # Apply restriction if tau is provided
  if (!is.null(tau)) {
    restricted <- restricted_data_gen(survival_time, status, tau)
    survival_time <- restricted$time
    status <- restricted$status
  }
  
  # Sort survival_time and status
  y.order <- order(survival_time)
  y <- survival_time[y.order]
  delta <- status[y.order]
  predicted_data <- predicted_data[y.order]
  if (!is.null(case_weight)) case_weight <- case_weight[y.order]
  
  # Calculate inverse probability of censoring weights (IPCW)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  km.censoring.minus <- pa_censoring_left_limit(fit.km.censoring, y)
  ratio.km <- delta / km.censoring.minus
  
  # Handle NaN and Inf values
  ratio.km[is.nan(ratio.km)] <- 0
  ratio.km[is.infinite(ratio.km)] <- 0
  
  # Calculate weights
  # weight.km <- ratio.km / sum(ratio.km)
  if (!is.null(case_weight)) {
    ratio.km.new <- ratio.km * case_weight
    #weight.new <- (case_weight / sum(case_weight)) * weight.km
  }else{
    ratio.km.new <- ratio.km
  }
  weight.new <- ratio.km.new / sum(ratio.km.new)
  
  
  # Fit weighted least squares (WLS) regression
  wls.fitted <- tryCatch({
    lm(y ~ predicted_data, weights = weight.new)
  }, error = function(e) {
    warning("WLS regression failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(wls.fitted)) {
    return(list(R.squared = NA, L.squared = NA, Psuedo.R = NA))
  }
  
  # Calculate fitted values
  calibrate.fitted <- predict(wls.fitted)
  
  # Calculate R^2 components
  weighted_mean_y <- sum(weight.new * y)
  num.rho2 <- sum(weight.new * (calibrate.fitted - weighted_mean_y)^2)
  denom.rho2 <- sum(weight.new * (y - weighted_mean_y)^2)
  R2 <- round(num.rho2 / denom.rho2, digits = 4)
  
  # Calculate L^2 components
  num.L2 <- sum(weight.new * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.new * (y - predicted_data)^2)
  L2 <- round(num.L2 / denom.L2, digits = 4)
  
  # Calculate Psuedo R^2
  SR <- round(R2 * L2, digits = 4)
  
  # Return results
  return(list(R_squared = R2, L_squared = L2, Pseudo_R_squared = SR))
}
# Helper function to restrict survival data
restricted_data_gen <- function(time, status, tau) {
  output <- data.frame(
    time = ifelse(time > tau, tau, time),
    status = ifelse(time > tau, 1, status)
  )
  return(output)
}
