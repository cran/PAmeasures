#' Find the mean of censoring times
#'
#' This function is an internal helper used to find the shift `c_mu` needed
#' to achieve a desired censoring proportion.
#'
#' @param c_mu A numeric value representing the shift to be applied to censoring times.
#' @param censor A numeric value representing the desired censoring proportion.
#' @param cens.times A vector of censoring times.
#' @param event.times A vector of event times.
#' @param f.event A vector representing the event type.
#' @param n The number of observations.
#' @return A numeric value representing the difference between the observed censoring
#' proportion and the desired censoring proportion. This is used for root-finding.
#' @keywords internal
#' @noRd
find_mu_c <- function(c_mu, censor, cens.times, event.times, f.event, n) {
  # Shift censoring times by the provided mu_c value
  cens.times <- cens.times + c_mu
  
  # Determine observed times and events based on censoring
  obs.times <- round(pmin(event.times, cens.times), 7)
  obs.event <- (event.times <= cens.times) * f.event
  
  # Calculate the difference between the observed censoring proportion and the desired one
  result <- mean(obs.event == 0) - censor
  
  # Handle cases where the result is not a single scalar
  if (length(result) > 1) {
    warning("Result is not a scalar, taking mean.")
    result <- mean(result)
  }
  return(result)
}


#' Simulate data from a two-cause Fine-Gray model
#'
#' This function generates a simulated dataset for a competing risks model with two causes,
#' following a Fine-Gray-like model specification.
#'
#' @importFrom magrittr %>%
#'
#' @param n An integer specifying the number of subjects to simulate.
#' @param v A numeric value for the Weibull shape parameter.
#' @param beta1 A numeric vector of coefficients for cause 1.
#' @param beta2 A numeric vector of coefficients for cause 2.
#' @param lambda1 A numeric value for the Weibull scale parameter for cause 1.
#' @param X A matrix of covariates, if not provided, covariates are generated from a
#'   standard normal distribution.
#' @param mu A numeric value for the mean of the generated covariates.
#' @param p A numeric value for the marginal probability of a cause 1 event.
#' @param c_scale A numeric value to scale the standard deviation of censoring times.
#' @param censor A numeric value representing the desired censoring proportion.
#' @param sd.time A numeric value for the standard deviation of censoring times.
#'   If NULL, it's calculated based on event times.
#' @param mu.c.time A numeric value for the mean of censoring times. If NULL and `censor > 0`,
#'   it is calculated using `uniroot` to achieve the desired censoring.
#' @param independent_c A logical value. If TRUE, censoring times are independent of covariates.
#' @param report.mu_and_sd A logical value. If TRUE, a summary of censoring parameters
#'   and event time quantiles is returned instead of the full dataset.
#' @param seed Optional integer. If supplied, the random number generator is
#'   seeded with it for reproducibility. The default \code{NULL} uses the current random number stream
#'   without resetting it; random draws advance its state.
#'
#' @return A data frame containing the simulated data with columns `obs.times`, `obs.event`, and the covariates `X`.
#'   Alternatively, if `report.mu_and_sd` is TRUE, it returns a data frame with censoring
#'   parameters and event time quantiles.
#'
#' @examples
#' # Simulate a dataset with 100 subjects
#' sim_data <- pa_simulate_fine_gray(
#'   n = 100, v = 1, beta1 = c(0.5, -0.2), beta2 = c(-0.3, 0.4)
#' )
#' head(sim_data)
#'
#' # Simulate with a specific censoring proportion and report parameters
#' params_summary <- pa_simulate_fine_gray(
#'   n = 1000, v = 1, beta1 = c(0.5, -0.2), beta2 = c(-0.3, 0.4),
#'   censor = 0.2, report.mu_and_sd = TRUE
#' )
#' params_summary
#' @export
pa_simulate_fine_gray <- function (n, v, beta1, beta2, lambda1 = 1,
                                           X = NULL, mu = 0, p = 0.7,
                                           c_scale = 1, censor = 0,
                                           sd.time = NULL, mu.c.time = NULL,
                                           independent_c = TRUE,
                                           report.mu_and_sd = FALSE,
                                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta1) != length(beta2))
    stop("Dimension of beta1 and beta2 should be the same")
  ncovs <- length(beta1)
  if (is.null(X)) {
    X <- matrix(rnorm(n * ncovs, mean = mu), nrow = n)
  }
  c.ind <- 1 + rbinom(n, 1, prob = (1 - p)^exp(X %*% beta1))
  ftime <- numeric(n)
  eta1 <- X[c.ind == 1, ] %*% beta1
  eta2 <- X[c.ind == 2, ] %*% beta2
  u1 <- runif(length(eta1))
  t1 <- (-log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^
                     (1/exp(eta1)))/p)/lambda1)^(1/v)
  t2 <- rweibull(length(eta2), shape = v, scale = exp(eta2)^(-1/v))
  
  if(is.null(sd.time)) sd.time <- sd(log(t1)) * c_scale
  if(independent_c==TRUE){
    cens.times <- rnorm(n, mean = 0, sd = sd.time)
  }else{
    temp <- apply(X, 1, sum)
    cens.times <- sapply(temp, function(x){
      rnorm(n = 1, mean = x, sd = sd.time)})
  }
  
  ftime[c.ind == 1] <- t1
  ftime[c.ind == 2] <- t2
  
  if(censor>0){
    if(is.null(mu.c.time)){
      result <- uniroot(find_mu_c,
                        interval = c(-1e5, 1e5),
                        censor = censor,
                        cens.times = cens.times,
                        event.times = log(ftime),
                        f.event = c.ind,
                        n = n)
      mu.c.time <- result$root
    }
    # This requires the `magrittr` package for the %>% operator.
    # It is good practice to explicitly import it in the Roxygen comments with @importFrom magrittr %>%
    obs.times <- pmin(log(ftime), cens.times) %>% exp
    obs.event <- c(log(ftime) <= cens.times) * c.ind
  }else{
    mu.c.time <- NA
    obs.times <- ftime
    obs.event <- c.ind
  }
  
  if(report.mu_and_sd == FALSE){
    final_data <- data.frame(obs.times, obs.event, X)
    return(final_data)
  }else{
    mu_and_sd <- data.frame(sd.time, mu.c.time,
                            t(quantile(obs.times, c(0.1, 0.3, 0.5, 0.7, 0.9))))
    return(mu_and_sd)
  }
}
