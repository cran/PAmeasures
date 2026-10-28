#' Simulate Cox survival data with optional interaction terms and targeted censoring
#'
#' @description
#' Generates survival data under a Cox proportional hazards model with a Weibull
#' baseline cumulative hazard \eqn{H_0(t) = (0.5 t)^v}. Event times are sampled via
#' inverse transform:
#' \deqn{Y = H_0^{-1}(-\log U \cdot \exp(-\beta^\top Z)), \quad U \sim \mathrm{Unif}(0,1),}
#' where \eqn{H_0^{-1}(s) = 2 s^{1/v}}. 
#' 
#' Covariates are simulated from a standard normal distribution.
#' When \code{interact = FALSE}, the design matrix \eqn{Z} contains only main effects.
#' When \code{interact = TRUE}, the function automatically detects whether \code{length(beta)} 
#' equals 3 or 7:
#' \itemize{
#'   \item If \code{length(beta) = 3}, two covariates (\eqn{x_1, x_2}) and their 
#'         two-way interaction (\eqn{x_1 x_2}) are included.
#'   \item If \code{length(beta) = 7}, three covariates (\eqn{x_1, x_2, x_3}), 
#'         their three two-way interactions (\eqn{x_1 x_2, x_1 x_3, x_2 x_3}), 
#'         and their three-way interaction (\eqn{x_1 x_2 x_3}) are included.
#' }
#'
#' Right-censoring times are generated as \eqn{C_i = \exp(\mu + r_i)} with 
#' \eqn{r_i \sim \mathcal{N}(0, \mathrm{sd}^2)}. By default, \code{sd} is set to 
#' \eqn{\mathrm{sd}(\log Y)}, and \code{mu} is calibrated so that the expected censoring 
#' proportion equals \code{pi_c}. If \code{mu} and/or \code{sd} are supplied, calibration 
#' is skipped accordingly:
#' \itemize{
#'   \item \strong{mu and sd provided:} use them directly (no calibration).
#'   \item \strong{sd only:} solve for \code{mu} using the given \code{sd}.
#'   \item \strong{mu only:} use the given \code{mu} and estimate \eqn{\mathrm{sd}(\log Y)}.
#'   \item \strong{neither:} estimate both \eqn{\mathrm{sd}(\log Y)} and \code{mu}.
#' }
#'
#' @param n Integer. Sample size.
#' @param pi_c Numeric in \eqn{[0, 0.99]}. Target censoring proportion. Ignored if both
#'   \code{mu} and \code{sd} are provided.
#' @param v Positive numeric. Weibull shape in \eqn{H_0(t) = (0.5 t)^v}.
#' @param beta Numeric vector of regression coefficients. 
#'   If \code{interact = FALSE}, its length defines the number of main effects. 
#'   If \code{interact = TRUE}, its length must be either 3 (two main + one interaction) 
#'   or 7 (three main + four interactions).
#' @param mu Optional numeric. Mean of \eqn{\log C}. If supplied with \code{sd}, 
#'   no calibration is performed.
#' @param sd Optional positive numeric. Standard deviation used in 
#'   \eqn{\log C = \mu + r_i}. If \code{NULL}, defaults to \eqn{\mathrm{sd}(\log Y)}.
#' @param seed Optional integer seed for reproducibility.
#' @param nonlinear Logical. If `TRUE`, generates two covariates plus their
#'   interaction and a quadratic term; `beta` must then have length 4.
#' @param interact Logical (default \code{FALSE}). 
#'   If \code{TRUE}, include pre-defined interaction terms as described above.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{time}: observed time \eqn{T_i = \min(Y_i, C_i)}.
#'   \item \code{status}: event indicator \eqn{1\{Y_i \le C_i\}}.
#'   \item \code{x1, x2, ...}: simulated covariates.
#'   \item \code{x1:x2, x1:x3, ...}: interaction terms (if \code{interact = TRUE}).
#'   \item \code{y_true}: latent event time \eqn{Y_i}.
#'   \item \code{cens_time}: censoring time \eqn{C_i}.
#' }
#'
#' Attributes include:
#' \itemize{
#'   \item \code{mu}: censoring shift used.
#'   \item \code{sd_log}: standard deviation used in \eqn{\log C}.
#'   \item \code{sd_logY}: \eqn{\mathrm{sd}(\log Y)} from the generated events.
#'   \item \code{pi_c_target}: requested censoring rate.
#'   \item \code{pi_c_observed}: realized censoring rate.
#'   \item \code{v}, \code{beta}: model parameters.
#'   \item \code{interact}: logical flag indicating whether interactions were included.
#'   \item \code{nonlinear}: logical flag indicating whether polynomial and interactions term were included.
#' }
#'
#' @examples
#' # (1) Main effects only
#' set.seed(1)
#' d0 <- pa_sim_cox_weibull(500, 0.3, v = 1.2, beta = c(0.4, -0.2))
#'
#' # (2) Two covariates + their interaction
#' d1 <- pa_sim_cox_weibull(500, 0.3, v = 1.2,
#'                                beta = c(0.5, -0.3, 0.2), 
#'                                interact = TRUE)
#'
#' # (3) Three covariates + all two-way and three-way interactions
#' b <- c(0.4, -0.2, 0.1, 0.15, -0.05, 0.08, 0.02)
#' d2 <- pa_sim_cox_weibull(500, 0.3, v = 1.0, beta = b, interact = TRUE)
#'
#' @export
pa_sim_cox_weibull <- function(n, pi_c, v, beta, mu = NULL,
                                     sd = NULL, seed = NULL, 
                                     interact = FALSE, nonlinear = FALSE) {
  stopifnot(length(n) == 1L, n > 0, is.finite(n))
  stopifnot(length(v) == 1L, v > 0, is.finite(v))
  stopifnot(length(pi_c) == 1L, is.finite(pi_c), pi_c >= 0, pi_c < 1)
  stopifnot(is.numeric(beta), length(beta) >= 1L, all(is.finite(beta)))
  if (!is.null(mu)) stopifnot(length(mu) == 1L, is.finite(mu))
  if (!is.null(sd)) stopifnot(length(sd) == 1L, is.finite(sd), sd > 0)
  
  if (!is.null(seed)) set.seed(seed)
  
  # ----- Generate covariates -----
  if ((!interact) & (!nonlinear)) {
    p <- length(beta)
    X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
    colnames(X) <- paste0("x", seq_len(p))
    Z <- X
  } else if(nonlinear){
    if (length(beta) == 4) {
      X <- matrix(rnorm(n * 2, mean = 0, sd = 1), nrow = n, ncol = 2)
      colnames(X) <- paste0("x", 1:2)
      Z <- cbind(
        X,
        "x1:x2" = X[,1] * X[,2],
        "x1^2" = X[,1] * X[,1]
      )
    } else {
      stop("When nonlinear = TRUE, length(beta) must be 4 (X1, X2, x1:x2, X2^2).")
    }
  } else {
    if (length(beta) == 3) {
      # 2 covariates + 1 two-way interaction
      X <- matrix(rnorm(n * 2, mean = 0, sd = 1), nrow = n, ncol = 2)
      colnames(X) <- c("x1", "x2")
      Z <- cbind(X, "x1:x2" = X[,1] * X[,2])
    } else if (length(beta) == 7) {
      # 3 covariates + all two-way + three-way interactions
      X <- matrix(rnorm(n * 3, mean = 0, sd = 1), nrow = n, ncol = 3)
      colnames(X) <- c("x1", "x2", "x3")
      Z <- cbind(
        X,
        "x1:x2" = X[,1] * X[,2],
        "x1:x3" = X[,1] * X[,3],
        "x2:x3" = X[,2] * X[,3],
        "x1:x2:x3" = X[,1] * X[,2] * X[,3]
      )
    } else {
      stop("When interact = TRUE, length(beta) must be 3 (2 vars + 1 interaction) or 7 (3 vars + 4 interactions).")
    }
  }
  
  # Event times: Y = H0^{-1}(-log U * exp(-eta)), H0^{-1}(s) = 2 * s^{1/v}
  U <- stats::runif(n)
  eta <- as.numeric(Z %*% beta)
  s <- -log(U) * exp(-eta)
  Y <- 2 * (s)^(1 / v)
  
  logY <- log(Y)
  sd_logY <- stats::sd(logY)
  if (!is.finite(sd_logY) || sd_logY <= .Machine$double.eps) sd_logY <- 1e-6
  
  if(pi_c == 0){
    out <- data.frame(
      time = Y,
      status = 1,
      X
    )
    return(out)
  }
  # Determine sd to use for log-censoring noise
  sd_use <- if (is.null(sd)) sd_logY else sd
  
  # Determine mu to use
  mu_use <- mu
  if (is.null(mu_use)) {
    # Solve for mu so that mean(Phi((logY - mu)/sd_use)) = 1 - pi_c
    target_obs <- 1 - pi_c
    f_mu <- function(m) {
      mean(stats::pnorm((m - logY) / sd_use)) - target_obs
    }
    lo <- min(logY) - 10 * sd_use - 5
    hi <- max(logY) + 10 * sd_use + 5
    f_lo <- f_mu(lo); f_hi <- f_mu(hi)
    expand <- 0
    while (f_lo * f_hi > 0 && expand < 5) {
      lo <- lo - 10; hi <- hi + 10
      f_lo <- f_mu(lo); f_hi <- f_mu(hi)
      expand <- expand + 1
    }
    mu_use <- tryCatch(
      stats::uniroot(f_mu, lower = lo, upper = hi)$root,
      error = function(e) stats::quantile(logY, probs = target_obs, names = FALSE)
    )
  }
  
  # Generate censoring and observed
  r <- stats::rnorm(n, mean = 0, sd = sd_use)
  C <- exp(mu_use + r)
  
  time <- pmin(Y, C)
  status <- as.integer(Y <= C)
  pi_c_obs <- mean(status == 0)
  
  out <- data.frame(
    time = time,
    status = status,
    X,
    y_true = Y,
    cens_time = C,
    check.names = FALSE
  )
  
  attr(out, "mu") <- as.numeric(mu_use)
  attr(out, "sd_log") <- as.numeric(sd_use)
  attr(out, "sd_logY") <- as.numeric(sd_logY)
  attr(out, "pi_c_target") <- as.numeric(pi_c)
  attr(out, "pi_c_observed") <- as.numeric(pi_c_obs)
  attr(out, "v") <- as.numeric(v)
  attr(out, "beta") <- as.numeric(beta)
  
  out
}

