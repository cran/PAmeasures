# Deterministic fixtures shared by all characterization tests.
# Seeds are fixed; never change them without regenerating every snapshot.

# survival is deliberately NOT attached here. PAmeasures calls concordancefit()
# unqualified without importing it (findings.md #11), and a global attach would
# hide that defect from the whole suite. Tests that genuinely need it use
# withr::local_package("survival"), which scopes the attach to one test and
# detaches on exit. test-eval-survival.R reproduces the failure in a clean
# subprocess.

fx_covs <- function() c("x1", "x2")

fx_surv <- function() {
  d <- pa_sim_cox_weibull(
    n = 200, pi_c = 0.3, v = 2, beta = c(0.5, -0.5), seed = 1001
  )
  # pi_c > 0 also returns y_true and cens_time; drop them so model
  # formulas built with "." cannot pick them up as covariates.
  d[, c("time", "status", "x1", "x2")]
}

fx_surv_uncensored <- function() {
  pa_sim_cox_weibull(
    n = 200, pi_c = 0, v = 2, beta = c(0.5, -0.5), seed = 1002
  )[, c("time", "status", "x1", "x2")]
}

fx_cr <- function() {
  pa_simulate_fine_gray(
    n = 200, v = 2, beta1 = c(0.5, -0.5), beta2 = c(-0.3, 0.3),
    censor = 0.3, seed = 1003
  )
}

fx_cox <- function() {
  survival::coxph(
    survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), x = TRUE, y = TRUE
  )
}

fx_survreg <- function() {
  # x = TRUE, y = TRUE is required: pa_predict_survreg errors without it
  survival::survreg(
    survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), dist = "weibull", x = TRUE, y = TRUE
  )
}
