weighted_calibration_data <- function() {
  withr::local_seed(1929)
  n <- 180
  rate <- .12 * exp(.5 * stats::rnorm(n))
  event <- stats::rexp(n, rate)
  censor <- stats::rexp(n, .06)
  list(time = pmin(event, censor), status = as.integer(event <= censor),
       p = 1-exp(-rate*3), weight = rep(c(1, 3, 8), length.out = n))
}

test_that("unit weights reproduce ordinary right-censored calibration", {
  d <- weighted_calibration_data()
  ordinary <- pa_calibration(d$p, d$time, d$status, 3)
  weighted <- pa_calibration(d$p, d$time, d$status, 3,
                              case_weights = rep(1, length(d$p)))
  expect_equal(weighted$curve, ordinary$curve, tolerance = 1e-9)
  expect_equal(weighted$metrics, ordinary$metrics, tolerance = 1e-9)
  expect_identical(weighted$method, "weighted_cox_rcs")
  expect_length(weighted$knots, 5L)
})

test_that("sampling weights affect both curve fitting and error summaries", {
  d <- weighted_calibration_data()
  ordinary <- pa_calibration(d$p, d$time, d$status, 3)
  result <- pa_calibration(d$p, d$time, d$status, 3, case_weights = d$weight)
  error <- abs(d$p - result$curve$calibrated)
  expect_equal(unname(result$metrics['ici']), sum(d$weight*error)/sum(d$weight))
  expect_equal(unname(result$metrics[c('e50','e90')]),
    unname(Hmisc::wtd.quantile(error, d$weight, probs = c(.5,.9), normwt = TRUE)))
  expect_gt(max(abs(result$curve$calibrated-ordinary$curve$calibrated)), 1e-4)
  posthoc <- sum(d$weight*abs(d$p-ordinary$curve$calibrated))/sum(d$weight)
  expect_gt(abs(unname(result$metrics['ici'])-posthoc), 1e-4)
  expect_identical(result$case_weights, d$weight)

  scaled <- pa_calibration(d$p, d$time, d$status, 3, case_weights = 17*d$weight)
  expect_equal(scaled$metrics, result$metrics, tolerance = 1e-10)
  expect_equal(scaled$curve, result$curve, tolerance = 1e-10)
  o <- rev(seq_along(d$p))
  permuted <- pa_calibration(d$p[o], d$time[o], d$status[o], 3,
                              case_weights = d$weight[o])
  expect_equal(permuted$metrics, result$metrics, tolerance = 1e-9)
  expect_equal(permuted$curve$calibrated, rev(result$curve$calibrated),
                tolerance = 1e-9)
})

test_that("weighted constant-risk calibration agrees with hand KM risk sets", {
  time <- 1:6
  status <- c(1, 0, 1, 0, 1, 0)
  w <- c(1, 2, 4, 1, 3, 2)
  # At 4: the weighted KM survival is (1-1/13)*(1-4/10).
  risk <- 1-(1-1/13)*(1-4/10)
  result <- pa_calibration(rep(.2, 6), time, status, 4, case_weights = w)
  expect_equal(result$curve$calibrated, rep(risk, 6))
  expect_equal(unname(result$metrics), rep(risk-.2, 3))
  expect_identical(result$method, "weighted_kaplan_meier_constant")
  expect_length(result$knots, 0)
})

test_that("weighted Cox risk matches independently solved score and baseline", {
  # Two risk values give one calibration coefficient. Unique event times let
  # us solve the weighted partial-likelihood score and baseline without Cox.
  time <- 1:12
  status <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0)
  group <- rep(0:1, 6)
  w <- c(1, 2, 3, 1, 4, 2, 1, 5, 2, 1, 3, 4)
  events <- which(status == 1)
  score <- function(beta) {
    risk <- w*exp(beta*group)
    sum(vapply(events, function(i) {
      at <- time >= time[i]
      w[i]*(group[i]-sum(risk[at]*group[at])/sum(risk[at]))
    }, numeric(1)))
  }
  beta <- uniroot(score, c(-10,10), tol = 1e-12)$root
  h0 <- sum(vapply(events[time[events] <= 7.5], function(i) {
    w[i]/sum(w[time >= time[i]]*exp(beta*group[time >= time[i]]))
  }, numeric(1)))
  expected <- -expm1(-h0*exp(beta*group))
  p <- ifelse(group == 0, .2, .6)
  result <- pa_calibration(p, time, status, 7.5, case_weights = w)
  expect_equal(result$curve$calibrated, expected, tolerance = 1e-8)
  expect_equal(unname(result$metrics['ici']), sum(w*abs(p-expected))/sum(w),
                tolerance = 1e-8)
  expect_identical(result$method, "weighted_cox_linear")
})

test_that("two-phase calibration uses aligned weights and the past grid value", {
  d <- weighted_calibration_data()
  grid <- sort(d$time)
  rate <- -log1p(-d$p)/3
  surv <- outer(rate, grid, function(a,t) exp(-a*t))
  t <- stats::median(grid)
  p <- 1-surv[, findInterval(t, grid)]
  reference <- pa_calibration(p, d$time, d$status, t,
                              nk = 3, case_weights = d$weight)$metrics
  # No model, mean time prediction or separate censoring fit is needed.
  pred <- list(times = d$time, status = d$status, surv_prob = surv)
  result <- pa_evaluate_two_phase(pred, case_weights = d$weight,
    metrics = c('ici','e50','e90'), t_star = t, calibration_args = list(nk = 3))
  expect_equal(setNames(result$Value, result$Metric), round(reference, 4))
  summary <- pa_sample_design(list(external = pred), case_weights = d$weight,
    metrics = c('ici','e50','e90'), t_star = t, calibration_args = list(nk = 3),
    digits = 4)
  expect_equal(summary$external, unname(round(reference, 4)))
  expect_identical(summary$Metric, c('ici','e50','e90'))
  o <- rev(seq_along(d$time))
  other <- list(times = d$time[o], status = d$status[o],
                surv_prob = surv[o, , drop = FALSE])
  permuted <- pa_evaluate_two_phase(other, case_weights = d$weight[o],
    metrics = c('ici','e50','e90'), t_star = t, calibration_args = list(nk = 3))
  expect_equal(permuted, result)
})

test_that("sampling calibration validates weights and preserves boundary risks", {
  d <- weighted_calibration_data()
  for (w in list(d$weight[-1], replace(d$weight, 1, NA_real_),
                 replace(d$weight, 1, Inf), replace(d$weight, 1, 0),
                 replace(d$weight, 1, -1), matrix(d$weight, ncol = 1))) {
    expect_error(pa_calibration(d$p, d$time, d$status, 3, case_weights = w),
                   'positive finite sampling-weight vector')
  }
  p <- d$p; p[1:2] <- c(0, 1)
  expect_warning(result <- pa_calibration(p, d$time, d$status, 3,
                                          case_weights = d$weight),
                   'Boundary probabilities')
  expect_identical(result$curve$predicted, p)
  expect_equal(unname(result$metrics['ici']),
    sum(d$weight*abs(p-result$curve$calibrated))/sum(d$weight))
  pred <- list(times = d$time, status = d$status,
                surv_prob = outer(d$p, sort(d$time), function(p,t) (1-p)^(t/3)))
  expect_error(pa_evaluate_two_phase(pred, case_weights = d$weight,
    metrics = 'ici', calibration_args = list(case_weights = d$weight)),
    "supply sampling weights")
  expect_warning(unavailable <- pa_evaluate_two_phase(pred,
    case_weights = d$weight, metrics = 'ici', calibration_args = list(nk = 2)),
    'Calibration unavailable')
  expect_true(is.na(unavailable$Value))
})

test_that("ordinary evaluation cannot silently misalign sampling weights", {
  d <- weighted_calibration_data()
  pred <- outer(d$p, sort(d$time), function(p,t) (1-p)^(t/3))
  expect_error(pa_survival_eval(event_time = d$time, status = d$status,
    predicted_probability = pred, metrics = 'ici',
    calibration_args = list(case_weights = d$weight)),
    'For sampling-adjusted calibration use pa_sample_design')
  expect_error(pa_fit_and_eval(fx_surv(), fx_covs(), metrics = 'ici',
    calibration_args = list(case_weights = rep(1, nrow(fx_surv())))),
    'For sampling-adjusted calibration use pa_sample_design')
})
