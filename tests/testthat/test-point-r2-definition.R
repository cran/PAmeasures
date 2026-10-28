test_that("CR point R2 includes a target event exactly at the horizon", {
  time <- c(1, 2, 3, 5, 6, 7)
  status <- c(1, 2, 1, 1, 2, 1)
  probability <- c(.85, .1, .7, .6, .2, .15)
  grid <- c(1, 3, 5, 7)
  cif <- outer(c(.2, .5, 1, 1), probability)
  # Independent equal-weight moments for I(T <= 5, D = 1), not I(T < 5, D = 1).
  expected <- c(r2_point = .9293, l2_point = .3076, pseudo_r2_point = .2859)
  for (method in c("km", "cox")) {
    actual <- pa_survival_eval_cr(cif, time, grid, status,
      tau = 5, t_star = 5, censoring = method, metrics = names(expected))
    expect_equal(actual$Value, unname(expected[actual$Metric]))
    prediction <- list(cif_pred = cbind(grid, cif), times = time, status = status)
    summary <- pa_summarize_cr(list(model = prediction), tau = 5, t_star = 5,
      censoring = method, metrics = names(expected), digits = 4)
    expect_equal(summary$model, unname(expected[summary$Metric]))
  }
})

test_that("ordinary point R2 compares the event indicator with F rather than S", {
  time <- c(1, 2, 3, 6, 7, 8)
  probability <- c(.85, .7, .6, .2, .15, .1)
  # Constant time columns isolate probability direction from grid selection.
  survival <- matrix(rep(1 - probability, length(time)), nrow = length(time))
  expected <- c(r2_point = .9293, l2_point = .3076, pseudo_r2_point = .2859)
  actual <- pa_survival_eval(model = NULL, event_time = time,
    predicted_probability = survival, pred_mean_survival = rep(4, length(time)),
    status = rep(1, length(time)), covariates = NULL, tau = 4, t_star = 4,
    metrics = names(expected))
  expect_equal(actual$Value, unname(expected[actual$Metric]))
  prediction <- list(times = time, status = rep(1, length(time)),
    surv_prob = survival, pred = rep(4, length(time)))
  summary <- pa_summarize(list(model = prediction), tau = 4, t_star = 4,
    metrics = names(expected), digits = 4)
  expect_equal(summary$model, unname(expected[summary$Metric]))
})

test_that("ordinary point R2 includes an event at the evaluation time", {
  time <- c(1, 2, 4, 6, 7, 8)
  probability <- c(.85, .7, .6, .2, .15, .1)
  survival <- matrix(rep(1 - probability, length(time)), nrow = length(time))
  actual <- pa_survival_eval(model = NULL, event_time = time,
    predicted_probability = survival, pred_mean_survival = rep(4, length(time)),
    status = rep(1, length(time)), covariates = NULL, tau = 4, t_star = 4,
    metrics = c("r2_point", "l2_point", "pseudo_r2_point"))
  expected <- c(r2_point = .9293, l2_point = .3076, pseudo_r2_point = .2859)
  expect_equal(actual$Value, unname(expected[actual$Metric]))
})
