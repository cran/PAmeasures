simulation_calibration_helpers <- function() {
  root <- skip_without_source_tree()
  path <- file.path(root, "analysis", "simulation_calibration.R")
  skip_if_not(file.exists(path), "simulation helpers are excluded from installed builds")
  e <- new.env(parent = asNamespace("PAmeasures"))
  sys.source(path, envir = e)
  e
}

test_that("simulation Cox calibration uses the exact horizon in subject order", {
  e <- simulation_calibration_helpers()
  fit <- fx_cox()
  d <- fx_surv_uncensored()[seq(1, 200, by = 5), ]
  pred <- pa_predict_coxph(fit, fx_covs(), new_data = d)
  # Choose an actual training jump strictly between validation grid points.
  baseline <- survival::basehaz(fit)
  candidates <- baseline$time[baseline$hazard > c(0, head(baseline$hazard, -1))]
  candidates <- candidates[candidates > min(d$time) & candidates < max(d$time)]
  tau <- candidates[which.max(candidates - pred$times[findInterval(candidates, pred$times)])]
  ordered <- d[order(d$time), ]
  sf <- survival::survfit(fit, newdata = ordered, se.fit = FALSE)
  expected <- 1 - as.numeric(summary(sf, times = tau, extend = TRUE)$surv)
  actual <- e$paper_simulation_risk_at(fit, ordered, tau)
  expect_equal(actual, expected, tolerance = 1e-12)
  old <- 1 - pred$surv_prob[, findInterval(tau, pred$times)]
  expect_gt(max(abs(actual - old)), 1e-6)
  cal <- e$paper_simulation_calibration(pred, tau, list(nk = 5), fit, d)
  reference <- pa_calibration(expected, pred$times, pred$status, tau)
  expect_equal(cal$ici, unname(reference$metrics[["ici"]]), tolerance = 1e-12)
  expect_equal(cal$diagnostics$prediction_grid_time, tau)
  expect_identical(cal$diagnostics$error, "")

  changed <- ordered
  changed$time <- rev(changed$time)
  changed$status <- 1 - changed$status
  expect_identical(e$paper_simulation_risk_at(fit, changed, tau), actual)
})

test_that("simulation AFT calibration requests the exact distribution probability", {
  e <- simulation_calibration_helpers()
  fit <- survival::survreg(survival::Surv(time, status) ~ x1 + x2,
    data = fx_surv(), dist = "lognormal", x = TRUE, y = TRUE)
  d <- fx_surv_uncensored()[seq(1, 200, by = 5), ]
  pred <- pa_predict_survreg(fit, fx_covs(), new_data = d)
  ordered <- d[order(d$time), ]
  tau <- median(d$time)
  expected <- stats::pnorm((log(tau) - predict(fit, newdata = ordered,
                                              type = "lp")) / fit$scale)
  expect_equal(e$paper_simulation_risk_at(fit, ordered, tau),
               as.numeric(expected), tolerance = 1e-12)
  cal <- e$paper_simulation_calibration(pred, tau, list(nk = 5), fit, d)
  expect_equal(cal$ici, unname(pa_calibration(as.numeric(expected),
    pred$times, pred$status, tau)$metrics[["ici"]]), tolerance = 1e-12)
  expect_equal(cal$diagnostics$prediction_grid_time, tau)
})

test_that("simulation calibration rejects misaligned outcomes and retains support failures", {
  e <- simulation_calibration_helpers()
  d <- fx_surv_uncensored()[seq(1, 200, by = 5), ]
  fit <- fx_cox()
  pred <- pa_predict_coxph(fit, fx_covs(), new_data = d)
  bad <- pred; bad$times <- rev(bad$times)
  expect_warning(result <- e$paper_simulation_calibration(
    bad, median(d$time), list(nk = 5), fit, d), "not aligned")
  expect_true(is.na(result$ici))
  expect_warning(result <- e$paper_simulation_calibration(
    pred, max(d$time) + 1, list(nk = 5), fit, d), "within observed follow-up")
  expect_true(is.na(result$ici))
})
