calibration_data <- function() {
  withr::local_seed(29)
  n <- 180
  x <- stats::rnorm(n)
  a <- .12 * exp(.5 * x)
  b <- .08 * exp(-.3 * x)
  target <- stats::rexp(n, a)
  competing <- stats::rexp(n, b)
  censor <- stats::rexp(n, .06)
  list(time = pmin(target, censor), status = as.integer(target <= censor),
       cr_time = pmin(target, competing, censor),
       cr_status = ifelse(censor < pmin(target, competing), 0L,
                           ifelse(target < competing, 1L, 2L)),
       rate = a, cr_rate = b)
}

test_that("right-censored calibration matches the selected public backend", {
  d <- calibration_data()
  p <- 1 - exp(-d$rate * 3)
  reference <- pmcalibration::pmcalibration(
    y = survival::Surv(d$time, d$status), p = p, time = 3,
    smooth = "rcs", nk = 5, ci = "none", plot = FALSE)
  result <- pa_calibration(p, d$time, d$status, 3)
  expect_equal(unname(result$metrics), unname(reference$metrics[c("Eavg", "E50", "E90")]))
  expect_identical(names(result$metrics), c("ici", "e50", "e90"))
  expect_equal(result$curve$predicted, p)
  expect_equal(result$curve$calibrated, reference$p_c$p_c)
  expect_equal(result$t_star, 3)
  expect_true(all(is.finite(result$metrics)))
  expect_true(all(result$metrics >= 0 & result$metrics <= 1))
})

test_that("constant right-censored risks use an independent KM reference", {
  time <- 1:6
  status <- c(1, 0, 1, 0, 1, 0)
  # KM survival at 4 = (1 - 1/6) * (1 - 1/4) = 5/8.
  result <- pa_calibration(rep(.2, 6), time, status, 4)
  expect_equal(result$curve$calibrated, rep(3/8, 6))
  expect_equal(unname(result$metrics), rep(3/8 - .2, 3))
  expect_equal(result$method, "kaplan_meier_constant")
})

test_that("competing-risk metrics average over subjects, not distinct risks", {
  # No censoring before 5: target incidences are exactly 2/10 and 3/5.
  # Unequal group sizes and tied probabilities detect averaging over a grid
  # or unique prediction values instead of the original subjects.
  time <- c(rep(1, 2), rep(2, 3), rep(8, 5), rep(1, 3), 2, 8)
  status <- c(rep(7L, 2), rep(9L, 3), rep(0L, 5), rep(7L, 3), 9L, 0L)
  p <- c(rep(.1, 10), rep(.8, 5))
  result <- pa_calibration_cr(p, time, status, 5, event_type = 7)
  expect_equal(result$curve$calibrated, c(rep(.2, 10), rep(.6, 5)))
  expect_equal(result$metrics, c(ici = (10*.1 + 5*.2)/15, e50 = .1, e90 = .2))

  # More than two causes are pooled correctly; changing their labels or
  # subject order cannot change a target-cause calibration assessment.
  status[c(3, 4)] <- 12L
  o <- rev(seq_along(p))
  other <- pa_calibration_cr(p[o], time[o], status[o], 5, event_type = 7)
  expect_equal(other$metrics, result$metrics)
  expect_equal(other$curve$calibrated, rev(result$curve$calibrated))
  boundary <- pa_calibration_cr(rep(c(0, 1), c(10, 5)), time, status, 5,
                                event_type = 7)
  expect_equal(boundary$metrics, c(ici = 4/15, e50 = .2, e90 = .4))
})

test_that("competing-risk calibration includes target events at t_star", {
  p <- rep(.25, 6)
  time <- c(1, 2, 3, 5, 6, 7)
  status <- c(1, 2, 2, 1, 0, 0)
  before <- pa_calibration_cr(p, time, status, 4)
  at <- pa_calibration_cr(p, time, status, 5)
  expect_equal(before$curve$calibrated, rep(1/6, 6))
  expect_equal(at$curve$calibrated, rep(2/6, 6))
})

test_that("boundary probabilities are retained in right-censored summaries", {
  d <- calibration_data()
  p <- 1 - exp(-d$rate * 3)
  p[1:2] <- c(0, 1)
  expect_warning(result <- pa_calibration(p, d$time, d$status, 3),
                  "Boundary probabilities")
  expect_equal(nrow(result$curve), length(p))
  expect_equal(result$curve$predicted, p)
  error <- abs(p - result$curve$calibrated)
  expect_equal(result$metrics, c(ici = mean(error), e50 = stats::median(error),
                                 e90 = unname(stats::quantile(error, .9))))
})

test_that("calibration input failures have explicit diagnostics", {
  d <- calibration_data()
  p <- 1 - exp(-d$rate * 3)
  expect_error(pa_calibration(p[-1], d$time, d$status, 3), "one per subject")
  expect_error(pa_calibration(replace(p, 1, NA), d$time, d$status, 3), "probability")
  expect_error(pa_calibration(p, d$time, d$status, max(d$time)+1), "follow-up")
  expect_error(pa_calibration(p, d$time, d$cr_status, 3), "only 0")
  expect_error(pa_calibration(p, d$time, d$status, 3, nk = 2), "at least 3")
  expect_error(pa_calibration_cr(p, d$cr_time, d$cr_status, 3, event_type = 0),
               "positive integer")
  expect_error(pa_calibration_cr(p, d$cr_time, d$cr_status, 3, bandwidth = -1),
               "bandwidth")
  expect_error(pa_calibration_cr(p, d$cr_time, rep(0, length(p)), 3),
               "observed target event")
})

test_that("matrix evaluators use past grid values at the requested horizon", {
  d <- calibration_data()
  grid <- sort(d$time)
  k <- length(grid) %/% 2
  t <- grid[k] + .8*(grid[k+1]-grid[k]) # nearer the future grid column
  survival <- outer(d$rate, grid, function(a, t) exp(-a*t))
  p <- 1 - survival[, k]
  reference <- pa_calibration(p, d$time, d$status, t, nk = 3)$metrics
  metrics <- c("ici", "e50", "e90")
  result <- pa_survival_eval(event_time = d$time, status = d$status,
    predicted_probability = survival, metrics = metrics, t_star = t,
    calibration_args = list(nk = 3))
  expect_equal(stats::setNames(result$Value, result$Metric), round(reference, 4))

  cr_grid <- c(1, 3, 5)
  cif <- outer(cr_grid, seq_along(d$rate), function(t, j) {
    a <- d$rate[j]; b <- d$cr_rate[j]
    a/(a+b) * (1-exp(-(a+b)*t))
  })
  cr <- pa_survival_eval_cr(cif, d$cr_time, cr_grid, d$cr_status,
    metrics = metrics, t_star = 4.8, calibration_args = list(bandwidth = .2))
  reference <- pa_calibration_cr(cif[2, ], d$cr_time, d$cr_status, 4.8,
                                  bandwidth = .2)$metrics
  expect_equal(stats::setNames(cr$Value, cr$Metric), round(reference, 4))
  o <- rev(seq_along(d$time))
  permuted <- pa_survival_eval_cr(cif[, o], d$cr_time[o], cr_grid, d$cr_status[o],
    metrics = metrics, t_star = 4.8, calibration_args = list(bandwidth = .2))
  expect_equal(permuted, cr)

  # A 1-row CIF grid remains a matrix; no integrated time prediction is needed.
  one <- pa_survival_eval_cr(cif[2, , drop = FALSE], d$cr_time, 3, d$cr_status,
                             metrics = metrics, t_star = 3)
  expect_true(all(is.finite(one$Value)))
  expect_error(pa_survival_eval_cr(cif, d$cr_time, cr_grid, d$cr_status,
                                  metrics = "ici", t_star = 6), "prediction grid")
})

test_that("default comparison tables add only ICI from the calibration family", {
  p <- pa_predict_coxph(fx_cox(), fx_covs(), new_data = fx_surv())
  tab <- pa_summarize(list(cox = p))
  expect_true("ici" %in% tab$Metric)
  expect_false(any(c("e50", "e90") %in% tab$Metric))
  all <- pa_summarize(list(cox = p), metrics = "all")
  expect_true(all(c("ici", "e50", "e90") %in% all$Metric))
  selected <- pa_summarize(list(cox = p), metrics = c("ici", "e90"), digits = 4)
  expect_identical(selected$Metric, c("ici", "e90"))
  expect_false(anyNA(selected$Metric))

  # Use externally supplied CIFs to exercise the model-free summary path.
  d <- calibration_data()
  grid <- sort(unique(d$cr_time))
  cif <- outer(grid, seq_along(d$rate), function(t, j) {
    a <- d$rate[j]; b <- d$cr_rate[j]
    a/(a+b) * (1-exp(-(a+b)*t))
  })
  model <- list(cif_pred = cbind(time = grid, cif), times = d$cr_time,
                status = d$cr_status)
  tab <- pa_summarize_cr(list(external = model))
  expect_true("ici" %in% tab$Metric)
  expect_false(any(c("e50", "e90") %in% tab$Metric))
  all <- pa_summarize_cr(list(external = model), metrics = "all")
  expect_true(all(c("ici", "e50", "e90") %in% all$Metric))
})

test_that("fit-and-evaluate selects calibration and keeps ICI alone by default", {
  d <- fx_surv()
  t <- stats::median(d$time)
  default <- pa_fit_and_eval(d, fx_covs())
  expect_true("ici" %in% names(default))
  expect_false(any(c("e50", "e90") %in% names(default)))
  selected <- pa_fit_and_eval(d, fx_covs(), models = c("coxph", "weibull"),
    metrics = c("e90", "ici", "e50"), t_star = t)
  expect_identical(names(selected), c("Model", "e90", "ici", "e50"))
  expect_true(all(is.finite(as.matrix(selected[-1]))))
  for (i in seq_len(nrow(selected))) {
    fit <- if (selected$Model[i] == "coxph") fx_cox() else fx_survreg()
    p <- if (selected$Model[i] == "coxph") {
      1-as.numeric(pec::predictSurvProb(fit, newdata = d, times = t))
    } else {
      # Independent Weibull event probability from the AFT parameterization.
      eta <- as.numeric(stats::predict(fit, d, type = "lp"))
      1-exp(-(t/exp(eta))^(1/fit$scale))
    }
    reference <- pa_calibration(p, d$time, d$status, t)$metrics
    expect_equal(unlist(selected[i, -1], use.names = FALSE),
                  round(unname(reference[c("e90", "ici", "e50")]), 2))
  }
})


test_that("three and four distinct risks use a quiet linear Cox calibration", {
  withr::local_seed(882)
  n <- 90L
  time <- rexp(n)
  status <- rbinom(n, 1, .7)
  for (distinct in 3:4) {
    p <- rep(seq(.1, .7, length.out = distinct), length.out = n)
    expect_silent(result <- pa_calibration(p, time, status, .7))
    expect_identical(result$method, "pmcalibration_none")
    expect_true(all(is.finite(result$metrics)))
    # Independent Cox fit and survfit prediction, rather than the backend's
    # own calibration summary, check the fallback's numerical meaning.
    z <- log(-log1p(-p))
    fit <- survival::coxph(survival::Surv(time, status) ~ z)
    curves <- survival::survfit(fit, newdata = data.frame(z = z))
    reference <- 1 - as.numeric(summary(curves, times = .7)$surv)
    expect_equal(result$curve$calibrated, reference, tolerance = 1e-8)
    error <- abs(p - reference)
    expect_equal(result$metrics, c(ici = mean(error), e50 = median(error),
      e90 = unname(quantile(error, .9))), tolerance = 1e-8)
    expect_silent(weighted <- pa_calibration(p, time, status, .7,
      case_weights = rep(1, n)))
    expect_equal(result$curve, weighted$curve, tolerance = 1e-8)
    expect_equal(result$metrics, weighted$metrics, tolerance = 1e-8)
  }
})
