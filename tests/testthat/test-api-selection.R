test_that("competing-risk AUC can be requested without Brier score", {
  d <- fx_cr()
  d$time <- d$obs.times
  d$status <- d$obs.event
  m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
                        data = d, x = TRUE, y = TRUE)
  m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
                        data = d, x = TRUE, y = TRUE)
  p <- pa_predict_cif(model1 = m1, model2 = m2, newdata = d,
                      covs = c("X1", "X2"))
  evaluate <- function(metrics) pa_survival_eval_cr(
    pred_cif = p$cif_pred[, -1, drop = FALSE], event_time = p$times,
    time.cif = p$cif_pred[, 1], status = p$status, metrics = metrics
  )
  all <- evaluate(NULL)
  auc <- evaluate("td_auc")
  expect_identical(auc$Metric, "td_auc")
  expect_equal(auc$Value, all$Value[all$Metric == "td_auc"])
})

test_that("fit-and-evaluate keeps metric values attached to their names", {
  all <- pa_fit_and_eval(fx_surv(), fx_covs(), metrics = "all")
  selected <- c("brier_score", "pseudo_r2", "harrell_c")
  subset <- pa_fit_and_eval(fx_surv(), fx_covs(), metrics = selected)
  expect_identical(names(subset), c("Model", selected))
  expect_equal(subset, all[, c("Model", selected)])
  expect_error(pa_fit_and_eval(fx_surv(), fx_covs(), metrics = "c_index"),
               "Invalid metrics")
  expect_error(pa_fit_and_eval(fx_surv(), fx_covs(), models = "unknown"),
               "Invalid models")
})

test_that("survreg prediction supports a single covariate in new data", {
  d <- fx_surv()
  fit <- survival::survreg(survival::Surv(time, status) ~ x1,
                           data = d, x = TRUE, y = TRUE, dist = "weibull")
  p <- pa_predict_survreg(fit, covs = "x1", new_data = d)
  expected <- stats::predict(fit, newdata = d, type = "lp")[order(d$time)]
  expect_equal(unname(p$linear.pred), unname(expected))
  expect_identical(dim(p$surv_prob), c(nrow(d), nrow(d)))
})

test_that("Cox validation probabilities depend on the fitted model, not labels", {
  train <- fx_surv()
  fit <- fx_cox()
  d <- fx_surv_uncensored()[1:80, ]
  p <- pa_predict_coxph(fit, fx_covs(), new_data = d)
  d$status[seq(1, nrow(d), by = 2)] <- 0
  changed <- pa_predict_coxph(fit, fx_covs(), new_data = d)
  expect_equal(changed$surv_prob, p$surv_prob)
  ref <- survival::survfit(fit, newdata = d[order(d$time), ], se.fit = FALSE)
  expected <- t(summary(ref, times = sort(d$time), extend = TRUE)$surv)
  expect_equal(unname(p$surv_prob), unname(expected), tolerance = 1e-12)
})

test_that("finite-horizon AFT means integrate survival and preserve outcomes", {
  d <- fx_surv()
  tau <- unname(stats::median(d$time))
  for (distribution in c("weibull", "lognormal")) {
    fit <- survival::survreg(survival::Surv(time, status) ~ x1 + x2,
                             data = d, x = TRUE, y = TRUE, dist = distribution)
    p <- pa_predict_survreg(fit, fx_covs(), new_data = d, tau = tau)
    eta <- stats::predict(fit, d[order(d$time), ], type = "lp")
    expected <- vapply(eta, function(mu) stats::integrate(function(t) {
      if (distribution == "weibull") exp(-(t / exp(mu))^(1 / fit$scale))
      else stats::pnorm((log(t) - mu) / fit$scale, lower.tail = FALSE)
    }, 0, tau)$value, numeric(1))
    expect_equal(unname(p$pred), unname(expected), tolerance = 1e-6)
    expect_true(all(p$pred >= 0 & p$pred <= tau))
    expect_equal(p$times, sort(d$time))
    expect_equal(p$status, d$status[order(d$time)])
  }
})

test_that("explained variation keeps predictions and weights aligned when sorting", {
  y <- c(1, 2, 4, 8, 10)
  p <- c(1.2, 2.3, 3.8, 7.5, 9.7)
  w <- c(1, 2, 1, 3, 2)
  i <- c(5, 3, 1, 4, 2)
  evaluate <- function(idx) PAmeasures:::pam.r2_metrics(
    p[idx], y[idx], rep(1, length(idx)), case_weight = w[idx])
  expect_equal(evaluate(i), evaluate(seq_along(y)))
})

test_that("NCC weights use the documented number of controls and matching strata", {
  expect_equal(pa_nested_case_control_weights(c(1, 2, 3), c(1, 0, 0), m = 1),
               c(1, 2, 2))
  expect_equal(pa_nested_case_control_weights(c(1, 2, 3), c(1, 0, 0), m = 2),
               c(1, 1, 1))
  # Events in another stratum must not add a selection opportunity.
  time <- c(1, 2, 3, 4, 5, 6)
  status <- c(1, 1, 0, 0, 0, 0)
  strata <- c("a", "b", "a", "a", "b", "b")
  expect_equal(pa_nested_case_control_weights(time, status, strata, m = 1),
               c(1, 1, 2, 2, 2, 2))
  expect_equal(pa_nested_case_control_weights(c(1, 1, 2), c(1, 0, 0), m = 1),
               c(1, 2, 2))
  # Two tied cases each draw one of two controls independently: P(selected)=3/4.
  expect_equal(pa_nested_case_control_weights(c(1, 1, 2, 3), c(1, 1, 0, 0), m = 1),
               c(1, 1, 4/3, 4/3))
  expect_error(pa_nested_case_control_weights(c(1, 2), c(1, 0), m = 0),
               "positive integer")
})

test_that("two-phase metric selection accepts all and rejects unknown names", {
  d <- fx_surv()
  p <- pa_predict_coxph(fx_cox(), fx_covs(), new_data = d)
  km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, d)
  evaluate <- function(metrics) pa_evaluate_two_phase(
    p, km_cens_fit = km, case_weights = rep(1, nrow(d)), metrics = metrics)
  expect_error(evaluate("c_index"), "Invalid metrics")
  expect_setequal(evaluate("all")$Metric,
                   c("pseudo_r2", "r_square", "l_square", "harrell_c",
                     "uno_c", "brier_score", "td_auc", "ici", "e50", "e90"))
  expect_equal(nrow(evaluate(NULL)), 6L)
})

test_that("point pseudo R-squared compares events with event probabilities", {
  d <- fx_surv_uncensored()
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         d, x = TRUE, y = TRUE)
  p <- pa_predict_coxph(fit, fx_covs(), new_data = d)
  # Use an observed event time to check inclusion of events at the horizon.
  j <- 100L
  t_star <- p$times[j]
  y <- as.numeric(p$times <= t_star)
  risk <- 1 - p$surv_prob[, j]
  calibration <- stats::lm(y ~ risk)
  r2 <- round(summary(calibration)$r.squared, 4)
  l2 <- round(sum(stats::residuals(calibration)^2) / sum((y - risk)^2), 4)
  result <- pa_survival_eval(fit, p$times, p$surv_prob, status = p$status,
    covariates = fx_covs(), metrics = "pseudo_r2_point", t_star = t_star)
  expect_equal(result$Value, round(r2 * l2, 4))
})

test_that("held-out AUC does not hide reversed discrimination", {
  d <- pa_sim_cox_weibull(120, 0, 2, c(1, -1), seed = 9191)
  d <- d[, c("time", "status", "x1", "x2")]
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         d, x = TRUE, y = TRUE)
  validation <- d
  validation$x1 <- -validation$x1
  validation$x2 <- -validation$x2
  t <- stats::median(validation$time)
  score <- stats::predict(fit, validation, type = "lp")
  cases <- validation$time <= t
  # Uncensored empirical AUC: probability a case has higher risk than a control.
  comparisons <- outer(score[cases], score[!cases], "-")
  reference <- mean((comparisons > 0) + 0.5 * (comparisons == 0))
  result <- pa_fit_and_eval(d, fx_covs(), predicted_data = validation,
                            metrics = "td_auc", t_star = t)
  expect_lt(result$td_auc, 0.5)
  expect_equal(result$td_auc, round(reference, 2))
})

tail_support_data <- function() {
  set.seed(99)
  train <- data.frame(time = seq_len(40),
    status = rep(c(1, 1, 0, 1), 10), x1 = stats::rnorm(40), x2 = stats::rnorm(40))
  train$status[40] <- 0
  evaluation <- train[1:12, ]
  evaluation$time[12] <- 45
  evaluation$status[12] <- 1
  list(train = train, evaluation = evaluation)
}

test_that("metric subsets avoid unused R-squared outside censoring support", {
  d <- tail_support_data()
  result <- expect_no_warning(pa_fit_and_eval(d$train, c("x1", "x2"),
    predicted_data = d$evaluation, models = "all",
    metrics = c("brier_score", "harrell_c"), t_star = 5))
  expect_identical(names(result), c("Model", "brier_score", "harrell_c"))
  expect_true(all(is.finite(as.matrix(result[, -1]))))
})

test_that("undefined censoring weights return missing R-squared with a diagnostic", {
  d <- tail_support_data()
  formula <- survival::Surv(time, status) ~ x1 + x2
  for (kind in c("coxph", "exponential", "weibull", "lognormal")) {
    fit <- if (kind == "coxph") {
      survival::coxph(formula, data = d$train, x = TRUE, y = TRUE)
    } else {
      survival::survreg(formula, data = d$train, dist = kind, x = TRUE, y = TRUE)
    }
    fun <- if (kind == "coxph") pa_predict_coxph else pa_predict_survreg
    expect_warning(result <- fun(fit, c("x1", "x2"), new_data = d$evaluation,
      predict = FALSE), "require positive censoring survival")
    expect_true(all(is.na(unlist(result))))
    prediction <- expect_no_warning(fun(fit, c("x1", "x2"),
      new_data = d$evaluation, predict = TRUE))
    expect_true(all(is.finite(prediction$surv_prob)))
  }
})
