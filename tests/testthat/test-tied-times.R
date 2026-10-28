# Tied follow-up times: the IPCW left limit G(t-) must be shared by every row
# tied at t, so metrics cannot depend on the order of tied rows.

tie_surv_data <- function(n = 150) {
  set.seed(7)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$time <- ceiling(stats::rexp(n, exp(.5 * d$x1)) * 4)
  d$status <- stats::rbinom(n, 1, .6)
  d
}

tie_surv_eval <- function(fit, d, metrics) {
  p <- pa_predict_coxph(fit, c("x1", "x2"), new_data = d)
  r <- pa_survival_eval(model = fit, event_time = p$times,
    predicted_probability = p$surv_prob, status = p$status,
    covariates = c("x1", "x2"), new_data = d, tau = 10e10, metrics = metrics)
  stats::setNames(r$Value, r$Metric)
}

test_that("censoring left limit is shared within tied times", {
  time <- c(3, 1, 2, 2, 2, 4, 2)
  status <- c(1, 0, 0, 1, 0, 1, 1)
  fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  # Censoring KM: 6/7 after t = 1, then 6/7 * 4/6 after t = 2 (2 of 6 at risk).
  expect_equal(PAmeasures:::pa_censoring_left_limit(fit, time),
    c(4 / 7, 1, 6 / 7, 6 / 7, 6 / 7, 4 / 7, 6 / 7))
})

test_that("right-censored metrics do not depend on the order of tied rows", {
  d <- tie_surv_data()
  expect_true(any(duplicated(d$time[d$status == 0])))
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = d,
    x = TRUE, y = TRUE)
  metrics <- c("pseudo_r2", "r_square", "l_square", "r2_point", "l2_point")
  base <- tie_surv_eval(fit, d, metrics)
  for (o in list(order(d$time, d$status), order(d$time, -d$status),
                 rev(seq_len(nrow(d))))) {
    expect_equal(tie_surv_eval(fit, d[o, ], metrics), base)
  }
})

test_that("point R2 matches the KM reference when censored at t*", {
  set.seed(7)
  n <- 101
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$time <- stats::rexp(n, exp(.5 * d$x1))
  d$status <- stats::rbinom(n, 1, .7)
  d$status[order(d$time)[51]] <- 0
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = d,
    x = TRUE, y = TRUE)
  p <- pa_predict_coxph(fit, c("x1", "x2"), new_data = d)
  t_star <- stats::median(p$times)
  expect_true(any(p$times == t_star & p$status == 0))

  # Restricted outcomes; KM left limit from censor counts and risk sets.
  time <- pmin(p$times, t_star)
  observed <- as.numeric(p$status == 1 | p$times > t_star)
  censor_times <- sort(unique(time[observed == 0]))
  factors <- vapply(censor_times, function(t)
    1 - sum(time == t & observed == 0) / sum(time >= t), numeric(1))
  G <- vapply(time, function(t) prod(factors[censor_times < t]), numeric(1))
  w <- (observed / G) / sum(observed / G)
  k <- PAmeasures:::pa_nearest_time_index(p$times, t_star)
  prediction <- 1 - p$surv_prob[, k]
  y <- as.numeric(p$times <= t_star & p$status == 1)
  ybar <- sum(w * y); pbar <- sum(w * prediction)
  vp <- sum(w * (prediction - pbar)^2)
  cv <- sum(w * (y - ybar) * (prediction - pbar))
  corrected <- ybar + cv / vp * (prediction - pbar)
  reference <- c(
    r2_point = round(cv^2 / (sum(w * (y - ybar)^2) * vp), 4),
    l2_point = round(sum(w * (y - corrected)^2) / sum(w * (y - prediction)^2), 4))

  r <- pa_survival_eval(model = fit, event_time = p$times,
    predicted_probability = p$surv_prob, status = p$status,
    covariates = c("x1", "x2"), new_data = d, tau = 10e10,
    metrics = c("r2_point", "l2_point"))
  expect_equal(stats::setNames(r$Value, r$Metric)[names(reference)], reference)
})

test_that("null Cox censoring equals exp(-Nelson-Aalen) with tied censoring", {
  time <- c(1, 2, 2, 2, 3, 4, 5, 6)
  status <- c(1, 0, 0, 1, 0, 1, 2, 1)
  ct <- sort(unique(time[status == 0]))
  jumps <- vapply(ct, function(t) sum(time == t & status == 0) / sum(time >= t),
    numeric(1))
  G <- vapply(time, function(t) exp(-sum(jumps[ct < t])), numeric(1))
  out <- PAmeasures:::pa_cr_event_ipcw(time, status, censoring = "cox")
  expect_equal(out$probability, G)
})

test_that("competing-risk metrics do not depend on the order of tied rows", {
  set.seed(11)
  n <- 120
  x <- stats::rnorm(n)
  time <- ceiling(stats::rexp(n, exp(.4 * x)) * 3)
  status <- sample(0:2, n, TRUE, prob = c(.35, .4, .25))
  expect_true(any(duplicated(time[status == 0])))
  grid <- sort(unique(time))
  cif <- outer(grid, x, function(t, xi) .6 * (1 - exp(-.3 * t * exp(.4 * xi))))
  ev <- function(o, censoring) {
    r <- pa_survival_eval_cr(pred_cif = cif[, o], event_time = time[o],
      time.cif = grid, status = status[o], event_type = 1,
      censoring = censoring, metrics = "all")
    stats::setNames(r$Value, r$Metric)
  }
  for (censoring in c("cox", "km")) {
    base <- ev(seq_len(n), censoring)
    expect_equal(ev(order(time, status), censoring), base)
    expect_equal(ev(order(time, -status), censoring), base)
  }
})

test_that("Cox prediction supports strata-by-covariate interactions", {
  set.seed(1)
  n <- 120
  d <- data.frame(x1 = stats::rnorm(n),
                  g = factor(sample(c("A", "B"), n, TRUE)))
  d$time <- stats::rexp(n, exp(.5 * d$x1))
  d$status <- stats::rbinom(n, 1, .7)
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 * strata(g),
    data = d, x = TRUE, y = TRUE)
  nd <- d[1:8, ]
  p <- expect_no_warning(pa_predict_coxph(fit, c("x1", "g"), new_data = nd))
  for (i in seq_len(nrow(nd))) {
    j <- which(p$times == nd$time[i] & p$status == nd$status[i])[1]
    sf <- survival::survfit(fit, newdata = nd[i, ])
    expected <- summary(sf, times = p$times, extend = TRUE)$surv
    expect_equal(p$surv_prob[j, ], expected)
  }
  expect_error(pa_predict_coxph(fit, c("x1", "g")), "pass the evaluation data")

  fit_int <- survival::coxph(survival::Surv(time, status) ~ x1 * g, data = d,
    x = TRUE, y = TRUE)
  expect_no_warning(pa_predict_coxph(fit_int, c("x1", "g"), new_data = nd))
})
