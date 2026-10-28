cr_event_fixture <- function() {
  list(time = c(1, 2, 3, 4, 6, 7), status = c(2, 0, 1, 1, 0, 1), tau = 5)
}

# Independent weighted-moment calculation, without fitting a regression.
cr_moment_reference <- function(y, p, w) {
  ybar <- sum(w * y)
  pbar <- sum(w * p)
  vy <- sum(w * (y - ybar)^2)
  vp <- sum(w * (p - pbar)^2)
  covariance <- sum(w * (y - ybar) * (p - pbar))
  corrected <- ybar + covariance / vp * (p - pbar)
  r <- round(covariance^2 / (vy * vp), 4)
  l <- round(sum(w * (y - corrected)^2) / sum(w * (y - p)^2), 4)
  c(r_square = r, l_square = l, pseudo_r2 = round(r * l, 4))
}

test_that("KM event weights use original times and exclude late censoring", {
  d <- cr_event_fixture()
  out <- PAmeasures:::pa_cr_event_ipcw(d$time, d$status, censoring = "km")
  # Censoring at 2: G falls to 4/5; censoring at 6: G falls to 2/5.
  # The competing event at 1 has already left the original risk set.
  expect_equal(out$probability, c(1, 1, .8, .8, .8, .4))
  expect_equal(out$weight, c(1, 0, 1.25, 1.25, 0, 2.5) / 6)
  expect_equal(out$observed_event, c(1, 0, 1, 1, 0, 1))
  expect_equal(out$weight[d$time == 6], 0)
})

test_that("late follow-up still affects event weights for a restricted outcome", {
  d <- cr_event_fixture()
  a <- PAmeasures:::pa_cr_event_ipcw(d$time, d$status, censoring = "km")
  changed <- d$time
  changed[5] <- 8 # Move a censoring time from before to after the last event.
  b <- PAmeasures:::pa_cr_event_ipcw(changed, d$status, censoring = "km")
  expect_false(isTRUE(all.equal(a$weight, b$weight)))
  expect_equal(pmin(d$time, d$tau), pmin(changed, d$tau))
})

test_that("interval and point CR metrics follow the event-weight moment identity", {
  d <- cr_event_fixture()
  grid <- c(1, 3, 5)
  cif <- outer(c(.1, .3, .5), seq(.3, .8, length.out = 6))
  w <- c(1, 0, 1.25, 1.25, 0, 2.5) / 6
  interval_prediction <- vapply(seq_len(ncol(cif)), function(i) {
    # The package's fixed step grid: [0,1), [1,3), [3,5).
    1 + 2 * (1 - cif[1, i]) + 2 * (1 - cif[2, i])
  }, numeric(1))
  outcome <- ifelse(d$status == 2, d$tau, pmin(d$time, d$tau))
  expected <- cr_moment_reference(outcome, interval_prediction, w)
  expected_point <- cr_moment_reference(
    as.numeric(d$time <= d$tau & d$status == 1), cif[3, ], w)
  res <- pa_survival_eval_cr(
    cif, d$time, grid, d$status, tau = d$tau, t_star = d$tau,
    metrics = c("r_square", "l_square", "pseudo_r2", "r2_point", "l2_point", "pseudo_r2_point"),
    censoring = "km")
  reference <- c(expected, pseudo_r2_point = unname(expected_point["pseudo_r2"]),
    r2_point = unname(expected_point["r_square"]),
    l2_point = unname(expected_point["l_square"]))
  expect_equal(res$Value, unname(reference[res$Metric]))
})

test_that("event weights and metrics do not depend on tied row order", {
  time <- c(1, 2, 2, 3, 6)
  status <- c(2, 0, 1, 1, 0)
  permutation <- c(1, 3, 2, 4, 5)
  for (method in c("cox", "km")) {
    a <- PAmeasures:::pa_cr_event_ipcw(time, status, censoring = method)
    b <- PAmeasures:::pa_cr_event_ipcw(
      time[permutation], status[permutation], censoring = method)
    expect_equal(a$weight[permutation], b$weight)
    expect_equal(a$probability[permutation], b$probability)
  }
})

test_that("conditional Cox probabilities match survfit and predictor shifts", {
  set.seed(4201)
  n <- 80
  x <- stats::rnorm(n, mean = 4)
  time_event <- stats::rexp(n, rate = .5)
  time_censor <- stats::rexp(n, rate = exp(.25 * x) / 4)
  time <- pmin(time_event, time_censor)
  status <- ifelse(time_event <= time_censor, rep(c(1, 2), length.out = n), 0)
  result <- PAmeasures:::pa_cr_event_ipcw(time, status, X = x)
  dat <- data.frame(time = time, censored = as.numeric(status == 0), x = x)
  fit <- survival::coxph(survival::Surv(time, censored) ~ x, data = dat)
  events <- which(status > 0)
  expected <- vapply(events, function(i) {
    sf <- survival::survfit(fit, newdata = dat[i, , drop = FALSE], se.fit = FALSE)
    as.numeric(summary(sf, times = time[i], extend = TRUE)$surv)
  }, numeric(1))
  expect_equal(result$probability[events], pmax(expected, 1e-4), tolerance = 1e-10)
  shifted <- PAmeasures:::pa_cr_event_ipcw(time, status, X = x + 100)
  expect_equal(result$probability, shifted$probability, tolerance = 1e-10)

  cif <- outer(c(.1, .3, .5), seq(.3, .8, length.out = n))
  grid <- c(.2, .7, 1.5)
  prediction <- list(cif_pred = cbind(grid, cif), times = time, status = status)
  direct <- pa_survival_eval_cr(cif, time, grid, status,
    metrics = c("r_square", "l_square", "pseudo_r2"), tau = 1.5, t_star = 1.5,
    censor_covariates = x)
  wide <- pa_summarize_cr(list(conditional = prediction),
    metrics = c("r_square", "l_square", "pseudo_r2"), tau = 1.5, t_star = 1.5,
    censor_covariates = x, digits = 4)
  expect_equal(wide$conditional, direct$Value[match(wide$Metric, direct$Metric)])
})

test_that("Cox and KM share the source probability floor and no-censor limit", {
  for (method in c("cox", "km")) {
    out <- PAmeasures:::pa_cr_event_ipcw(1:4, c(1, 2, 1, 2), censoring = method)
    expect_equal(out$probability, rep(1, 4))
    expect_equal(out$weight, rep(.25, 4))
    d <- cr_event_fixture()
    censored <- PAmeasures:::pa_cr_event_ipcw(d$time, d$status, censoring = method)
    expect_true(all(censored$probability >= 1e-4))
  }
  tail <- PAmeasures:::pa_cr_event_ipcw(1:10003,
    c(1, rep(0, 10001), 2), censoring = "km")
  expect_equal(tail$probability[10003], 1e-4)
  expect_equal(tail$weight[c(1, 10003)], c(1, 10000) / 10001)
})

test_that("summary forwards censoring choices and covariates without changing row alignment", {
  d <- cr_event_fixture()
  cif <- outer(c(.1, .3, .5), seq(.3, .8, length.out = 6))
  prediction <- list(cif_pred = cbind(c(1, 3, 5), cif), times = d$time, status = d$status)
  for (method in c("cox", "km")) {
    direct <- pa_survival_eval_cr(cif, d$time, c(1, 3, 5), d$status,
      metrics = c("pseudo_r2", "pseudo_r2_point"), tau = 5, t_star = 5,
      censoring = method)
    summary <- pa_summarize_cr(list(a = prediction, b = prediction),
      metrics = c("pseudo_r2", "pseudo_r2_point"), tau = 5, t_star = 5,
      censoring = method, digits = 4)
    expect_equal(summary$a, direct$Value[match(summary$Metric, direct$Metric)])
    expect_equal(summary$a, summary$b)
  }
  expect_error(PAmeasures:::pa_cr_event_ipcw(d$time, d$status, X = 1:3),
    "one row per subject")
  expect_error(PAmeasures:::pa_cr_event_ipcw(d$time, d$status, X = 1:6, censoring = "km"),
    "require censoring = 'cox'")
})
