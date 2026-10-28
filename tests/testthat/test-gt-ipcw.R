# Regression tests for finding 41: Gt() must be a correct Kaplan-Meier estimate
# of the censoring distribution, evaluated as a STEP function.
#
# Contract, per the standard IPCW construction (Graf et al. 1999):
#   Gt(object, t)              = G-hat(t),   right-continuous, for the
#                                evaluation-time weight
#   Gt(object, t, left = TRUE) = G-hat(t-),  the left limit, for
#                                subject-specific event-time weights
#
# No linear interpolation between jump times. Values are checked against
# pec::ipcw(), an independent implementation, rather than against values this
# package produced.
#
# These tests were written to FAIL against the pre-correction implementation.

skip_if_no_pec <- function() testthat::skip_if_not_installed("pec")

# fixture: ties, and an EVEN number of distinct event times, so the default
# t_star = median(distinct event times) is NOT an observed time and the
# evaluation genuinely exercises a non-observed timepoint.
tie_fixture <- function() {
  set.seed(4141)
  n <- 200
  d <- data.frame(
    time   = rep(1:10, each = 20),
    status = rep(c(1L, 1L, 1L, 0L), length.out = n),
    x1     = stats::rnorm(n),
    x2     = stats::rnorm(n)
  )
  d
}

cens_km <- function(d) {
  survival::survfit(survival::Surv(d$time, as.numeric(d$status == 0)) ~ 1)
}

test_that("the tie fixture really has an even number of distinct event times", {
  d <- tie_fixture()
  ev <- sort(unique(d$time[d$status == 1]))
  expect_true(length(ev) %% 2 == 0)
  t_star <- stats::median(ev)
  expect_false(t_star %in% d$time)   # median of an even count falls between times
})

test_that("Gt() is invariant to row permutation", {
  d <- tie_fixture()
  set.seed(808)
  dp <- d[sample(nrow(d)), ]
  sv  <- survival::Surv(d$time,  d$status)
  svp <- survival::Surv(dp$time, dp$status)

  ts <- c(1.5, 2.5, 4.5, 5.5, 7.5, 9.5, 3, 7)
  for (t in ts) {
    expect_identical(PAmeasures:::Gt(sv, t), PAmeasures:::Gt(svp, t),
                     info = paste("G(t) at t =", t))
    expect_identical(PAmeasures:::Gt(sv, t, left = TRUE),
                     PAmeasures:::Gt(svp, t, left = TRUE),
                     info = paste("G(t-) at t =", t))
  }
})

test_that("the Brier score is invariant to row permutation", {
  d <- tie_fixture()
  set.seed(808)
  dp <- d[sample(nrow(d)), ]
  f  <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = d,  x = TRUE, y = TRUE)
  fp <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = dp, x = TRUE, y = TRUE)

  expect_equal(as.numeric(PAmeasures:::pam.Brier(f,  d,  5.5)),
               as.numeric(PAmeasures:::pam.Brier(fp, dp, 5.5)),
               tolerance = 1e-12)
})

test_that("censoring survival values are in [0, 1] and non-increasing", {
  d <- tie_fixture()
  sv <- survival::Surv(d$time, d$status)
  # Restricted to where G is defined and positive. On this fixture the reverse
  # KM reaches zero at the last time (t = 10), where Gt() errors by design, and
  # is NA beyond it -- both are separately tested below.
  grid <- seq(0.5, 9.75, by = 0.25)

  g  <- vapply(grid, function(t) PAmeasures:::Gt(sv, t), numeric(1))
  gl <- vapply(grid, function(t) PAmeasures:::Gt(sv, t, left = TRUE), numeric(1))

  expect_true(all(g  >= 0 & g  <= 1), info = "G(t) outside [0, 1]")
  expect_true(all(gl >= 0 & gl <= 1), info = "G(t-) outside [0, 1]")
  expect_true(all(diff(g)  <= 1e-12), info = "G(t) is not non-increasing")
  expect_true(all(diff(gl) <= 1e-12), info = "G(t-) is not non-increasing")
  # the left limit is never below the right-continuous value at the same t
  expect_true(all(gl >= g - 1e-12))
})

test_that("with no censoring, G(t) is 1 across the supported range", {
  d <- data.frame(time = 1:10, status = rep(1L, 10))
  sv <- survival::Surv(d$time, d$status)
  for (t in c(1, 2, 5, 8, 9.5, 10)) {
    expect_equal(PAmeasures:::Gt(sv, t), 1,
                 info = paste("G(t) at t =", t, "with no censoring"))
  }
})

test_that("Gt() agrees with pec::ipcw IPCW.times for the evaluation-time weight", {
  skip_if_no_pec()
  d <- tie_fixture()
  sv <- survival::Surv(d$time, d$status)
  ts <- c(1, 2.5, 4, 5.5, 7, 8.5, 9)

  ref <- pec::ipcw(survival::Surv(time, status) ~ 1, data = d, method = "marginal",
                   times = ts, subjectTimes = d$time, subjectTimesLag = 1,
                   what = "IPCW.times")$IPCW.times
  obs <- vapply(ts, function(t) PAmeasures:::Gt(sv, t), numeric(1))

  expect_equal(obs, as.numeric(ref), tolerance = 1e-12)
})

test_that("Gt(left = TRUE) agrees with pec::ipcw IPCW.subjectTimes", {
  skip_if_no_pec()
  d <- tie_fixture()
  sv <- survival::Surv(d$time, d$status)

  ref <- pec::ipcw(survival::Surv(time, status) ~ 1, data = d, method = "marginal",
                   times = sort(unique(d$time)), subjectTimes = d$time,
                   subjectTimesLag = 1, what = "IPCW.subjectTimes")$IPCW.subjectTimes
  obs <- vapply(d$time, function(t) PAmeasures:::Gt(sv, t, left = TRUE), numeric(1))

  expect_equal(obs, as.numeric(ref), tolerance = 1e-12)
})

test_that("Gt() errors rather than substituting when censoring survival is zero", {
  # largest observation censored => the censoring KM drops to exactly 0 there.
  # The old code silently returned min(surv[surv != 0]); it must not.
  d <- data.frame(time = 1:10, status = c(rep(1L, 9), 0L))
  sv <- survival::Surv(d$time, d$status)

  expect_error(PAmeasures:::Gt(sv, 10), "censoring distribution")
  expect_error(PAmeasures:::pam.Brier(
    survival::coxph(survival::Surv(time, status) ~ 1, data = d, x = TRUE, y = TRUE),
    d, 10), "censoring distribution")
})

test_that("Gt() is NA beyond the last observed time, as pec::ipcw is", {
  skip_if_no_pec()
  d <- tie_fixture()
  sv <- survival::Surv(d$time, d$status)
  beyond <- max(d$time) + 0.5

  ref <- pec::ipcw(survival::Surv(time, status) ~ 1, data = d, method = "marginal",
                   times = beyond, subjectTimes = d$time, subjectTimesLag = 1,
                   what = "IPCW.times")$IPCW.times
  expect_true(is.na(ref))
  expect_true(is.na(PAmeasures:::Gt(sv, beyond)))
})

test_that("pa_fit_and_eval gives the corrected Brier score on the exposed fixture", {
  d <- tie_fixture()
  res <- suppressMessages(pa_fit_and_eval(
    train_data = d, covariates = c("x1", "x2"),
    models = "coxph", metrics = "brier_score", t_star = 5.5
  ))
  b <- as.numeric(res$brier_score)

  expect_true(is.finite(b))
  expect_gte(b, 0)
  expect_lte(b, 1)
  # the pre-correction implementation inflated this well above the correct
  # value on tied data at a non-observed t_star
  expect_lt(b, 0.35)
})
