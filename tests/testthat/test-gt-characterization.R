# Characterization of Gt(), the Kaplan-Meier estimate of the censoring
# distribution.
#
# SUPERSEDED AND INTENTIONALLY UPDATED, 2026-09-25.
#
# These expectations originally pinned the pre-correction implementation, to
# prove that removing the survminer dependency was behaviour-preserving. That
# purpose is served and finished. Finding 41 then established that the pinned
# behaviour was itself wrong in three independent ways:
#
#   1. the interpolation branch indexed the raw, unsorted input vector with
#      indices into the sorted summary table, so G(t) was neither monotone nor
#      invariant to row order;
#   2. the censoring indicator was derived as status == min(status), which
#      inverted the coding entirely on data with no censoring, making G decay
#      from 1 instead of staying at 1;
#   3. a zero censoring survival was silently replaced by the smallest positive
#      value, hiding an undefined IPCW weight behind a plausible number.
#
# The values below are therefore the CORRECTED ones, from the reverse-
# Kaplan-Meier step-function estimator. They were not regenerated to make a
# failing test pass: the correction was reviewed and approved, and old and new
# values were reported side by side first.
#
# The authoritative contract now lives in test-gt-ipcw.R, which checks against
# pec::ipcw() rather than against values this package produced. This file keeps
# the boundary, tie and validation coverage.

gt_ <- function(...) PAmeasures:::Gt(...)

test_that("Gt is pinned on ordinary right-censored data", {
  d  <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  qs <- unname(stats::quantile(d$time, c(0.1, 0.25, 0.5, 0.75, 0.9)))

  expect_equal(gt_(sv, qs[1]), 0.994897959184, tolerance = 1e-10)
  expect_equal(gt_(sv, qs[2]), 0.948234671222, tolerance = 1e-10)
  expect_equal(gt_(sv, qs[3]), 0.872896816372, tolerance = 1e-10)
  expect_equal(gt_(sv, qs[4]), 0.695781808710, tolerance = 1e-10)
  expect_equal(gt_(sv, qs[5]), 0.498466785981, tolerance = 1e-10)
})

test_that("Gt is pinned at timepoints falling exactly on an observed time", {
  d   <- fx_surv()
  sv  <- survival::Surv(d$time, d$status)
  f   <- survival::survfit(survival::Surv(d$time, as.numeric(d$status == 0)) ~ 1)
  on  <- f$time[c(5, 50, 150)]

  expect_equal(gt_(sv, on[1]), 0.994897959184, tolerance = 1e-10)
  expect_equal(gt_(sv, on[2]), 0.948234671222, tolerance = 1e-10)
  expect_equal(gt_(sv, on[3]), 0.695781808710, tolerance = 1e-10)
})

test_that("Gt is 1 before the first event and NA beyond the last observation", {
  d  <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  expect_equal(gt_(sv, min(d$time) / 2), 1, tolerance = 1e-12)
  # previously 0.2385 from a silent fallback; the censoring distribution is not
  # identified beyond the last observed time, so NA is correct (cf. pec::ipcw)
  expect_true(is.na(gt_(sv, max(d$time) * 2)))
})

test_that("Gt is pinned under heavily tied event times and stays monotone", {
  # 40 observations on 5 distinct times, alternating status. Under the previous
  # implementation this gave G(3) = 0.65625 > G(2.5) = 0.246094, impossible for
  # a survival function. It is now monotone.
  tt <- rep(c(1, 2, 3, 4, 5), each = 8)
  ss <- rep(c(1, 0), 20)
  sv <- survival::Surv(tt, ss)

  g <- c(gt_(sv, 1), gt_(sv, 2.5), gt_(sv, 3), gt_(sv, 4.5))
  expect_equal(g[1], 0.888888888889, tolerance = 1e-10)
  expect_equal(g[2], 0.761904761905, tolerance = 1e-10)
  expect_equal(g[3], 0.609523809524, tolerance = 1e-10)
  expect_equal(g[4], 0.406349206349, tolerance = 1e-10)
  expect_true(all(diff(g) <= 0))

  # the censoring distribution is exhausted at t = 5 on this fixture
  expect_error(gt_(sv, 5), "censoring distribution")
})

test_that("Gt is 1 throughout on the uncensored fixture", {
  # Previously a decaying sequence (0.7499, 0.6175, 0.3308) because
  # status == min(status) treated every observation as a censoring event.
  # With no censoring, G(t) is 1 by definition.
  du  <- fx_surv_uncensored()
  svu <- survival::Surv(du$time, du$status)
  qs  <- unname(stats::quantile(du$time, c(0.25, 0.5, 0.75)))

  for (q in qs) expect_equal(gt_(svu, q), 1, tolerance = 1e-12)
})

test_that("Gt rejects malformed input", {
  d  <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  expect_error(gt_(d$time, 1), "not of class Surv")
  expect_error(gt_(sv, -1), "must be positive")
  expect_error(gt_(sv, c(1, 2)), "single time point")
})
