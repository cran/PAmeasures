test_that("Gt returns a censoring survival probability in [0, 1]", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)

  res <- PAmeasures:::Gt(sv, stats::median(d$time))

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("Gt is non-increasing in the timepoint", {
  d <- fx_surv()
  sv <- survival::Surv(d$time, d$status)
  qs <- stats::quantile(d$time, c(0.2, 0.4, 0.6, 0.8))

  vals <- vapply(qs, function(tp) PAmeasures:::Gt(sv, tp), numeric(1))

  expect_true(all(diff(vals) <= 1e-8))
  expect_snapshot_value(snap_num(vals), style = "serialize")
})

test_that("Gt rejects non-Surv input", {
  expect_error(PAmeasures:::Gt(1:10, 5), "not of class Surv")
})

test_that("pam.Brier on a coxph fit returns a finite scalar in [0, 1]", {
  # This is the chain that reaches pec::predictSurvProb unqualified.
  d <- fx_surv()
  fit <- fx_cox()

  res <- PAmeasures:::pam.Brier(fit, d, stats::median(d$time))

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_gte(res, 0)
  expect_lte(res, 1)
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("pam.Brier default t_star is the median of TRAINING event times", {
  # source R/pam.Brier.R:59-63 --
  #   distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))
  #   if (t_star0 <= 0) t_star0 <- median(distime)
  # This is the median of the model's own event times, NOT median(d$time).
  d <- fx_surv()
  fit <- fx_cox()

  distime <- sort(unique(as.vector(fit$y[fit$y[, 2] == 1])))
  expected_t <- stats::median(distime)

  expect_equal(
    PAmeasures:::pam.Brier(fit, d),
    PAmeasures:::pam.Brier(fit, d, expected_t),
    tolerance = 1e-6
  )

  # and confirm this genuinely differs from the naive test-data median,
  # so the distinction is protected against a future "simplification"
  expect_false(isTRUE(all.equal(expected_t, stats::median(d$time),
                                tolerance = 1e-6)))
})

test_that("pec::predictSurvProb drives pam.Brier and is reachable", {
  skip_if_not_installed("pec")
  d <- fx_surv()
  fit <- fx_cox()
  t_star <- stats::median(d$time)

  probs <- pec::predictSurvProb(fit, d, t_star)

  expect_identical(NROW(probs), 200L)
  expect_true(all(probs >= 0 & probs <= 1))
  # fingerprint, not the full 200-element matrix
  expect_snapshot_value(mat_fingerprint(probs), style = "serialize")
})
