# Finding 42: the evaluation-time index must be chosen deterministically.
#
# `which.min(abs(event_time - t_star))` broke when t_star was the median of an
# even-length vector: the two middle order statistics are mathematically
# equidistant from it, so a one-ulp rounding difference decided the winner. That
# rounding differs between 80-bit-long-double platforms (x86_64 Linux, Windows)
# and plain-double ones (arm64 macOS), so the same data yielded different
# point-in-time metrics on different machines.
#
# The rule is now: nearest observed time, ties broken toward the EARLIER event
# time -- which also makes the choice independent of input row order.

nti <- function(...) PAmeasures:::pa_nearest_time_index(...)

test_that("an exact midpoint between two times resolves to the earlier one", {
  et <- c(1, 2, 3, 4)
  # exactly halfway between 2 and 3: mathematically tied
  expect_identical(nti(et, 2.5), 2L)
  expect_identical(et[nti(et, 2.5)], 2)
})

test_that("a near midpoint resolves to the genuinely nearer time on each side", {
  et <- c(1, 2, 3, 4)
  expect_identical(et[nti(et, 2.5 - 1e-6)], 2)  # just below -> earlier
  expect_identical(et[nti(et, 2.5 + 1e-6)], 3)  # just above -> later
})

test_that("the tie rule is independent of input row order", {
  set.seed(4242)
  et <- c(1, 2, 3, 4)
  for (i in 1:20) {
    p <- sample(length(et))
    # the chosen TIME must be identical regardless of row order
    expect_identical(et[p][nti(et[p], 2.5)], 2,
                     info = paste("permutation", i))
  }
})

test_that("duplicated times do not make the choice order-dependent", {
  et <- c(2, 2, 3, 3)
  expect_identical(et[nti(et, 2.5)], 2)
  set.seed(7)
  for (i in 1:10) {
    p <- sample(length(et))
    expect_identical(et[p][nti(et[p], 2.5)], 2)
  }
})

test_that("an exactly observed t_star selects that time", {
  et <- c(1, 2, 3, 4)
  for (t in et) expect_identical(et[nti(et, t)], t)
})

# --- end-to-end, through the public evaluator -------------------------------

tie_data <- function() {
  # n = 200 is even, so the default t_star = quantile(time, 0.5) falls exactly
  # between the 100th and 101st order statistics: the tie this rule fixes.
  d <- fx_surv()
  stopifnot(nrow(d) %% 2 == 0)
  d
}

eval_with <- function(d, t_star = NULL) {
  fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         data = d, x = TRUE, y = TRUE)
  pred <- pa_predict_coxph(model = fit, covs = c("x1", "x2"),
                           new_data = d, tau = 10e10)
  suppressMessages(pa_survival_eval(
    model = fit, event_time = pred$times,
    predicted_probability = pred$surv_prob, status = pred$status,
    covariates = c("x1", "x2"), new_data = d, tau = 10e10,
    t_star = t_star, metrics = "all"))
}

test_that("metrics are invariant to input row order at the default t_star", {
  d <- tie_data()
  set.seed(99)
  a <- eval_with(d)
  b <- eval_with(d[sample(nrow(d)), ])

  expect_identical(a$Metric, b$Metric)
  expect_equal(as.numeric(a$Value), as.numeric(b$Value), tolerance = 1e-10)
})

test_that("metrics are invariant to input row order at an explicit t_star", {
  d <- tie_data()
  ts <- stats::quantile(d$time, 0.4)
  set.seed(101)
  a <- eval_with(d, t_star = ts)
  b <- eval_with(d[sample(nrow(d)), ], t_star = ts)

  expect_equal(as.numeric(a$Value), as.numeric(b$Value), tolerance = 1e-10)
})

test_that("the default t_star really does sit exactly between two times", {
  # guards the regression: if this stops being a tie, the test above stops
  # exercising the defect and should be rebuilt on data that is still tied.
  d <- tie_data()
  s <- sort(d$time)
  n <- length(s)
  ts <- stats::quantile(d$time, 0.5)
  lo <- s[n / 2]
  hi <- s[n / 2 + 1]
  expect_false(ts %in% d$time)
  # the two middle times are equidistant to within a few ulp
  expect_lt(abs(abs(lo - ts) - abs(hi - ts)), 1e-12)
  # and the rule picks the earlier of them
  expect_identical(PAmeasures:::pa_nearest_time_index(sort(d$time), ts),
                   as.integer(n / 2))
})
