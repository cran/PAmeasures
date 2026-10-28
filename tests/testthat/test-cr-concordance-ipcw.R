# Independent pair enumeration using hand-computed censoring probabilities.
cr_concordance_reference <- function(time, status, prediction, tau, left, right) {
  numerator <- denominator <- 0
  for (i in which(status == 1 & time <= tau)) {
    for (j in seq_along(time)) {
      probability <- if (time[i] < time[j]) left[i] * right[i] else
        if (status[j] == 2) left[i] * left[j] else NA_real_
      if (is.na(probability)) next
      denominator <- denominator + 1 / probability
      numerator <- numerator + (prediction[i] < prediction[j]) / probability
    }
  }
  numerator / denominator
}

test_that("CR concordance uses reciprocal products for both types of pair", {
  time <- 1:8
  status <- c(1, 0, 1, 2, 0, 1, 2, 1)
  prediction <- c(2, 4, 1, 8, 3, 5, 7, 6)
  # Censoring at 2 and 5 changes G by factors 6/7 and 3/4.
  left <- c(1, 1, rep(6/7, 3), rep(9/14, 3))
  right <- c(1, rep(6/7, 3), rep(9/14, 4))
  expected <- cr_concordance_reference(time, status, prediction, 7, left, right)
  actual <- PAmeasures:::C_cr(time, status, prediction, tau = 7)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(actual, 0.951123807367981, tolerance = 1e-7)
  # The old product weighting returned 0.9170195 for these same pairs.
  expect_gt(abs(actual - 0.9170195), .03)
})

test_that("CR concordance includes the horizon and ignores unused zero G", {
  time <- 1:8
  status <- c(1, 0, 1, 2, 0, 1, 2, 0)
  prediction <- c(2, 4, 1, 8, 3, 5, 7, 6)
  left <- c(1, 1, rep(6/7, 3), rep(9/14, 3))
  right <- c(1, rep(6/7, 3), rep(9/14, 3), 0)
  for (tau in c(3, 6, 7)) {
    expected <- cr_concordance_reference(time, status, prediction, tau, left, right)
    actual <- PAmeasures:::C_cr(time, status, prediction, tau = tau)
    expect_true(is.finite(actual))
    expect_equal(actual, expected, tolerance = 1e-12)
  }
})

test_that("uncensored CR concordance reduces to unweighted comparable pairs", {
  time <- c(1, 2, 2, 3, 4, 5)
  status <- c(1, 1, 2, 1, 2, 1)
  prediction <- c(2, 1, 4, 3, 6, 5)
  expected <- cr_concordance_reference(time, status, prediction, 3,
                                      rep(1, 6), rep(1, 6))
  expect_equal(PAmeasures:::C_cr(time, status, prediction, tau = 3), expected)
  order <- c(6, 3, 1, 5, 2, 4)
  expect_equal(PAmeasures:::C_cr(time[order], status[order], prediction[order], tau = 3),
               expected)
})

test_that("public CR evaluation and summary use the corrected concordance", {
  time <- 1:8
  status <- c(1, 0, 1, 2, 0, 1, 2, 0)
  grid <- c(1, 3, 6)
  cif <- outer(c(.1, .3, .5), c(.2, .4, .1, .8, .3, .5, .7, .6))
  prediction <- apply(cif, 2, PAmeasures:::m_cif, time.cif = grid, tau = 6)
  left <- c(1, 1, rep(6/7, 3), rep(9/14, 3))
  right <- c(1, rep(6/7, 3), rep(9/14, 3), 0)
  expected <- round(cr_concordance_reference(time, status, prediction, 6,
                                             left, right), 4)
  actual <- pa_survival_eval_cr(cif, time, grid, status,
    metrics = "c_index", tau = 6, t_star = 6)
  expect_equal(actual$Value, expected)
  model <- list(cif_pred = cbind(grid, cif), times = time, status = status)
  summary <- pa_summarize_cr(list(model = model), metrics = "c_index",
    tau = 6, t_star = 6, digits = 4)
  expect_equal(summary$model, expected)
})
