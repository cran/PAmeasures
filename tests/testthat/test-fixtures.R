test_that("fx_surv is deterministic and correctly shaped", {
  d1 <- fx_surv()
  d2 <- fx_surv()

  expect_s3_class(d1, "data.frame")
  expect_identical(names(d1), c("time", "status", "x1", "x2"))
  expect_identical(nrow(d1), 200L)
  expect_identical(d1, d2)
  expect_true(all(d1$status %in% c(0, 1)))
  expect_true(all(d1$time > 0))
  expect_gt(sum(d1$status == 0), 0)
})

test_that("fx_surv_uncensored has no censoring", {
  d <- fx_surv_uncensored()

  expect_identical(names(d), c("time", "status", "x1", "x2"))
  expect_identical(nrow(d), 200L)
  expect_true(all(d$status == 1))
})

test_that("fx_cr produces competing-risk codes", {
  d <- fx_cr()

  expect_s3_class(d, "data.frame")
  expect_identical(nrow(d), 200L)
  expect_true(all(c("obs.times", "obs.event") %in% names(d)))
  expect_true(all(d$obs.event %in% c(0, 1, 2)))
  expect_gt(sum(d$obs.event == 2), 0)
  expect_true(all(d$obs.times > 0))
})

test_that("fixture models fit and are reproducible", {
  m <- fx_cox()
  expect_s3_class(m, "coxph")
  expect_identical(names(coef(m)), c("x1", "x2"))
  expect_equal(coef(m), coef(fx_cox()), tolerance = 1e-12)
  expect_false(is.null(m$x))
  expect_false(is.null(m$y))

  s <- fx_survreg()
  expect_s3_class(s, "survreg")
  expect_identical(names(coef(s)), c("(Intercept)", "x1", "x2"))
})

test_that("fixture numbers are snapshot-stable", {
  d <- fx_surv()

  expect_snapshot_value(snap_num(head(d$time, 10)), style = "serialize")
  expect_snapshot_value(sum(d$status), style = "serialize")
  expect_snapshot_value(snap_num(unname(coef(fx_cox()))), style = "serialize")
  expect_snapshot_value(snap_num(unname(coef(fx_survreg()))), style = "serialize")
  expect_snapshot_value(table(fx_cr()$obs.event), style = "serialize")
})
