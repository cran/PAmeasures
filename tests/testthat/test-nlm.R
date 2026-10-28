test_that("continuous-outcome measures agree with a hand-calculated example", {
  # Centered cross-product = 7; response SS = 10; prediction SS = 5.
  # Corrected SSE = 1/5 and original SSE = 2.
  out <- pa_nlm(c(1, 2, 4, 5), c(1, 2, 3, 4))
  expect_named(out, c("R.squared", "L.squared"))
  expect_type(out$R.squared, "double")
  expect_type(out$L.squared, "double")
  expect_equal(out$R.squared, 49/50, tolerance = 1e-14)
  expect_equal(out$L.squared, 1/10, tolerance = 1e-14)
})

test_that("supported pam.nlm results reproduce the archived implementation", {
  archived <- new.env(parent = asNamespace("stats"))
  sys.source(test_path("fixtures", "pam.nlm-0.1.0.R"), envir = archived)
  withr::local_seed(7201)
  y <- rnorm(80)
  predicted <- y^3 + rnorm(80, sd = .4)
  expect_identical(pam.nlm(y, predicted), archived$pam.nlm(y, predicted))
  expect_identical(pam.nlm(c(1, 2, 4, 5), 1:4),
                   list(R.squared = "0.9800", L.squared = "0.1000"))
  data("moore", package = "PAmeasures", envir = environment())
  fit <- stats::glm(log2(count) ~ time, data = moore, family = stats::gaussian())
  prediction <- 2^stats::predict(fit, type = "response")
  expect_identical(pam.nlm(moore$count, prediction),
                   archived$pam.nlm(moore$count, prediction))
  expect_warning(pam.nlm(y, predicted), NA)
})

test_that("the numeric interface retains precision and pairing", {
  withr::local_seed(7202)
  y <- rnorm(40)
  predicted <- sin(y) + rnorm(40, sd = .2)
  out <- pa_nlm(y, predicted)
  index <- sample.int(length(y))
  expect_equal(pa_nlm(y[index], predicted[index]), out, tolerance = 1e-13)
  expect_equal(pa_nlm(10*y + 3, 10*predicted + 3), out, tolerance = 1e-13)
  expect_gt(abs(out$R.squared - round(out$R.squared, 4)), 1e-8)
  expect_identical(pam.nlm(y, predicted), lapply(out, function(x)
    format(round(x, 4), nsmall = 4)))
})

test_that("constant predictions retain the intercept-only correction", {
  out <- pa_nlm(1:4, rep(4, 4))
  expect_identical(out$R.squared, 0)
  expect_equal(out$L.squared, 5/14, tolerance = 1e-14)
})

test_that("large finite integer inputs do not overflow during subtraction", {
  y <- c(-2147483647L, 2147483647L)
  expect_warning(out <- pa_nlm(y, rev(y)), NA)
  expect_equal(out$R.squared, 1)
  expect_equal(out$L.squared, 0, tolerance = 1e-14)
})

test_that("undefined measures warn without hiding the other defined measure", {
  expect_warning(constant <- pa_nlm(rep(3, 4), 1:4), "R.squared.*constant")
  expect_identical(constant$R.squared, NA_real_)
  expect_identical(constant$L.squared, 0)
  expect_warning(perfect <- pa_nlm(1:4, 1:4), "L.squared.*zero")
  expect_equal(perfect$R.squared, 1)
  expect_identical(perfect$L.squared, NA_real_)
  expect_warning(both <- pa_nlm(rep(3, 4), rep(3, 4)), "R.squared.*L.squared")
  expect_identical(both, list(R.squared = NA_real_, L.squared = NA_real_))
  expect_warning(legacy <- pam.nlm(1:4, 1:4), "L.squared")
  expect_identical(legacy$L.squared, "NA")
})

test_that("invalid inputs cannot be silently recycled or omitted", {
  for (bad in list("a", factor(1:4), list(1, 2), matrix(1:4, 2), 1+1i))
    expect_error(pa_nlm(bad, 1:4), "real numeric vectors")
  expect_error(pa_nlm(1:3, 1:4), "equal lengths")
  expect_error(pa_nlm(numeric(), numeric()), "at least two")
  expect_error(pa_nlm(1, 1), "at least two")
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(pa_nlm(c(1, bad), 1:2), "finite")
    expect_error(pa_nlm(1:2, c(1, bad)), "finite")
  }
  expect_error(pa_nlm(c(-1e200, 0, 1e200), rep(1, 3)), "overflow")
  expect_error(pa_nlm(c(1, 2, 3)*1e-200, rep(0, 3)), "underflow")
})
