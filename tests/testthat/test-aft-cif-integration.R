# Independent distribution identities catch the lost probability mass at zero.
aft_cif_fixture <- function() {
  set.seed(7183)
  data.frame(time = exp(rnorm(80, 7, 1)),
             status = rep(c(0L, 1L, 2L, 1L), 20), x = rnorm(80))
}

test_that("Weibull CIFs match the closed form, including scale five", {
  d <- aft_cif_fixture()
  nd <- data.frame(time = c(0, 0.01, 1, 100, 4380), status = 0L)
  for (sc in c(0.5, 1, 5)) {
    a <- survival::survreg(survival::Surv(time, status == 1) ~ 1,
      data = d, dist = "weibull", scale = sc)
    b <- survival::survreg(survival::Surv(time, status == 2) ~ 1,
      data = d, dist = "weibull", scale = sc)
    p <- pa_predict_cif(model1 = a, model2 = b, newdata = nd,
                        covs = character(), tau = 4380)
    rate1 <- exp(-unname(stats::coef(a)) / sc)
    rate2 <- exp(-unname(stats::coef(b)) / sc)
    expected <- rate1 / (rate1 + rate2) *
      (-expm1(-(rate1 + rate2) * p$cif_pred[, 1]^(1 / sc)))
    expect_equal(unname(p$cif_pred[, -1]),
      matrix(expected, nrow = length(expected), ncol = nrow(nd)), tolerance = 1e-8)
    expect_true(all(p$cif_pred[, -1] >= 0 & p$cif_pred[, -1] <= 1))
  }
})

test_that("AFT CIF values do not depend on intervening evaluation times", {
  d <- aft_cif_fixture()
  a <- survival::survreg(survival::Surv(time, status == 1) ~ x,
    data = d, dist = "weibull", scale = 5)
  b <- survival::survreg(survival::Surv(time, status == 2) ~ x,
    data = d, dist = "weibull", scale = 0.7)
  coarse <- data.frame(time = c(100, 4380), status = 0L, x = c(-0.5, 0.5))
  fine <- rbind(coarse, data.frame(time = seq(1, 4300, length.out = 50), status = 0L, x = 0))
  p <- pa_predict_cif(model1 = a, model2 = b, newdata = coarse, covs = "x")
  q <- pa_predict_cif(model1 = a, model2 = b, newdata = fine, covs = "x")
  expect_equal(p$cif_pred[, 2:3], q$cif_pred[match(p$cif_pred[, 1], q$cif_pred[, 1]), 2:3],
               tolerance = 1e-8)
})

test_that("lognormal and loglogistic CIFs obey the equal-hazard identity", {
  d <- aft_cif_fixture()
  nd <- data.frame(time = c(0, 1, 100, 4380), status = 0L, x = c(-1, 0, 1, 2))
  for (dist in c("lognormal", "loglogistic")) {
    a <- survival::survreg(survival::Surv(time, status == 1) ~ x,
      data = d, dist = dist)
    p <- pa_predict_cif(model1 = a, model2 = a, newdata = nd, covs = "x")
    lp <- stats::predict(a, newdata = nd, type = "lp")
    z <- outer(log(p$cif_pred[, 1]), lp, "-") / a$scale
    s <- if (dist == "lognormal") stats::pnorm(z, lower.tail = FALSE) else
      stats::plogis(z, lower.tail = FALSE)
    expect_equal(unname(p$cif_pred[, -1]), unname(0.5 * (1 - s^2)), tolerance = 1e-8)
  }
})

test_that("unequal Weibull hazards agree with integration on the time scale", {
  d <- aft_cif_fixture()
  a <- survival::survreg(survival::Surv(time, status == 1) ~ x,
    data = d, dist = "weibull", scale = 0.5)
  b <- survival::survreg(survival::Surv(time, status == 2) ~ x,
    data = d, dist = "weibull", scale = 0.8)
  nd <- data.frame(time = c(100, 4380), status = 0L, x = c(-1, 1))
  p <- pa_predict_cif(model1 = a, model2 = b, newdata = nd, covs = "x")
  sc1 <- exp(stats::predict(a, newdata = nd, type = "lp"))
  sc2 <- exp(stats::predict(b, newdata = nd, type = "lp"))
  expected <- sapply(seq_len(nrow(nd)), function(i) vapply(p$cif_pred[, 1], function(t) {
    stats::integrate(function(u) stats::dweibull(u, 1/a$scale, sc1[i]) *
      stats::pweibull(u, 1/b$scale, sc2[i], lower.tail = FALSE), 0, t,
      rel.tol = 1e-10)$value
  }, numeric(1)))
  expect_equal(unname(p$cif_pred[, -1]), unname(expected), tolerance = 1e-8)
})
