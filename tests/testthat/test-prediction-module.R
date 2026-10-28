# Competing-risks inputs shared by the pa_predict_cif tests.
# pa_predict_cif reads newdata$time and newdata$status, so the fx_cr columns
# obs.times / obs.event must be aliased before it can be used.
pc_data <- function() {
  d <- fx_cr()
  d$time <- d$obs.times
  d$status <- d$obs.event
  d
}

pc_covs <- function() c("X1", "X2")

pc_cause_model <- function(cause) {
  dd <- pc_data()
  survival::coxph(
    survival::Surv(time, status == cause) ~ X1 + X2,
    data = dd, x = TRUE, y = TRUE
  )
}

test_that("pa_predict_coxph returns pred, times, status, surv_prob", {
  d <- fx_surv()

  res <- pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                              new_data = d, tau = 10e10)

  expect_type(res, "list")
  # component names the evaluation functions depend on; note "times", not "time"
  expect_true(all(c("pred", "times", "status", "surv_prob") %in% names(res)))
  expect_identical(length(res$times), 200L)
  expect_identical(length(res$status), 200L)
  expect_identical(length(res$pred), 200L)
  expect_true(all(res$status %in% c(0, 1)))
  expect_true(all(is.finite(res$times)))
  expect_identical(NROW(res$surv_prob), 200L)
  expect_true(all(res$surv_prob >= 0 & res$surv_prob <= 1, na.rm = TRUE))

  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
  expect_snapshot_value(mat_fingerprint(res$surv_prob), style = "serialize")
})

test_that("pa_predict_coxph requires time and status in new_data", {
  d <- fx_surv()

  expect_error(
    pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                         new_data = d[, c("x1", "x2")], tau = 10e10),
    "new_data require time and status columns"
  )
})

test_that("pa_predict_survreg returns the same component structure", {
  d <- fx_surv()

  res <- pa_predict_survreg(model = fx_survreg(), covs = fx_covs(),
                                new_data = d, tau = 10e10)

  expect_type(res, "list")
  expect_true(all(c("pred", "times", "status", "surv_prob") %in% names(res)))
  expect_identical(length(res$times), 200L)
  expect_identical(NROW(res$surv_prob), 200L)

  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
  expect_snapshot_value(mat_fingerprint(res$surv_prob), style = "serialize")
})

test_that("pa_predict_survreg validates model class and covariates", {
  d <- fx_surv()

  expect_error(
    pa_predict_survreg(model = fx_cox(), covs = fx_covs(), new_data = d),
    "model must be an object of class 'survreg'"
  )
  expect_error(
    pa_predict_survreg(model = fx_survreg(), covs = c("x1", "nope"),
                           new_data = d),
    "All covariates must be present in new_data"
  )
  expect_error(
    pa_predict_survreg(model = fx_survreg(), covs = fx_covs(),
                           new_data = d[, c("x1", "x2")]),
    "new_data require time and status columns"
  )
})

test_that("pa_predict_cif dispatches on a pair of cause-specific coxph fits", {
  dd <- pc_data()

  res <- pa_predict_cif(
    model1 = pc_cause_model(1), model2 = pc_cause_model(2),
    newdata = dd, covs = pc_covs(), event.type = 1, tau = max(dd$time)
  )

  expect_type(res, "list")
  expect_identical(sort(names(res)),
                   c("cif_pred", "linear.pred", "pred", "status", "times"))
  expect_identical(length(res$times), 200L)
  expect_identical(length(res$status), 200L)
  expect_identical(length(res$pred), 200L)
  expect_true(all(is.finite(res$pred)))

  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
  expect_snapshot_value(mat_fingerprint(res$cif_pred), style = "serialize")
})

test_that("cif_pred is times-by-subjects with the time grid in column 1", {
  # Documents the undocumented orientation that every downstream competing-risks
  # function assumes: m_pred <- cbind(time_grid, CIF), so rows are event times
  # and columns 2..n+1 are subjects.
  dd <- pc_data()

  res <- pa_predict_cif(
    model1 = pc_cause_model(1), model2 = pc_cause_model(2),
    newdata = dd, covs = pc_covs(), event.type = 1, tau = max(dd$time)
  )
  cif <- res$cif_pred

  expect_true(is.matrix(cif))
  # one column of times plus one column per subject
  expect_identical(ncol(cif), 201L)
  expect_gt(nrow(cif), 1)
  # column 1 is a strictly increasing time grid
  expect_true(all(diff(cif[, 1]) > 0))
  # the remaining columns are cumulative incidence: non-decreasing within [0, 1]
  expect_true(all(cif[, -1] >= 0 & cif[, -1] <= 1))
  expect_true(all(apply(cif[, -1], 2, function(col) all(diff(col) >= -1e-8))))
})

test_that("pa_predict_cif dispatches on a survreg model1", {
  # Both cause-specific models use the supported AFT prediction path.
  dd <- pc_data()
  m1 <- survival::survreg(
    survival::Surv(time, status == 1) ~ X1 + X2, data = dd, dist = "weibull",
    x = TRUE, y = TRUE
  )
  m2 <- survival::survreg(
    survival::Surv(time, status == 2) ~ X1 + X2, data = dd, dist = "weibull",
    x = TRUE, y = TRUE
  )

  res <- pa_predict_cif(
    model1 = m1, model2 = m2,
    newdata = dd, covs = pc_covs(), event.type = 1, tau = max(dd$time)
  )

  expect_type(res, "list")
  expect_true(all(c("times", "status", "cif_pred", "pred") %in% names(res)))
  expect_identical(length(res$pred), 200L)
  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
})

test_that("pa_predict_cif rejects an unrecognised model type", {
  dd <- pc_data()

  expect_error(
    pa_predict_cif(newdata = dd, covs = pc_covs(), event.type = 1),
    "Unknown model type"
  )
  expect_error(
    pa_predict_cif(model1 = "not a model", newdata = dd, covs = pc_covs()),
    "Unknown model type"
  )
})
