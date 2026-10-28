# Competing-risks inputs. pa_predict_cif's return value is exactly the shape
# pa_summarize_cr requires (cif_pred, times, status), so the two are used as the
# designed pair rather than with synthetic CIF matrices.
cr_data <- function() {
  d <- fx_cr()
  d$time <- d$obs.times
  d$status <- d$obs.event
  d
}

cr_pred <- function() {
  dd <- cr_data()
  m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
                        data = dd, x = TRUE, y = TRUE)
  m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
                        data = dd, x = TRUE, y = TRUE)
  pa_predict_cif(model1 = m1, model2 = m2, newdata = dd,
                 covs = c("X1", "X2"), event.type = 1, tau = max(dd$time))
}

test_that("m_cif reduces a CIF column to a scalar risk", {
  cif <- seq(0, 0.6, length.out = 50)
  times <- seq(1, 50, length.out = 50)

  res <- PAmeasures:::m_cif(cif, time.cif = times, tau = 50)

  expect_length(res, 1L)
  expect_true(is.finite(res))
  expect_snapshot_value(snap_num(res), style = "serialize")
})

test_that("m_cif is monotone in the CIF magnitude", {
  times <- seq(1, 50, length.out = 50)

  low  <- PAmeasures:::m_cif(seq(0, 0.3, length.out = 50),
                             time.cif = times, tau = 50)
  high <- PAmeasures:::m_cif(seq(0, 0.9, length.out = 50),
                             time.cif = times, tau = 50)

  expect_false(isTRUE(all.equal(low, high, tolerance = 1e-6)))
  expect_snapshot_value(snap_num(c(low, high)), style = "serialize")
})

test_that("pa_survival_eval_cr returns a Metric/Value table", {
  p <- cr_pred()

  res <- pa_survival_eval_cr(
    pred_cif   = p$cif_pred[, -1],
    event_time = p$times,
    time.cif   = p$cif_pred[, 1],
    status     = p$status,
    event_type = 1
  )

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "Value"))
  # FINDING 14 fixed: Value is numeric here, matching the right-censored
  # pam.predicted_survial_eval. is.finite() now works without coercion.
  expect_type(res$Value, "double")
  expect_true(all(is.finite(res$Value)))
  # Keep historical metric snapshots unchanged; calibration has independent
  # reference and default-selection tests in test-calibration.R.
  res <- res[res$Metric != "ici", ]
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("both evaluation entry points return Value as double (FINDING 14 fixed)", {
  # The two entry points previously disagreed on the type of the Value column,
  # so is.finite() was silently FALSE for every competing-risks metric.
  p <- cr_pred()
  cr <- pa_survival_eval_cr(
    pred_cif = p$cif_pred[, -1], event_time = p$times,
    time.cif = p$cif_pred[, 1], status = p$status, event_type = 1
  )

  pred <- pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                               new_data = fx_surv(), tau = 10e10)
  rc <- pa_survival_eval(
    model = fx_cox(), event_time = pred$times,
    predicted_probability = pred$surv_prob, status = pred$status,
    covariates = fx_covs(), new_data = fx_surv(), tau = 10e10
  )

  expect_type(cr$Value, "double")
  expect_type(rc$Value, "double")
  expect_true(all(is.finite(cr$Value)))
  expect_true(all(is.finite(rc$Value)))
})

test_that("pred_cif is times-by-subjects with time.cif supplied separately", {
  # pa_summarize_cr splits pa_predict_cif's cif_pred as cif_pred[, 1] -> time.cif
  # and cif_pred[, -1] -> pred_cif, so pred_cif rows are times and columns are
  # subjects. Documented here because no roxygen states it.
  p <- cr_pred()
  pred_cif <- p$cif_pred[, -1]
  time_cif <- p$cif_pred[, 1]

  expect_identical(ncol(pred_cif), 200L)
  expect_identical(nrow(pred_cif), length(time_cif))
  expect_true(all(diff(time_cif) > 0))
})

test_that("competing-risks default metrics use C_index, not Harrell/Uno", {
  # default_metrics is c("pseudo_r2", "pseudo_r2_point", "c_index",
  #                      "brier_score", "td_auc")
  # This differs from the right-censored set and must survive standardization.
  p <- cr_pred()

  res <- pa_survival_eval_cr(
    pred_cif = p$cif_pred[, -1], event_time = p$times,
    time.cif = p$cif_pred[, 1], status = p$status, event_type = 1
  )

  expect_true("c_index" %in% res$Metric)
  expect_true("pseudo_r2" %in% res$Metric)
  expect_false("r_sh" %in% res$Metric)
  expect_false("r_e" %in% res$Metric)
  expect_false("harrell_c" %in% res$Metric)
})

test_that("pa_survival_eval_cr rejects a withdrawn metric name", {
  p <- cr_pred()

  expect_error(
    pa_survival_eval_cr(
      pred_cif = p$cif_pred[, -1], event_time = p$times,
      time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
      metrics = "r_sh"
    ),
    "withdrawn before the first CRAN release"
  )
})

test_that("pa_survival_eval_cr rejects an unknown metric name", {
  p <- cr_pred()

  expect_error(
    pa_survival_eval_cr(
      pred_cif = p$cif_pred[, -1], event_time = p$times,
      time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
      metrics = "not_a_metric"
    ),
    "Invalid metrics"
  )
})

test_that("pa_summarize_cr consumes pa_predict_cif output directly", {
  res <- pa_summarize_cr(list(csh = cr_pred()), event_type = 1)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "csh"))
  # default digits = 2
  expect_equal(res$csh, round(res$csh, 2), tolerance = 1e-12)
  res <- res[res$Metric != "ici", ]
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$csh), style = "serialize")
})

test_that("pa_summarize_cr names each model entry's required components", {
  # required <- c("cif_pred", "times", "status")
  p <- cr_pred()

  expect_error(
    pa_summarize_cr(list(bad = p[c("times", "status")]), event_type = 1),
    "must contain: cif_pred, times, status"
  )
  expect_error(pa_summarize_cr(list()), "must be a non-empty named list")
  expect_error(pa_summarize_cr("not a list"), "must be a non-empty named list")
})

test_that("pa_summarize_cr puts one column per model", {
  p <- cr_pred()

  res <- pa_summarize_cr(list(a = p, b = p), event_type = 1)

  expect_identical(names(res), c("Metric", "a", "b"))
  # identical inputs must give identical columns
  expect_equal(res$a, res$b, tolerance = 1e-12)
})

test_that("competing-risks status codes include a third level", {
  # Documents that CR functions accept {0,1,2}, unlike single-event functions.
  # Spec step 6 makes validation function-specific on this basis.
  d <- fx_cr()
  expect_setequal(sort(unique(d$obs.event)), c(0, 1, 2))
})
