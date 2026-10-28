# Shared prediction input for the evaluation tests.
ev_pred <- function() {
  pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                       new_data = fx_surv(), tau = 10e10)
}

ev_result <- function() {
  pred <- ev_pred()
  pa_survival_eval(
    model = fx_cox(),
    event_time = pred$times,
    predicted_probability = pred$surv_prob,
    status = pred$status,
    covariates = fx_covs(),
    new_data = fx_surv(),
    tau = 10e10
  )
}

test_that("pa_survival_eval returns a Metric/Value data frame", {
  res <- ev_result()

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "Value"))
  expect_type(res$Value, "double")
  # Preserve the pre-calibration snapshots; the new default ICI and its
  # independent numerical references are covered in test-calibration.R.
  res <- res[res$Metric != "ici", ]
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$Value), style = "serialize")
})

test_that("the default metric set appears in the Metric column", {
  # Metric names live in res$Metric; names(res) is c("Metric", "Value").
  # r_sh and r_e were withdrawn before the first CRAN release; the assertions
  # that they are absent and rejected live in test-metric-tombstone.R.
  res <- ev_result()

  expect_true("brier_score" %in% res$Metric)
  expect_true("pseudo_r2" %in% res$Metric)
  expect_true("td_auc" %in% res$Metric)
  expect_true("harrell_c" %in% res$Metric)
  expect_true("uno_c" %in% res$Metric)
  # the legacy spellings are no longer emitted
  expect_false("Harrells_C" %in% res$Metric)
  expect_false("Harrell\u2019s C" %in% res$Metric)
  expect_false("Pesudo_R" %in% res$Metric)
})

test_that("pa_survival_eval rejects an unknown metric name", {
  pred <- ev_pred()

  expect_error(
    pa_survival_eval(
      model = fx_cox(),
      event_time = pred$times,
      predicted_probability = pred$surv_prob,
      status = pred$status,
      covariates = fx_covs(),
      new_data = fx_surv(),
      tau = 10e10,
      metrics = "not_a_real_metric"
    ),
    "Invalid metrics"
  )
})

test_that("concordancefit is imported (FINDING 11 fixed)", {
  # R/pam.predicted_survial_eval.R:162,169 call concordancefit() unqualified.
  # It is now imported explicitly, so the call resolves from the namespace
  # rather than depending on survival happening to be attached.
  imported <- unlist(getNamespaceImports("PAmeasures"), use.names = FALSE)

  expect_true("concordancefit" %in% imported)
})

test_that("a clean session without survival attached now succeeds (FINDING 11)", {
  # The definitive check: a fresh R subprocess that loads PAmeasures and nothing
  # else -- the state a first-time user is in. This previously died with
  # "could not find function concordancefit"; it must now produce metrics.
  skip_on_cran()
  skip_if_covr()
  pkg_root <- skip_without_source_tree()

  script <- sprintf('
    suppressWarnings(pkgload::load_all(%s, quiet = TRUE, attach_testthat = FALSE))
    stopifnot(!"package:survival" %%in%% search())
    d <- pa_sim_cox_weibull(n = 50, pi_c = 0.3, v = 2,
                                  beta = c(0.5, -0.5), seed = 1001)
    d <- d[, c("time", "status", "x1", "x2")]
    m <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         data = d, x = TRUE, y = TRUE)
    p <- pa_predict_coxph(model = m, covs = c("x1", "x2"),
                              new_data = d, tau = 10e10)
    r <- try(pa_survival_eval(
      model = m, event_time = p$times, predicted_probability = p$surv_prob,
      status = p$status, covariates = c("x1", "x2"), new_data = d, tau = 10e10
    ), silent = TRUE)
    cat(if (inherits(r, "try-error")) as.character(r) else "NO ERROR")
  ', shQuote(pkg_root))

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", "-e", shQuote(script)),
    stdout = TRUE, stderr = TRUE
  ))
  out <- paste(out, collapse = "\n")

  expect_no_match(out, "could not find function")
  expect_no_match(out, "Error")
  expect_match(out, "NO ERROR")   # sentinel printed only on success
})

test_that("pa_summarize pivots one model into a Metric column table", {
  res <- pa_summarize(list(value = ev_pred()), tau = 10e10)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "value"))
  res <- res[res$Metric != "ici", ]
  expect_snapshot_value(res$Metric, style = "serialize")
  expect_snapshot_value(snap_num(res$value), style = "serialize")
})

test_that("pa_summarize puts one column per model and rounds to digits", {
  d <- fx_surv()
  p_cox <- pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                                new_data = d, tau = 10e10)
  p_reg <- pa_predict_survreg(model = fx_survreg(), covs = fx_covs(),
                                  new_data = d, tau = 10e10)

  res <- pa_summarize(list(cox = p_cox, weibull = p_reg), tau = 10e10)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "cox", "weibull"))
  # default digits = 2
  expect_equal(res$cox, round(res$cox, 2), tolerance = 1e-12)
  res <- res[res$Metric != "ici", ]
  expect_snapshot_value(snap_num(res$cox), style = "serialize")
  expect_snapshot_value(snap_num(res$weibull), style = "serialize")
})

test_that("pa_summarize rejects a non-list or empty models argument", {
  expect_error(pa_summarize(list()), "must be a non-empty named list")
  expect_error(pa_summarize("not a list"), "must be a non-empty named list")
})

test_that("pa_fit_and_eval fits and evaluates from raw data (FINDING 12 fixed)", {
  # R/pam.survial_eval.R:106,112 previously passed covariates=/newdata= to
  # functions taking covs=/new_data=, so every call failed. They now also pass
  # predict = FALSE, which is the branch returning R.squared / L.squared that
  # the caller indexes as r_l_list[1] and [2].
  d <- fx_surv()

  res <- pa_fit_and_eval(train_data = d, covariates = fx_covs(),
                           models = "coxph", metrics = "all")

  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 1L)
  expect_true("Model" %in% names(res))
  expect_identical(res$Model, "coxph")
  # wide layout: one column per metric, unlike the long Metric/Value frame
  # returned by pa_survival_eval
  expect_true(all(c("pseudo_r2", "r_square", "l_square",
                    "brier_score") %in% names(res)))
  expect_snapshot_value(sort(setdiff(names(res), c("ici", "e50", "e90"))),
                          style = "serialize")
  expect_snapshot_value(snap_num(res[["r_square"]]), style = "serialize")
})

test_that("pa_fit_and_eval honours an explicit model and metric subset", {
  d <- fx_surv()

  res <- pa_fit_and_eval(
    train_data = d, covariates = fx_covs(),
    models  = c("weibull", "lognormal"),
    metrics = c("r_square", "l_square", "brier_score")
  )

  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 2L)
  expect_setequal(res$Model, c("weibull", "lognormal"))
  expect_false("harrell_c" %in% names(res))
})

test_that("pa_fit_and_eval requires train_data and covariates", {
  expect_error(
    pa_fit_and_eval(),
    "Please provide 'train_data', 'time_var', 'status_var', and 'covariates'"
  )
})
