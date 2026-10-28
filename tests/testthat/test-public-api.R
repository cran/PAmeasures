test_that("the supported API uses pa_ with pam.nlm as its sole legacy exception", {
  exports <- getNamespaceExports("PAmeasures")
  expected <- c(
    "pa_nlm", "pa_survival_eval", "pa_survival_eval_cr",
    "pa_evaluate_two_phase", "pa_fit_and_eval", "pa_metric_names",
    "pa_predict_coxph", "pa_predict_survreg", "pa_predict_cif",
    "pa_summarize", "pa_summarize_cr", "pa_sample_design",
    "pa_case_cohort_weights", "pa_nested_case_control_weights",
    "pa_plot_pred", "pa_plot_summary", "pa_sim_cox_weibull",
    "pa_simulate_fine_gray", "pa_calibration", "pa_calibration_cr",
    "pam.nlm"
  )
  expect_setequal(exports, expected)
  expect_identical(exports[!startsWith(exports, "pa_")], "pam.nlm")
})

test_that("retired TimeMetric script names are absent from the namespace", {
  retired <- c(
    "pam.predicted_survial_eval", "pam.predicted_survial_eval_cr",
    "pam.survival_eval", "pam.coxph_restricted", "pam.surverg_restricted",
    "pam.predict_cr", "pam.summary", "pam.summary_cr", "pam.sample_design",
    "cc_weights", "ncc_weights", "plot_pred", "summary_pred_plot",
    "sim_cox_weibull_censored", "simulateTwoCauseFineGrayModel"
  )
  ns <- asNamespace("PAmeasures")
  expect_false(any(startsWith(getNamespaceExports("PAmeasures"), "tm_")))
  expect_false(any(vapply(retired, exists, logical(1), envir = ns,
                          inherits = FALSE)))
  expect_true(is.function(getExportedValue("PAmeasures", "pam.nlm")))
})
