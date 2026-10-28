test_that("package exports the pa_ API", {
  exports <- getNamespaceExports("PAmeasures")

  expect_true(all(c(
    "pa_survival_eval", "pa_survival_eval_cr", "pa_evaluate_two_phase",
    "pa_fit_and_eval", "pa_predict_coxph", "pa_predict_survreg",
    "pa_predict_cif", "pa_summarize", "pa_summarize_cr", "pa_sample_design",
    "pa_case_cohort_weights", "pa_nested_case_control_weights",
    "pa_plot_pred", "pa_plot_summary", "pa_sim_cox_weibull",
    "pa_simulate_fine_gray", "pa_calibration", "pa_calibration_cr"
  ) %in% exports))
  # Every primary function uses pa_; pam.nlm is the retained compatibility entry.
  new_api <- grep("^pa_", exports, value = TRUE)
  expect_length(new_api, 20L)   # nineteen survival interfaces plus pa_nlm
  expect_true("pa_metric_names" %in% exports)
})

test_that("the deleted Cluster C duplicates are gone", {
  # pam.Brier_metric, pam.rsh_metric and pam.rsph_metric were removed after the
  # equivalence gate and the R_E audit. See docs/superpowers/ for the evidence.
  ns <- asNamespace("PAmeasures")

  expect_false(exists("pam.Brier_metric", envir = ns, inherits = FALSE))
  expect_false(exists("pam.rsh_metric", envir = ns, inherits = FALSE))
  expect_false(exists("pam.rsph_metric", envir = ns, inherits = FALSE))
})

test_that("functions this suite reaches with ::: are genuinely internal", {
  exports <- getNamespaceExports("PAmeasures")

  expect_false(any(c(
    "Gt", "pam.Brier", "m_cif"
  ) %in% exports))
  # pa_evaluate_two_phase was promoted from internal to public, so the
  # case-cohort and NCC functionality the paper advertises is now reachable
  expect_true("pa_evaluate_two_phase" %in% exports)
})
