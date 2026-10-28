# r_e and r_sh were withdrawn before the first CRAN release. Every spelling
# that used to resolve to them must now raise an error naming the withdrawal,
# rather than a generic "invalid metric" error, so a user with an existing
# script learns why it broke.
#
# The check lives in pa_normalize_metrics(), which pa_survival_eval(),
# pa_fit_and_eval(), pa_survival_eval_cr(), pa_evaluate_two_phase() and
# pa_sample_design() all call, so one implementation covers every public
# entry point.

tombstone_re <- "withdrawn before the first CRAN release"

test_that("every removed spelling raises the tombstone error", {
  for (nm in c("r_e", "r_sh", "R_E", "R_sh", "R_sph", "R_SH", "r e", "R_sPh")) {
    expect_error(PAmeasures:::pa_normalize_metrics(nm), tombstone_re,
                 info = paste0("spelling: ", nm))
  }
})

test_that("the tombstone fires even when mixed with valid metrics", {
  expect_error(
    PAmeasures:::pa_normalize_metrics(c("harrell_c", "r_e", "brier_score")),
    tombstone_re
  )
})

test_that("surviving metric names still normalise unchanged", {
  expect_identical(PAmeasures:::pa_normalize_metrics("harrell_c"), "harrell_c")
  expect_identical(
    suppressWarnings(PAmeasures:::pa_normalize_metrics("Harrells_C")),
    "harrell_c"
  )
  expect_null(PAmeasures:::pa_normalize_metrics(NULL))
})

test_that("an unrelated unknown name is NOT claimed by the tombstone", {
  # falls through unchanged, for the caller's own validation to report
  expect_identical(
    PAmeasures:::pa_normalize_metrics("not_a_metric", warn = FALSE),
    "not_a_metric"
  )
})

test_that("pa_metric_names() no longer offers the removed metrics", {
  nms <- pa_metric_names()
  expect_false("r_e" %in% nms)
  expect_false("r_sh" %in% nms)
  expect_length(nms, 14L)
  expect_identical(nms, sort(c(
    "brier_score", "c_index", "harrell_c", "l2_point", "l_square",
    "pseudo_r2", "pseudo_r2_point", "r2_point", "r_square", "td_auc", "uno_c",
    "ici", "e50", "e90"
  )))
})

# --- the tombstone must reach every public entry point -----------------------
# pa_survival_eval_cr() used to carry a duplicate validation block that ran
# BEFORE pa_normalize_metrics(). It rejected every legacy spelling on this path
# ("Harrells_C", "C_index", "Brier Score") while they worked everywhere else,
# and it intercepted r_e / r_sh before the withdrawal message could be raised.
# The duplicate was removed; these tests pin both halves of that fix.

cr_tomb_pred <- function() {
  dd <- fx_cr()
  dd$time <- dd$obs.times
  dd$status <- dd$obs.event
  m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
                        data = dd, x = TRUE, y = TRUE)
  m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
                        data = dd, x = TRUE, y = TRUE)
  pa_predict_cif(model1 = m1, model2 = m2, newdata = dd,
                 covs = c("X1", "X2"), event.type = 1, tau = max(dd$time))
}

cr_tomb_eval <- function(metrics) {
  p <- cr_tomb_pred()
  pa_survival_eval_cr(
    pred_cif = p$cif_pred[, -1], event_time = p$times,
    time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
    metrics = metrics
  )
}

test_that("pa_survival_eval_cr raises the tombstone, not a generic error", {
  expect_error(cr_tomb_eval("r_sh"), tombstone_re)
  expect_error(cr_tomb_eval("R_E"), tombstone_re)
})

test_that("pa_survival_eval_cr accepts legacy spellings of surviving metrics", {
  res <- suppressWarnings(cr_tomb_eval("C_index"))
  expect_identical(res$Metric, "c_index")

  res2 <- suppressWarnings(cr_tomb_eval("Brier Score"))
  expect_identical(res2$Metric, "brier_score")
})

test_that("pa_survival_eval_cr still rejects a genuinely unknown metric", {
  expect_error(cr_tomb_eval("not_a_metric"), "Invalid metrics")
})

# --- the evaluators must no longer produce the withdrawn metrics -------------

tomb_surv_eval <- function(metrics = NULL) {
  pred <- pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                           new_data = fx_surv(), tau = 10e10)
  pa_survival_eval(
    model = fx_cox(), event_time = pred$times,
    predicted_probability = pred$surv_prob, status = pred$status,
    covariates = fx_covs(), new_data = fx_surv(), tau = 10e10,
    metrics = metrics
  )
}

test_that("pa_survival_eval no longer emits the removed metrics", {
  res <- tomb_surv_eval(NULL)

  expect_false("r_e" %in% res$Metric)
  expect_false("r_sh" %in% res$Metric)
  expect_true(all(c("pseudo_r2", "harrell_c", "uno_c",
                    "brier_score", "td_auc") %in% res$Metric))
})

test_that("pa_survival_eval with metrics = 'all' omits the removed metrics", {
  res <- tomb_surv_eval("all")
  expect_false(any(c("r_e", "r_sh") %in% res$Metric))
})

test_that("pa_fit_and_eval no longer produces the removed metric columns", {
  res <- pa_fit_and_eval(train_data = fx_surv(), covariates = fx_covs(),
                         models = "coxph", metrics = "all")

  expect_false("r_e" %in% names(res))
  expect_false("r_sh" %in% names(res))
  expect_true(all(c("pseudo_r2", "r_square", "l_square",
                    "brier_score") %in% names(res)))
})

test_that("pa_fit_and_eval rejects the removed metrics by name", {
  expect_error(
    pa_fit_and_eval(train_data = fx_surv(), covariates = fx_covs(),
                    models = "coxph", metrics = "r_sh"),
    tombstone_re
  )
})

test_that("the borrowed implementations are gone from the namespace", {
  ns <- asNamespace("PAmeasures")
  for (nm in c("pam.rsph", "pam.rsph.coxph", "pam.rsph.aareg",
               "pam.rsph.survreg", "print.rsph", "summary.rsph",
               "pam.schemper", "my.survfit", "pam.re")) {
    expect_false(exists(nm, envir = ns, inherits = FALSE),
                 info = paste0("still present: ", nm))
  }
})
