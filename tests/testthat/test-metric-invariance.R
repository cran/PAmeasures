# Numerical invariance gate for the r_sh / r_e removal.
# Audit amendment 2026-10-11: four two_phase_ncc expectations and three
# right-censored point-metric expectations were updated after correcting
# m-controls probabilities and the event/survival probability mismatch. See
# test-api-selection.R for independent reference checks and
# docs/validation/results/metric-baseline-before.csv for the original values.
#
# fixtures/metric-baseline.csv holds committed literal values generated from
# the package BEFORE either metric was removed, at commit 675c99f, using the
# deterministic fixtures in helper-simdata.R. Every metric that survives the
# removal must keep its exact value on every public evaluator path.
#
# Nine scenarios cover all four public entry points that resolve metric names
# through pa_normalize_metrics(), because the tombstone added by the removal is
# shared by all of them:
#
#   eval_default      pa_survival_eval()      metrics = NULL
#   eval_all          pa_survival_eval()      metrics = "all"
#   fit_and_eval_all  pa_fit_and_eval()       metrics = "all"
#   cr_default        pa_survival_eval_cr()   metrics = NULL
#   cr_all            pa_survival_eval_cr()   metrics = "all"
#   two_phase_cc      pa_evaluate_two_phase() case-cohort weights
#   two_phase_ncc     pa_evaluate_two_phase() nested case-control weights
#   two_phase_unit    pa_evaluate_two_phase() unit weights (reference path)
#   sample_design_cc  pa_sample_design()      case-cohort
#
# Between them the scenarios pin all 11 canonical post-removal metric names.
# The three two-phase weightings are pinned separately because they genuinely
# diverge -- test-two-phase.R already asserts cc and ncc differ -- so a single
# weighting would not detect a change confined to one of them.
#
# If this test fails, the removal perturbed a metric it should not have
# touched. Investigate before regenerating the baseline: regenerating it to
# make the test pass destroys the only evidence that the values held.
# The paper branch intentionally replaces CR generalized R-squared weighting.
# Preserve the original CSV: those six CR metrics now use the independent
# original-time event-IPCW moment reference below. All other comparisons still
# use the original removal baseline unchanged, except for ordinary point R2:
# its event-probability direction and inclusive endpoint are now corrected.
# These three values use independent KM risk sets and weighted moments below;
# the historical CSV remains untouched.
# CR C-index now uses the independently enumerated reciprocal-IPCW reference
# below, following the intentional correction of Wolbers equation (3.4).

baseline_rows <- function() {
  path <- testthat::test_path("fixtures", "metric-baseline.csv")
  testthat::skip_if(!file.exists(path), "baseline fixture not present")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# --- scenario reconstruction -------------------------------------------------
# Each builder returns a named numeric vector of metric -> value.

inv_named <- function(metric, value) {
  stats::setNames(round(as.numeric(value), 6), as.character(metric))
}

inv_surv_pred <- function() {
  pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                   new_data = fx_surv(), tau = 10e10)
}

inv_eval <- function(metrics) {
  pred <- inv_surv_pred()
  res <- pa_survival_eval(
    model = fx_cox(), event_time = pred$times,
    predicted_probability = pred$surv_prob, status = pred$status,
    covariates = fx_covs(), new_data = fx_surv(), tau = 10e10,
    metrics = metrics
  )
  inv_named(res$Metric, res$Value)
}

inv_survival_point_reference <- function(p = inv_surv_pred()) {
  t_star <- median(p$times)
  # Reconstruct the restricted outcome's observation mechanism and KM weights
  # from censor counts/risk sets, without the package's censoring/R2 helpers.
  time <- pmin(p$times, t_star)
  observed <- as.numeric(p$status == 1 | p$times > t_star)
  censor_times <- sort(unique(time[observed == 0]))
  factors <- vapply(censor_times, function(t)
    1 - sum(time == t & observed == 0) / sum(time >= t), numeric(1))
  G <- vapply(time, function(t) prod(factors[censor_times < t]), numeric(1))
  q <- observed / G
  w <- q / sum(q)
  k <- PAmeasures:::pa_nearest_time_index(p$times, t_star)
  prediction <- 1 - p$surv_prob[, k]
  y <- as.numeric(p$times <= t_star & p$status == 1)
  ybar <- sum(w * y)
  pbar <- sum(w * prediction)
  vy <- sum(w * (y - ybar)^2)
  vp <- sum(w * (prediction - pbar)^2)
  covariance <- sum(w * (y - ybar) * (prediction - pbar))
  corrected <- ybar + covariance / vp * (prediction - pbar)
  r <- round(covariance^2 / (vy * vp), 4)
  l <- round(sum(w * (y - corrected)^2) / sum(w * (y - prediction)^2), 4)
  c(r2_point = r, l2_point = l, pseudo_r2_point = round(r * l, 4))
}

inv_cr_pred <- function() {
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

inv_eval_cr <- function(metrics) {
  p <- inv_cr_pred()
  res <- pa_survival_eval_cr(
    pred_cif = p$cif_pred[, -1], event_time = p$times,
    time.cif = p$cif_pred[, 1], status = p$status, event_type = 1,
    metrics = metrics
  )
  inv_named(res$Metric, res$Value)
}

inv_cr_event_reference <- function() {
  p <- inv_cr_pred()
  time <- p$times
  status <- p$status
  censor_times <- sort(unique(time[status == 0]))
  # No-covariate Breslow Cox = exp(-Nelson-Aalen), tied times included.
  # Construct hazard jumps directly from original censor counts/risk sets,
  # without the package weight helper.
  jumps <- vapply(censor_times, function(t) sum(time == t & status == 0) /
    sum(time >= t), numeric(1))
  G <- vapply(time, function(t) exp(-sum(jumps[censor_times < t])), numeric(1))
  q <- as.numeric(status > 0) / pmax(G, 1e-4)
  w <- q / sum(q)
  moment <- function(y, prediction) {
    a <- sum(w * y)
    b <- sum(w * prediction)
    vy <- sum(w * (y - a)^2)
    vp <- sum(w * (prediction - b)^2)
    covariance <- sum(w * (y - a) * (prediction - b))
    corrected <- a + covariance / vp * (prediction - b)
    r <- round(covariance^2 / (vy * vp), 4)
    l <- round(sum(w * (y - corrected)^2) / sum(w * (y - prediction)^2), 4)
    c(r = r, l = l, pseudo = round(r * l, 4))
  }
  tau <- max(time)
  grid <- p$cif_pred[, 1]
  cif <- p$cif_pred[, -1]
  k <- which(grid <= tau)
  prediction <- colSums(diff(c(0, grid[k], tau)) * rbind(1, 1 - cif[k, ]))
  interval <- moment(ifelse(status == 1, pmin(time, tau), tau), prediction)
  t_star <- median(time)
  point <- moment(as.numeric(time <= t_star & status == 1),
    cif[max(which(grid <= t_star)), ])
  # C-index uses marginal KM rather than the Cox probabilities above.
  left <- vapply(time, function(t) prod(1 - jumps[censor_times < t]), numeric(1))
  right <- vapply(time, function(t) prod(1 - jumps[censor_times <= t]), numeric(1))
  pairs <- expand.grid(i = which(status == 1 & time <= tau), j = seq_along(time))
  a <- time[pairs$i] < time[pairs$j]
  keep <- a | status[pairs$j] == 2
  pairs <- pairs[keep, ]
  a <- a[keep]
  probability <- ifelse(a, left[pairs$i] * right[pairs$i],
                         left[pairs$i] * left[pairs$j])
  concordance <- sum((prediction[pairs$i] < prediction[pairs$j]) / probability) /
    sum(1 / probability)
  c(pseudo_r2 = unname(interval["pseudo"]), r_square = unname(interval["r"]),
    l_square = unname(interval["l"]), pseudo_r2_point = unname(point["pseudo"]),
    r2_point = unname(point["r"]), l2_point = unname(point["l"]),
    c_index = round(concordance, 4))
}

# Mirrors tp_inputs() in test-two-phase.R. set.seed(2001) is required: both
# weight constructors sample, so the weights are only reproducible under it.
inv_tp_inputs <- function(design) {
  d <- fx_surv()
  pr <- pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                         new_data = d, tau = 10e10)
  set.seed(2001)
  w <- switch(
    design,
    cc   = pa_case_cohort_weights(time = d$time, status = d$status,
                                  subcohort = stats::rbinom(nrow(d), 1, 0.4)),
    ncc  = pa_nested_case_control_weights(time = d$time, status = d$status, m = 2),
    unit = rep(1, nrow(d))
  )
  list(pred = pr, weights = w,
       km_cens = survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1))
}

inv_two_phase <- function(design) {
  inp <- inv_tp_inputs(design)
  res <- pa_evaluate_two_phase(pred_results = inp$pred,
                               km_cens_fit = inp$km_cens,
                               case_weights = inp$weights)
  inv_named(res$Metric, res$Value)
}

inv_sample_design <- function() {
  inp <- inv_tp_inputs("cc")
  res <- pa_sample_design(models = list(cc = inp$pred),
                          case_weights = inp$weights, km_cens = inp$km_cens)
  inv_named(res$Metric, res[["cc"]])
}

inv_builders <- list(
  eval_default     = function() inv_eval(NULL),
  eval_all         = function() inv_eval("all"),
  fit_and_eval_all = function() {
    res <- pa_fit_and_eval(train_data = fx_surv(), covariates = fx_covs(),
                           models = "coxph", metrics = "all")
    num <- vapply(res, is.numeric, logical(1))
    inv_named(names(res)[num], unlist(res[num]))
  },
  cr_default       = function() inv_eval_cr(NULL),
  cr_all           = function() inv_eval_cr("all"),
  two_phase_cc     = function() inv_two_phase("cc"),
  two_phase_ncc    = function() inv_two_phase("ncc"),
  two_phase_unit   = function() inv_two_phase("unit"),
  sample_design_cc = inv_sample_design
)

# --- the gate ----------------------------------------------------------------

for (scenario_name in names(inv_builders)) {
  local({
    scn <- scenario_name
    test_that(paste0(scn, ": surviving metric values are unchanged"), {
      baseline <- baseline_rows()
      expected <- baseline[baseline$Scenario == scn, ]
      expect_gt(nrow(expected), 0)

      observed <- inv_builders[[scn]]()
      cr_reference <- if (scn %in% c("cr_default", "cr_all"))
        inv_cr_event_reference() else numeric()
      point_reference <- if (scn %in% c("eval_default", "eval_all"))
        inv_survival_point_reference() else numeric()
      method_reference <- c(cr_reference, point_reference)

      for (i in seq_len(nrow(expected))) {
        m <- expected$Metric[i]
        expect_true(
          m %in% names(observed),
          info = paste0(scn, ": metric '", m, "' disappeared from the result")
        )
        reference <- if (m %in% names(method_reference)) method_reference[[m]] else expected$Value[i]
        expect_equal(unname(observed[[m]]), unname(reference), tolerance = 1e-6,
          info = paste0(scn, ": metric '", m, "' differs from its reference"))
      }
    })
  })
}

test_that("the baseline covers every canonical metric and no withdrawn one", {
  baseline <- baseline_rows()

  expect_false("r_e" %in% baseline$Metric)
  expect_false("r_sh" %in% baseline$Metric)

  # exactly the 11 canonical names that survive the removal
  expect_identical(sort(unique(baseline$Metric)), sort(c(
    "brier_score", "c_index", "harrell_c", "l2_point", "l_square",
    "pseudo_r2", "pseudo_r2_point", "r2_point", "r_square", "td_auc", "uno_c"
  )))

  # every scenario this file can rebuild is actually present in the fixture
  expect_setequal(unique(baseline$Scenario), names(inv_builders))
})
