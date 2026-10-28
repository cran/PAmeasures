# Integration tests for pa_predict_cif's two optional prediction backends.
# Both live in Suggests, so each test skips when its package is absent -- a
# visible skip, not a silent gap. The full-dependency CI job installs both, so
# these run for real there.

ob_data <- function() {
  d <- fx_cr()
  d$time <- d$obs.times
  d$status <- d$obs.event
  d
}

ob_covs <- function() c("X1", "X2")

test_that("both optional backends are declared in Suggests, not Imports", {
  root <- skip_without_source_tree()
  dcf <- read.dcf(file.path(root, "DESCRIPTION"))

  field <- function(f) {
    if (!f %in% colnames(dcf)) return(character(0))
    trimws(strsplit(dcf[1, f], ",")[[1]])
  }
  imports  <- sub(" .*", "", field("Imports"))
  suggests <- sub(" .*", "", field("Suggests"))

  expect_true("randomForestSRC" %in% suggests)
  expect_true("cmprsk" %in% suggests)
  expect_false("randomForestSRC" %in% imports)
  expect_false("cmprsk" %in% imports)
  # and the newly-required runtime deps really are in Imports
  expect_true(all(c("pec", "magrittr",
                    "purrr", "tibble") %in% imports))
  # rms was removed entirely: it served only the Schemper-Henderson estimator
  # and forced R >= 4.4.0 (findings.md #35, #36)
  expect_false("rms" %in% imports)
  # survminer was removed too: its single use, surv_summary() in Gt(), is now
  # built from the survfit object directly. Verified bitwise identical.
  expect_false("survminer" %in% imports)
  expect_false("survminer" %in% suggests)
})

test_that("neither optional backend is imported into the namespace", {
  # They are reached through their own predict() methods on the fitted object,
  # so no importFrom is needed -- and importing randomForestSRC previously made
  # the package impossible to install (findings.md #6).
  imported <- unlist(getNamespaceImports("PAmeasures"), use.names = FALSE)

  expect_false("predict.rfsrc" %in% imported)
  expect_false("crr" %in% imported)
})

test_that("fg_model backend computes CIF predictions via cmprsk", {
  skip_if_not_installed("cmprsk")
  dd <- ob_data()
  X <- as.matrix(dd[, ob_covs(), drop = FALSE])

  fg <- cmprsk::crr(ftime = dd$time, fstatus = dd$status, cov1 = X,
                    failcode = 1, cencode = 0)

  res <- pa_predict_cif(fg_model = fg, newdata = dd, covs = ob_covs(),
                        event.type = 1, tau = max(dd$time))

  expect_type(res, "list")
  expect_true(all(c("times", "status", "cif_pred", "pred", "linear.pred")
                  %in% names(res)))
  expect_identical(length(res$pred), 200L)
  expect_true(all(is.finite(res$pred)))
  # same orientation contract as the coxph path: time grid in column 1
  expect_identical(ncol(res$cif_pred), 201L)
  expect_true(all(diff(res$cif_pred[, 1]) > 0))
  expect_snapshot_value(snap_num(head(res$pred, 10)), style = "serialize")
})

test_that("fg_model output feeds pa_summarize_cr like the coxph path", {
  skip_if_not_installed("cmprsk")
  dd <- ob_data()
  X <- as.matrix(dd[, ob_covs(), drop = FALSE])
  fg <- cmprsk::crr(ftime = dd$time, fstatus = dd$status, cov1 = X,
                    failcode = 1, cencode = 0)

  pred <- pa_predict_cif(fg_model = fg, newdata = dd, covs = ob_covs(),
                         event.type = 1, tau = max(dd$time))
  res <- pa_summarize_cr(list(fg = pred), event_type = 1)

  expect_metric_table(res)
  expect_identical(names(res), c("Metric", "fg"))
  expect_true("c_index" %in% res$Metric)
})

test_that("cr_model backend computes CIF predictions via randomForestSRC", {
  # Skips visibly when randomForestSRC is absent. This is the path that the
  # stale NAMESPACE import gestured at but never actually exercised.
  skip_if_not_installed("randomForestSRC")
  dd <- ob_data()

  rf <- randomForestSRC::rfsrc(
    Surv(time, status) ~ X1 + X2, data = dd, ntree = 50, seed = -1003
  )

  res <- pa_predict_cif(cr_model = rf, newdata = dd, covs = ob_covs(),
                        event.type = 1, tau = max(dd$time))

  expect_type(res, "list")
  expect_true(all(c("times", "status", "cif_pred", "pred") %in% names(res)))
  expect_identical(length(res$pred), 200L)
  expect_true(all(is.finite(res$pred)))
})

test_that("fg_model errors informatively when cmprsk is unavailable", {
  # Exercises the requireNamespace guard by making the package look absent.
  # The guard must name both the argument and the package to install.
  dd <- ob_data()
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "cmprsk")) FALSE else TRUE
    },
    .package = "base"
  )

  expect_error(
    pa_predict_cif(fg_model = structure(list(), class = "crr"),
                   newdata = dd, covs = ob_covs(), event.type = 1),
    "requires the 'cmprsk' package"
  )
})

test_that("cr_model errors informatively when randomForestSRC is unavailable", {
  dd <- ob_data()
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "randomForestSRC")) FALSE else TRUE
    },
    .package = "base"
  )

  expect_error(
    pa_predict_cif(cr_model = structure(list(), class = "rfsrc"),
                   newdata = dd, covs = ob_covs(), event.type = 1),
    "requires the 'randomForestSRC' package"
  )
})
