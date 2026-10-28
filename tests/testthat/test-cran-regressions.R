test_that("missing evaluation rows produce a suppressible warning and consistent metrics", {
  train <- fx_surv()
  validation <- train
  validation$x1[1] <- NA_real_
  evaluate <- function(d) pa_fit_and_eval(train, fx_covs(), predicted_data = d)
  expect_warning(result <- evaluate(validation), "Dropped 1 incomplete row.*predicted_data")
  expect_equal(result, evaluate(validation[-1, ]))
  expect_output(suppressMessages(suppressWarnings(evaluate(validation))), NA)
})

test_that("missing training rows use the same complete sample for fitting and evaluation", {
  d <- fx_surv()
  d$x2[2] <- NA_real_
  expect_warning(result <- pa_fit_and_eval(d, fx_covs()),
                 "Dropped 1 incomplete row.*train_data")
  expect_equal(result, pa_fit_and_eval(d[-2, ], fx_covs()))
  expect_output(suppressWarnings(pa_fit_and_eval(d, fx_covs())), NA)
  d$x1 <- NA_real_
  expect_error(pa_fit_and_eval(d, fx_covs()), "No complete rows remain")
})

test_that("WLS failure is reported as a warning without console output", {
  fit <- fx_survreg()
  d <- fx_surv()
  local_mocked_bindings(lm = function(...) stop("fixture WLS failure"),
                       .package = "PAmeasures")
  evaluate <- function() pa_predict_survreg(fit, fx_covs(), new_data = d,
                                            predict = FALSE)
  expect_warning(evaluate(), "Error in WLS fitting: fixture WLS failure")
  expect_output(suppressWarnings(evaluate()), NA)
})

test_that("CIF fit failures use suppressible warning conditions", {
  fit <- fx_cox()
  d <- fx_surv()
  local_mocked_bindings(survfit = function(...) stop("fixture CIF failure"),
                       .package = "PAmeasures")
  evaluate <- function() PAmeasures:::get_CIF(fit, list(fit), d)
  expect_warning(result <- evaluate(), "fixture CIF failure")
  expect_null(result)
  expect_output(suppressWarnings(evaluate()), NA)
})

test_that("default Fine-Gray simulation respects the caller's random stream", {
  simulate <- function(...) pa_simulate_fine_gray(
    80, 2, c(.5, -.5), c(-.3, .3), censor = .3, ...)
  withr::local_seed(721)
  before <- .Random.seed
  first <- simulate()
  expect_false(identical(before, .Random.seed))
  second <- simulate()
  expect_false(identical(first, second))
  set.seed(721)
  expect_identical(simulate(), first)
  expect_identical(simulate(seed = 213), simulate(seed = 213))
})

test_that("plot sampling respects the caller's stream and explicit seeds", {
  pred <- pa_predict_coxph(fx_cox(), fx_covs())
  selected <- NULL
  capture_plot <- function(data, sample_index, ...) {
    selected <<- sample_index
    ggplot2::ggplot()
  }
  draw <- function(...) pa_plot_summary(list(cox = pred),
                                        plot_fun = capture_plot, ...)
  withr::local_seed(123)
  expected <- sample(seq_along(pred$times), 30)
  expected_state <- .Random.seed
  set.seed(123)
  invisible(draw(sample_size = 30))
  expect_identical(selected, expected)
  expect_identical(.Random.seed, expected_state)
  invisible(draw(sample_size = 30, seed = 812))
  explicit <- selected
  invisible(draw(sample_size = 30, seed = 812))
  expect_identical(selected, explicit)
  before <- .Random.seed
  invisible(draw())
  expect_identical(.Random.seed, before)
})

test_that("two-phase default AUC has a required runtime dependency", {
  desc <- utils::packageDescription("PAmeasures")
  imports <- trimws(strsplit(desc$Imports, ",", fixed = TRUE)[[1]])
  expect_true("tidyr" %in% sub("[[:space:]]*\\(.*", "", imports))
  d <- fx_surv()
  p <- pa_predict_coxph(fx_cox(), fx_covs(), new_data = d)
  km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, d)
  result <- pa_evaluate_two_phase(p, km_cens_fit = km, case_weights = rep(1,nrow(d)))
  auc <- result$Value[result$Metric == "td_auc"]
  expect_length(auc, 1L)
  expect_true(is.finite(auc) && auc >= 0 && auc <= 1)
})
