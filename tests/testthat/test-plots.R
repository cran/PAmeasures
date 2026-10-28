plot_input <- function() {
  pa_predict_coxph(model = fx_cox(), covs = fx_covs(),
                       new_data = fx_surv(), tau = 10e10)
}

test_that("pa_plot_pred requires linear.pred, pred, times and status", {
  pred <- plot_input()

  expect_true(all(c("linear.pred", "pred", "times", "status") %in% names(pred)))
  expect_error(
    pa_plot_pred(pred[c("pred", "times", "status")]),
    "Input data must contain columns: linear.pred, pred, times, status"
  )
})

test_that("pa_plot_pred with the default sample_index draws every subject (FINDING 18 fixed)", {
  # sample_index defaults to NULL, which previously reached x_var[NULL] and
  # produced a zero-length data frame -- a silently blank chart on the
  # documented default call. NULL now means "plot every subject".
  p <- pa_plot_pred(plot_input())
  built <- ggplot2::ggplot_build(p)

  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), 200L)
  expect_identical(nrow(built$data[[1]]), 200L)
  expect_identical(nrow(built$data[[2]]), 200L)
})

test_that("pa_plot_pred sample_index thins the plot to the chosen subjects", {
  built <- ggplot2::ggplot_build(pa_plot_pred(plot_input(), sample_index = 1:50))

  expect_identical(nrow(built$data[[1]]), 50L)
})

test_that("pa_plot_pred documents the arguments it actually has (FINDING 20)", {
  # man/plot_pred.Rd previously documented a sample_size argument that did not
  # exist, while the real sample_index -- whose NULL default blanked the plot --
  # was undocumented.
  rd <- paste(readLines(skip_without_rd("pa_plot_pred.Rd"), warn = FALSE),
              collapse = "\n")

  expect_true(grepl("item{sample_index}", rd, fixed = TRUE))
  expect_true(grepl("item{restrict_time}", rd, fixed = TRUE))
  expect_false(grepl("sample_size", rd, fixed = TRUE))
})

test_that("pa_plot_pred builds exactly two layers", {
  p <- pa_plot_pred(plot_input())
  built <- ggplot2::ggplot_build(p)

  expect_identical(length(p$layers), 2L)
  expect_snapshot_value(
    vapply(built$data, nrow, integer(1)), style = "serialize"
  )
})

test_that("pa_plot_pred uses the supplied labels", {
  p <- pa_plot_pred(plot_input(), title = "characterization", xlab = "RS", ylab = "Days")

  expect_identical(p$labels$title, "characterization")
  expect_identical(p$labels$x, "RS")
  expect_identical(p$labels$y, "Days")
})

test_that("pa_plot_pred defaults xlab to Risk Score and ylab to Days", {
  p <- pa_plot_pred(plot_input())

  expect_identical(p$labels$x, "Risk Score")
  expect_identical(p$labels$y, "Days")
})

test_that("pa_plot_pred restrict_time caps the plotted times", {
  cap <- 1
  p <- pa_plot_pred(plot_input(), restrict_time = cap)

  expect_true(all(p$data$times <= cap))
  expect_gt(sum(p$data$times == cap), 0)
})

test_that("pa_plot_summary combines one panel per input", {
  local_null_device()
  pred <- plot_input()

  p <- pa_plot_summary(list(a = pred, b = pred), ncol = 2)

  expect_s3_class(p, "patchwork")
  expect_identical(length(p$patches$plots) + 1L, 2L)
  expect_no_error(suppressWarnings(print(p)))
})

test_that("pa_plot_summary accepts explicit panel titles", {
  local_null_device()
  pred <- plot_input()

  p <- pa_plot_summary(list(a = pred, b = pred),
                         titles = c("first", "second"), ncol = 2)

  expect_s3_class(p, "patchwork")
  expect_no_error(suppressWarnings(print(p)))
})
