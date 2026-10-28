test_that("FHS figure consumers resolve new fits and retain historical runs", {
  root <- skip_without_source_tree()
  path <- file.path(root, "figures", "original_format.R")
  skip_if_not(file.exists(path), "figure helpers are excluded from installed builds")
  helper <- new.env(parent = globalenv())
  sys.source(path, envir = helper)
  out <- tempfile("fhs_fit_")
  dir.create(file.path(out, "framingham"), recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_error(helper$paper_framingham_fit(out), "Missing FHS training fit")
  old <- file.path(out, "framingham", "framingham_visit3_10y.rds")
  current <- file.path(out, "framingham", "framingham_visit3_12y.rds")
  saveRDS(list(horizon = 10), old)
  expect_identical(helper$paper_framingham_fit(out), old)
  saveRDS(list(horizon = 12), current)
  expect_identical(helper$paper_framingham_fit(out), current)
  expect_identical(readRDS(old), list(horizon = 10))

  # Exercise the active-horizon checksum against the new fit, while avoiding
  # rendering: an old hard-coded 10y path would reject this valid selection.
  selected <- file.path(out, "framingham", "selected.rds")
  saveRDS(list(settings = list(tau = 4380), predictions = list(model = list())), selected)
  utils::write.csv(data.frame(Horizon_years = 12,
    Fit_source_md5 = unname(tools::md5sum(current)),
    Source_predictions = normalizePath(selected)),
    file.path(out, "framingham", "active_horizon.csv"), row.names = FALSE)
  helper$paper_original_real_plot <- function(dataset, predictions, tau) {
    expect_equal(tau, 4380)
    list(list(data = data.frame(times = 1)))
  }
  helper$paper_original_save <- function(plot, dir, stem, width, height) {
    for (ext in c("pdf", "png")) writeLines("test", file.path(dir, paste0(stem, ".", ext)))
  }
  record <- helper$paper_original_real(root, out, datasets = "framingham")
  expect_identical(record$Source_predictions, normalizePath(selected))
  saveRDS(list(horizon = 12, changed = TRUE), current)
  expect_error(helper$paper_original_real(root, out, datasets = "framingham"),
               "FHS fit changed")
})
