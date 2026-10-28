# Reduce values to something small and platform-stable before snapshotting.
# Never snapshot a fitted model or a full prediction matrix.

snap_num <- function(x, digits = 6) round(as.numeric(x), digits)

mat_fingerprint <- function(m, digits = 6) {
  m <- as.matrix(m)
  list(
    dim  = dim(m),
    min  = round(min(m, na.rm = TRUE), digits),
    max  = round(max(m, na.rm = TRUE), digits),
    mean = round(mean(m, na.rm = TRUE), digits),
    n_na = sum(is.na(m))
  )
}

# Structural assertions common to every Metric/Value result table.
expect_metric_table <- function(res) {
  testthat::expect_true(is.data.frame(res))
  testthat::expect_true("Metric" %in% names(res))
  testthat::expect_gt(nrow(res), 0)
  testthat::expect_type(res$Metric, "character")
  testthat::expect_false(any(duplicated(res$Metric)))
  invisible(res)
}

# Path to the package SOURCE root, or NULL when tests run against an installed
# package (as under covr or R CMD check), where no source tree is present.
# Tests that read DESCRIPTION/NAMESPACE or re-load the source must skip then.
pkg_source_root <- function() {
  root <- tryCatch(
    normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
    error = function(e) NULL
  )
  if (is.null(root)) return(NULL)
  if (!file.exists(file.path(root, "DESCRIPTION"))) return(NULL)
  if (!dir.exists(file.path(root, "R"))) return(NULL)
  root
}

skip_without_source_tree <- function() {
  root <- pkg_source_root()
  testthat::skip_if(is.null(root), "no package source tree (installed-package run)")
  root
}

# The clean-subprocess reproductions re-load the package source in a fresh R
# session. Under covr the source is instrumented in a temporary library, so a
# plain pkgload::load_all() there fails for reasons unrelated to what the test
# asserts. Skip them during coverage runs only; they run normally otherwise.
skip_if_covr <- function() {
  testthat::skip_if(
    nzchar(Sys.getenv("R_COVR")),
    "clean-subprocess reload is incompatible with covr instrumentation"
  )
}

# Path to an .Rd file in the package source, skipping when man/ is absent.
# Under covr and R CMD check the package is installed to a temporary library
# where Rd sources are replaced by a help database, so man/ does not exist.
skip_without_rd <- function(name) {
  root <- skip_without_source_tree()
  rd <- file.path(root, "man", name)
  testthat::skip_if(!file.exists(rd), paste0("man/", name, " not present"))
  rd
}


# Printing a ggplot/patchwork object in a non-interactive session with no
# device open makes R open the default device, which writes an Rplots.pdf into
# the working directory and leaves it behind. Tests that print a plot call this
# first: it opens a pdf device pointed at the null file, scoped to the calling
# test, so rendering is still exercised but nothing is written to disk.
local_null_device <- function(.local_envir = parent.frame()) {
  withr::local_pdf(nullfile(), .local_envir = .local_envir)
}
