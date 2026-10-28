framingham_analysis <- function() {
  root <- skip_without_source_tree()
  path <- file.path(root, "analysis", "framingham.R")
  testthat::skip_if_not(file.exists(path), "paper analysis is excluded from installed builds")
  e <- new.env(parent = globalenv())
  sys.source(path, envir = e)
  e
}

third_exam_fixture <- function() {
  baseline <- data.frame(
    RANDID = 1:8, PERIOD = 1, TIME = 0, PREVHYP = c(0, 1, 0, 0, 0, 0, 0, 0),
    DEATH = c(0, 0, 0, 0, 1, 0, 0, 0),
    TIMEDTH = c(8766, 8766, 8766, 8766, 6500, 8766, 8766, 8766),
    HYPERTEN = c(1, 1, 1, 1, 1, 0, 0, 0),
    TIMEHYP = c(6000, 0, 4000, 4400, 6500, 8766, 8766, 8766))
  covs <- c("SEX", "TOTCHOL", "AGE", "SYSBP", "DIABP", "CURSMOKE",
            "CIGPDAY", "BMI", "DIABETES", "HEARTRTE", "GLUCOSE", "educ")
  for (v in covs) baseline[[v]] <- seq_len(nrow(baseline))
  third <- baseline[1:7, ]
  third$PERIOD <- 3
  third$TIME <- c(rep(4400, 6), 4700)
  third$PREVHYP <- c(0, 1, 1, 1, 0, 0, 0)
  third$GLUCOSE[6] <- NA
  baseline$GLUCOSE[7] <- NA # A missing baseline predictor must not exclude id 7.
  rbind(baseline, third)
}

test_that("third-examination cohort uses current eligibility and actual examination dates", {
  e <- framingham_analysis()
  p <- e$paper_framingham_cohort(third_exam_fixture())
  # Excludes baseline disease, prior/concurrent disease, missing current data,
  # and absent third visits; allows missing baseline predictors.
  expect_identical(p$cohort$RANDID, c(1L, 5L, 7L))
  expect_equal(p$cohort$time, c(1600, 2100, 4066))
  expect_equal(p$cohort$status, c(1, 2, 0)) # Death wins the endpoint tie for id 5.
  expect_equal(p$cohort$landmark_time, c(4400, 4400, 4700))
  expect_equal(p$flow$N, c(8, 7, 7, 6, 4, 3))
  expect_identical(p$reduced_covs, c("SEX", "TOTCHOL", "AGE", "DIABP"))
})

test_that("moving future endpoints does not select another examination or predictors", {
  e <- framingham_analysis()
  d <- third_exam_fixture()
  before <- e$paper_framingham_cohort(d)$cohort
  selected <- d$RANDID %in% before$RANDID
  d$TIMEHYP[selected] <- d$TIMEHYP[selected] + 500
  d$TIMEDTH[selected] <- d$TIMEDTH[selected] + 500
  # Additional later visits must never replace examination 3.
  fourth <- d[d$PERIOD == 3 & selected, ]
  fourth$PERIOD <- 4
  fourth$TIME <- 5500
  fourth$GLUCOSE <- 999
  after <- e$paper_framingham_cohort(rbind(d, fourth))$cohort
  keep <- setdiff(names(before), c("time", "endpoint_time"))
  expect_equal(after[, keep], before[, keep])
  expect_equal(after$time, before$time + 500)
})
