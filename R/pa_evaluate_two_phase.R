
# Helper (internal)
#' @keywords internal
#' @noRd
km_surv <- function(t, km_cens) {
  if (is.na(t)) return(NA_real_)
  summary(km_cens, times = t, extend = TRUE)$surv
}

#' Evaluate Survival Model Metrics for Two-Phase Sampling Designs
#'
#' This function constructs a standardized evaluation dataset (`dat1`) from the
#' analysis dataset, predicted survival probabilities, and sampling weights, and
#' computes key performance metrics in a tidy table format.
#' Works for common two-phase sampling designs:
#'   - Case-cohort (unstratified or stratified)
#'   - Nested case-control (matched or unmatched)
#'
#' @param case_weights Numeric sampling weights in the same subject order as
#'   \code{pred_results}. Prediction helpers sort subjects by observed time;
#'   reorder weights using \code{order(new_data$time)} as well.
#' @param pred_results A prediction list, such as the output of
#'   \code{pa_predict_coxph()}, with \code{times} (or legacy \code{time}),
#'   \code{status} (0 = censored, 1 = event), \code{pred} (predicted mean
#'   survival times), and \code{surv_prob} (subjects by sorted observed times).
#' @param t_star (Optional) Evaluation time for Brier score, time-dependent
#'   AUC, and calibration.
#'   If NULL, it defaults to the median survival time.
#' @param tau (Optional) A numeric value specifying the truncation time for calculating explained variation metrics. 
#'   If NULL, it defaults to the maximum observed survival time.
#' @param km_cens_fit Kaplan-Meier estimate of the censoring distribution, a
#'   `survival::survfit` object. Not required when only calibration is requested.
#' @param metrics Character vector of metrics. \code{NULL} uses the default
#'   metrics, including ICI; \code{"all"} additionally returns \code{r_square},
#'   \code{l_square}, \code{e50}, and \code{e90}.
#' @param calibration_args Optional named list of calibration settings,
#'   e.g. \code{list(nk = 3)}. Sampling weights are supplied by
#'   \code{case_weights}, not through this list.
#'
#' @details Calibration uses the event risk from the last prediction-grid
#'   column at or before \code{t_star}, with the exact requested outcome
#'   horizon. Both the Cox spline calibration model and absolute-error
#'   summaries use sampling weights; see [pa_calibration()].
#'   Calibration-only inputs need \code{times}, \code{status}, and
#'   \code{surv_prob}, without fitted models or mean time predictions.
#'
#' @return A tibble with columns:
#'   - `Metric`: name of metric
#'   - `Value`: numeric value
#'
#'
#' @export
#' @importFrom dplyr mutate case_when
#' @importFrom tibble tibble
#' @importFrom tidyr unnest
#' @examples
#' d <- pa_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' d <- d[, c("time", "status", "x1", "x2")]
#' fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
#'                        data = d, x = TRUE, y = TRUE)
#' pred <- pa_predict_coxph(model = fit, covs = c("x1", "x2"), new_data = d)
#'
#' # case-cohort sampling weights, and the censoring distribution they need
#' set.seed(2025)
#' in_subcohort <- stats::rbinom(nrow(d), 1, 0.4) == 1
#' w <- pa_case_cohort_weights(d$time, d$status, subcohort = in_subcohort)
#' keep <- which(in_subcohort | d$status == 1)
#' keep <- keep[order(d$time[keep])]
#' pred <- pa_predict_coxph(fit, c("x1", "x2"), new_data = d[keep, ])
#' w <- w[keep]
#' km <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
#'
#' pa_evaluate_two_phase(pred_results = pred, km_cens_fit = km,
#'                       case_weights = w)
pa_evaluate_two_phase <- function(pred_results,
                                                 t_star = NULL, tau = 10e10, 
                                                 km_cens_fit, case_weights, 
                                                 metrics = c("pseudo_r2", "harrell_c", "uno_c", "brier_score", "td_auc", "ici"),
                                                 calibration_args = list()) {
  metrics <- pa_normalize_metrics(metrics)
  default_metrics <- c("pseudo_r2", "harrell_c", "uno_c", "brier_score", "td_auc", "ici")
  valid_metrics <- c("pseudo_r2", "r_square", "l_square", "harrell_c",
                     "uno_c", "brier_score", "td_auc", "ici", "e50", "e90")
  if (is.null(metrics)) metrics <- default_metrics
  invalid <- setdiff(metrics, c("all", valid_metrics))
  if (length(invalid) > 0 || length(metrics) == 0) {
    stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  if ("all" %in% metrics) metrics <- valid_metrics
  # Accept the canonical prediction field without partial `$time` matching.
  if (is.null(pred_results[["time"]])) pred_results$time <- pred_results[["times"]]
  
  if (is.null(t_star)) t_star <- quantile(pred_results$time, 0.5)
  if (is.null(tau)) tau <- max(pred_results$time)

  calibration_metrics <- intersect(c("ici", "e50", "e90"), metrics)
  calibration_results <- list()
  if (length(calibration_metrics)) {
    pa_validate_calibration_weights(case_weights, length(pred_results$time))
    if (!is.list(calibration_args) || (length(calibration_args) &&
        (is.null(names(calibration_args)) || any(names(calibration_args) == "") ||
         anyDuplicated(names(calibration_args)))) ||
        any(!names(calibration_args) %in% "nk")) {
      stop("'calibration_args' must be a named list containing only 'nk'; ",
           "supply sampling weights in 'case_weights'.", call. = FALSE)
    }
    p <- pa_calibration_grid_risk(pred_results$surv_prob,
      sort(pred_results$time), t_star, length(pred_results$time), survival = TRUE)
    calibration <- pa_calibration_values(pa_calibration, c(list(
      p = p, event_time = pred_results$time, status = pred_results$status,
      t_star = t_star, case_weights = case_weights), calibration_args))
    calibration_results <- as.list(round(calibration[calibration_metrics], 4))
    if (all(metrics %in% calibration_metrics)) {
      return(data.frame(Metric = names(calibration_results),
        Value = unlist(calibration_results, use.names = FALSE),
        stringsAsFactors = FALSE))
    }
  }
  
  dat1_intermediate <- tibble(time = pred_results$time, 
                                  status = pred_results$status,
                                  pred = pred_results$surv_prob[, max(which(pred_results$time <= t_star))],
                                  pred.t = pred_results$pred,
                                  case_weights = case_weights)
  dat1_intermediate <- dat1_intermediate %>%
    mutate(
      surv_obj = Surv(pred_results$time, 
                      pred_results$status),
      # Determine the time for IPCW calculation based on event status and time
      weight_time = case_when(
        status == 1 & time <= t_star ~ pmax(time, 0), # Category 1: Events before t_star
        time > t_star              ~ t_star,    # Category 2: Censored after t_star
        TRUE                          ~ NA_real_      # Category 3: Events after t_star (no weight)
      )
    )
  
  # Now, create the nested '.pred' column using the columns from the step above.
  dat1 <- dat1_intermediate %>%
    mutate(
      .pred = purrr::map2(
        pred, weight_time,
        ~ tibble(
          .eval_time = t_star,
          .pred_survival = .x,
          # Calculate Inverse Probability of Censoring Weight (IPCW)
          .weight_censored = if (!is.na(.y)) 1 / km_surv(.y, km_cens_fit) else 0
        )
      )
    )
    results <- list()
    
    if (any(c("r_square", "l_square", "pseudo_r2") %in% metrics)) {
      r_l_list <- pam.r2_metrics(pred_results$pred, 
                                 pred_results$time, 
                                 pred_results$status,
                                 tau = tau,
                                 case_weight = case_weights)
      
      if ("r_square" %in% metrics) results <- append(results, list("r_square" = round(as.numeric(r_l_list$R_square), 2)))
      if ("l_square" %in% metrics) results <- append(results, list("l_square" = round(as.numeric(r_l_list$L_square), 2)))
      if ("pseudo_r2" %in% metrics) results <- append(results, list("pseudo_r2" = round(as.numeric(r_l_list$Pseudo_R_squared), 2)))
    }
    
    if ("harrell_c" %in% metrics) {
      c_index <- survival::concordance(
        dat1$surv_obj ~ dat1$pred.t,
        weights = dat1$case_weights,
      )$concordance
      results <- append(results, list("harrell_c" = round(c_index, 4)))
    }
    
    if ("uno_c" %in% metrics) {
      c_index <- survival::concordance(
        dat1$surv_obj ~ dat1$pred.t,
        weights = dat1$case_weights,
        timewt = "n/G2",
      )$concordance
      results <- append(results, list("uno_c" = round(c_index,4 )))
    }
    
    # --- Brier score
    if ("brier_score" %in% metrics) {
      brier <- yardstick::brier_survival(
        data         = dat1,
        truth        = surv_obj,
        .pred,
        case_weights = case_weights
      )$.estimate
      results <- append(results, list("brier_score" = round(brier, 4)))
    }
    
    # --- Time-dependent AUC
    if ("td_auc" %in% metrics) {
      # yardstick uses tidyr::unnest here but only suggests tidyr itself.
      # Import that backend dependency because AUC is in our default output.
      auc <- yardstick::roc_auc_survival(
        data         = dat1,
        truth        = surv_obj,
        .pred,
        case_weights = case_weights
      )$.estimate
      results <- append(results, list("td_auc" = round(auc, 4)))
    }
    

    results <- c(results, calibration_results)
    result_df <- data.frame(
      Metric = names(results),
      Value = unlist(results, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    
    result_df
    
}


#' Summarize multiple two-phase (case-cohort) models into a wide table
#'
#' @description
#' For a list of models (each providing `time`, `status`, `surv_prob`, `pred`),
#' calls `pa_evaluate_two_phase()` for each and returns a wide
#' comparison table with metrics as rows and model names as columns.
#'
#' @param models A **named list** in which each element corresponds to a fitted model.  
#'   Each model entry must itself be a list containing:
#'   \itemize{
#'     \item \code{times} -- numeric vector of observed follow-up times.
#'     \item \code{surv_prob} -- an \eqn{n \times K} matrix (or data frame) of
#'           subject-specific predicted survival probabilities on a common time grid.
#'     \item \code{status} -- event indicator (1 = event, 0 = censored).
#'     \item \code{pred} -- (optional) predicted mean survival time (restricted or unrestricted).
#'     \item \code{new_data} -- (optional) dataset used for prediction.
#'     \item \code{covs} -- character vector of covariate names used for prediction.
#'     \item \code{model} -- (optional) the underlying fitted survival model object.
#'   }
#'   Note: If \code{pred} is provided, it will be used to calculate R2 and concordance measure.
#' @param case_weights Numeric vector of inverse-inclusion sampling weights,
#'   aligned with each model's evaluation subjects.
#' @param km_cens A Kaplan-Meier fit for censoring (based on training data) used to compute IPCW.
#' @param metrics Character vector of metrics to compute (default shown below).
#' @param t_star Optional numeric scalar, specify the time point to evaluate Brier score and AUC.Default is median of observation time.
#' @param tau Optional numeric scalar, specify the max time horizon for R2 measure and concordance measure (default = 10e10).
#' @param digits Integer number of decimal places to round values (default 2).
#' @param calibration_args Optional named list passed to [pa_calibration()],
#'   e.g. \code{list(nk = 3)}. Weights come from \code{case_weights}.
#'
#' @return A data.frame: rows = metrics, columns = model names.
#' @export
#' @examples
#' d <- pa_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' d <- d[, c("time", "status", "x1", "x2")]
#' fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
#'                        data = d, x = TRUE, y = TRUE)
#' pred <- pa_predict_coxph(model = fit, covs = c("x1", "x2"), new_data = d)
#'
#' # case-cohort sampling weights, and the censoring distribution they need
#' set.seed(2025)
#' in_subcohort <- stats::rbinom(nrow(d), 1, 0.4) == 1
#' w <- pa_case_cohort_weights(d$time, d$status, subcohort = in_subcohort)
#' keep <- which(in_subcohort | d$status == 1)
#' keep <- keep[order(d$time[keep])]
#' pred <- pa_predict_coxph(fit, c("x1", "x2"), new_data = d[keep, ])
#' w <- w[keep]
#' km <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
#'
#' # the wide summary of the same two-phase evaluation
#' pa_sample_design(models = list(cc = pred), case_weights = w, km_cens = km)
pa_sample_design <- function(models,
                              case_weights,
                              km_cens,
                              metrics = c("pseudo_r2", "harrell_c", "uno_c",
                                          "brier_score", "td_auc", "ici"),
                              t_star = NULL, tau = NULL,
                              digits = 2, calibration_args = list()) {
  metrics <- pa_normalize_metrics(metrics)
  
  if (!is.list(models) || length(models) == 0)
    stop("'models' must be a non-empty named list.")
  if (is.null(names(models)) || any(names(models) == ""))
    names(models) <- paste0("Model_", seq_along(models))
  
  required <- c("times", "status", "surv_prob")
  needs_time_prediction <- is.null(metrics) || any(metrics %in% c(
    "all", "pseudo_r2", "r_square", "l_square", "harrell_c", "uno_c"))
  if (needs_time_prediction) required <- c(required, "pred")
  
  eval_one <- function(mod, name) {
    if (!all(required %in% names(mod))) {
      stop(sprintf("Model '%s' must contain: %s",
                   name, paste(required, collapse = ", ")))
    }
    
    # Build the single-model input structure expected by your evaluator
    pred_results <- list(
      time      = mod$times,
      status    = mod$status,
      surv_prob = mod$surv_prob,
      pred      = mod$pred
    )
    
    res <- pa_evaluate_two_phase(
      pred_results = pred_results,
      t_star       = t_star,    # global (NULL allowed)
      tau          = tau,       # global (NULL allowed)
      km_cens_fit  = km_cens,
      case_weights = case_weights,
      metrics      = metrics,
      calibration_args = calibration_args
    )
    
    if (!all(c("Metric", "Value") %in% names(res)))
      stop(sprintf("Unexpected result from pa_evaluate_two_phase() for '%s'", name))
    
    # Normalize minor label variant
    res$Metric <- gsub("^Time Dependent AUC$", "td_auc", res$Metric)
    
    res$Model <- name
    res[, c("Model", "Metric", "Value")]
  }
  
  # Evaluate all models
  pieces <- mapply(eval_one, models, names(models), SIMPLIFY = FALSE)
  long   <- do.call(rbind, pieces)
  rownames(long) <- NULL
  
  # Pivot to wide (rows = Metric, cols = Model)
  wide <- reshape(
    long,
    idvar = "Metric",
    timevar = "Model",
    direction = "wide"
  )
  names(wide)    <- sub("^Value\\.", "", names(wide))
  rownames(wide) <- NULL
  wide           <- wide[, c("Metric", setdiff(names(wide), "Metric")), drop = FALSE]
  
  # Preferred ordering (keep present ones)
  preferred <- c("pseudo_r2", "r_square", "l_square",
                 "harrell_c", "uno_c",
                 "brier_score", "td_auc", "ici", "e50", "e90")
  present <- intersect(preferred, wide$Metric)
  others  <- setdiff(wide$Metric, preferred)
  wide$Metric <- factor(wide$Metric, levels = c(present, sort(others)))
  wide <- wide[order(wide$Metric), ]
  wide$Metric <- as.character(wide$Metric)
  
  # Round numeric columns
  num_cols <- setdiff(names(wide), "Metric")
  wide[num_cols] <- lapply(wide[num_cols], function(x) if (is.numeric(x)) round(x, digits) else x)
  
  wide
}
