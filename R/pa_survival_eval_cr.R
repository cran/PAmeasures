#' Evaluate Survival Prediction Performance with Competing Risks
#'
#' This function evaluates the predictive performance of a survival model in the presence of competing risks.
#' It calculates various performance metrics including Pseudo R-squared, R-squared, L-squared, C-index, 
#' Brier Score, and Time-Dependent AUC.
#'
#' @param pred_cif Numeric matrix of cumulative incidence functions with one
#'   row per value of \code{time.cif} and one column per subject.
#' @param event_time A numeric vector of observed event or censoring times for subjects.
#' @param time.cif A numeric vector of time points corresponding to the predicted CIFs.
#' @param status A numeric vector indicating the event type for each subject: by default
#'   \itemize{
#'     \item `0` for censored cases.
#'     \item `1` for the primary event of interest.
#'     \item `2` (or other values) for competing events.
#'   }
#' @param metrics A character vector specifying the performance metrics to calculate. 
#'   The default \code{NULL} computes pseudo R-squared (integrated and point),
#'   C-index, Brier score, time-dependent AUC, and ICI. Use \code{"all"} to
#'   include the components and E50/E90 as well.
#'   Valid options include:
#'   \itemize{
#'     \item `"pseudo_r2"`
#'     \item `"r_square"`
#'     \item `"l_square"`
#'     \item `"pseudo_r2_point"`
#'     \item `"r2_point"`
#'     \item `"l2_point"`
#'     \item `"c_index"`
#'     \item `"brier_score"`
#'     \item `"td_auc"`
#'     \item `"ici"`, `"e50"`, `"e90"`
#'   }
#' @param t_star Optional numeric value specifying a reference time for evaluating point-version Pseudo_R_square, 
#'   the Brier score, time-dependent AUC, and calibration. Defaults to the median
#'   of `event_time`.
#' @param tau Optional numeric value specifying a restricted time horizon for model evaluation. 
#'   Defaults to the maximum observed follow-up time.
#' @param event_type An integer indicating the primary event type of interest (default = `1`).
#' @param censoring Censoring estimator for the generalized R-squared metrics:
#'   `"cox"` (default, optionally conditional on `censor_covariates`) or `"km"`.
#'   Other metrics retain their own censoring adjustments.
#' @param censor_covariates Optional numeric vector, matrix or data frame of
#'   censoring-model covariates, in the same subject order as `event_time` and
#'   the columns of `pred_cif`. Requires `censoring = "cox"`. With `NULL`, the
#'   Cox censoring model has no covariates. The Cox censoring model uses
#'   Breslow ties, so without covariates it equals exp(-Nelson-Aalen).
#' @param calibration_args Optional named list passed to [pa_calibration_cr()],
#'   e.g. \code{list(bandwidth = 0.2)}. These settings affect calibration only.
#'
#' @details
#' The generalized R-squared metrics use event-observation IPCW. Censoring is
#' estimated on the original follow-up times, treating any observed event type
#' as an observed event. The unnormalized weight is
#' \eqn{I(status_i > 0)/\widehat G(Y_i-\mid X_i)}; weights are normalized to sum
#' to one. A censored subject has zero weight even if follow-up exceeds `tau`.
#' Restricting the working outcome therefore does not limit censoring support
#' requirements to `tau`: support at contributing original event times is needed.
#' Cox weights require conditional independent censoring and a valid censoring
#' model; marginal weights require the corresponding marginal independence.
#' Censoring survival probabilities are floored at \eqn{10^{-4}}.
#' The interval outcome is the target-event time capped at `tau`, with competing
#' events assigned to `tau`. The point outcome is `I(T <= t_star, D = event_type)`,
#' including target events at the evaluation time, as in the CIF. Normalized ratio
#' estimates are not asserted to be exactly unbiased in finite samples.
#' Calibration uses the target-event CIF at or before \code{t_star} on the
#' supplied grid and evaluates the outcome at \code{t_star}. It uses the local
#' Aalen--Johansen estimator described in [pa_calibration_cr()]; \code{censoring}
#' and \code{censor_covariates} do not control this estimator. A failed curve
#' fit warns and returns missing calibration metrics without discarding others.
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item `"Metric"`: Name of the evaluated metric.
#'     \item `"Value"`: Computed value of the corresponding metric.
#'   }
#'
#' 
#' 
#' @export
#' @examples
#' # Competing risks: cause 1 is the event of interest, cause 2 the competing
#' # risk. Both cause-specific hazards are modelled.
#' cr <- pa_simulate_fine_gray(n = 150, v = 2, beta1 = c(0.5, -0.5),
#'                             beta2 = c(-0.3, 0.3), censor = 0.3, seed = 2025)
#' cr$time <- cr$obs.times
#' cr$status <- cr$obs.event
#' m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
#'                       data = cr, x = TRUE, y = TRUE)
#' m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
#'                       data = cr, x = TRUE, y = TRUE)
#' p <- pa_predict_cif(model1 = m1, model2 = m2, newdata = cr,
#'                     covs = c("X1", "X2"), event.type = 1, tau = max(cr$time))
#'
#' pa_survival_eval_cr(
#'   pred_cif = p$cif_pred[, -1], event_time = p$times,
#'   time.cif = p$cif_pred[, 1], status = p$status, event_type = 1
#' )
pa_survival_eval_cr <- function (pred_cif, event_time, time.cif, status,
                                           metrics = NULL,  t_star = NULL, 
                                           tau = NULL, event_type = 1,
                                           censoring = c("cox", "km"),
                                           censor_covariates = NULL,
                                           calibration_args = list())
  
{
  
  if (missing(event_time) || missing(pred_cif)) {
    stop("Please provide 'event_time',  and 'pred_cif' arguments.")
  }
  censoring <- match.arg(censoring)
  
  if(is.null(tau)) tau <- max(event_time)
  if(is.null(t_star)) t_star <- median(event_time)
  time_idx <- findInterval(t_star, time.cif)
  
  metrics_results <- list()
  
  valid_metrics <- c("pseudo_r2", "r_square", "l_square", 
                     "pseudo_r2_point", "r2_point", "l2_point",
                     "c_index", "brier_score", "td_auc", "ici", "e50", "e90")
  default_metrics <- c("pseudo_r2", "pseudo_r2_point", "c_index",
                       "brier_score", "td_auc", "ici")
  # Validation runs AFTER pa_normalize_metrics() only. A duplicate copy of this
  # block used to run before it, which rejected every legacy spelling
  # ("Harrells_C", "C_index", "Brier Score") on this path while they worked
  # everywhere else, and prevented the withdrawal message for r_e / r_sh from
  # ever being reached.
  metrics <- pa_normalize_metrics(metrics)
  


  if (is.null(metrics)){
    metrics <- default_metrics
  }
  else if ("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    invalid <- setdiff(metrics, valid_metrics)
    if (length(invalid) > 0) 
      stop("Invalid metrics: ", paste(invalid, collapse = ", "))
  }
  
  
  if("pseudo_r2" %in% metrics ||
     "r_square" %in% metrics || 
     "l_square" %in% metrics) {
    r_l_list <- pam.censor.cr(event_time, status, tau, 
                              pred_cif, time.cif, event_type,
                              X = censor_covariates, censoring = censoring)
  }
  
  if ("pseudo_r2" %in% metrics) {
    metrics_results$pseudo_r2 <- r_l_list$Pseudo.R
  } 
  if ("r_square" %in% metrics) {
    metrics_results$r_square <- r_l_list$R.squared
  }
  if ("l_square" %in% metrics) {
    metrics_results$l_square <- r_l_list$L.squared
  }
  
  if("pseudo_r2_point" %in% metrics ||
     "r2_point" %in% metrics || 
     "l2_point" %in% metrics) {
    r_p_list <- pam.censor.cr.point(ftime = event_time, 
                                    fstatus = status,
                                    tau = t_star, pred.cif = pred_cif[time_idx, ],
                                    event.type = event_type,
                                    X = censor_covariates, censoring = censoring)
  }
  
  if ("pseudo_r2_point" %in% metrics) {
    metrics_results$pseudo_r2_point <- r_p_list$Pseudo.R
  } 
  if ("r2_point" %in% metrics) {
    metrics_results$r2_point <- r_p_list$R.squared
  }
  if ("l2_point" %in% metrics) {
    metrics_results$l2_point <- r_p_list$L.squared
  }
  
  if ("c_index" %in% metrics) {
    pred_risk <- apply(pred_cif, 2, m_cif, time.cif=time.cif, tau=tau)
    C_index <- C_cr(
      time = event_time,
      status = status,
      predicted = pred_risk,
      tau = tau,
      Cause_int = event_type
    )
    metrics_results$"c_index" <- round(C_index, 4)
  }
  
  status.recode <- ifelse(status == event_type, -1, status)
  status.recode <- ifelse(status.recode > 0, 2, status.recode)
  status.recode <- ifelse(status.recode == -1, 1, status.recode)
  X <- if (time_idx == 0L) rep(0, length(event_time)) else pred_cif[time_idx, ]
  
  if ("brier_score" %in% metrics) {
    brier_result <- suppressMessages(tdROC::tdROC.cr(
      X = X,  
      Y = event_time,      
      delta = status.recode,
      tau = t_star, 
      nboot = 0))
    
    metrics_results$"brier_score" <- round(
      as.numeric(brier_result$calibration_res[1]), 4)
  }
  
  
  if ("td_auc" %in% metrics) {
    AUC_result <- suppressMessages(tdROC::tdROC.cr(
      X = X,  
      Y = event_time,      
      delta = status.recode,
      tau = t_star,       
      method = "both", 
      output = "AUC"   
    ))
    metrics_results$"td_auc" <- round(AUC_result$main_res$AUC.A.integral, 4)
    #metrics_results$"Time Dependent Auc Empirical" <- round(AUC_result$main_res$AUC.B.integral, 4)
  }
  
  calibration_metrics <- intersect(c("ici", "e50", "e90"), metrics)
  if (length(calibration_metrics)) {
    p <- pa_calibration_grid_risk(t(as.matrix(pred_cif)), time.cif, t_star,
                                  length(event_time))
    calibration <- pa_calibration_values(pa_calibration_cr, c(list(
      p = p, event_time = event_time, status = status, t_star = t_star,
      event_type = event_type), calibration_args))
    metrics_results[calibration_metrics] <- as.list(round(
      calibration[calibration_metrics], 4))
  }

  result_df <- data.frame(
    Metric = names(metrics_results),
    Value = as.numeric(unlist(metrics_results)),
    row.names = NULL,  
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}


##################################################
## calculate integral for int_0^tau t dF(t|x) 
##################################################
### input:
# x: input estimated CIF for subject i
# time.cif: time point of observation event 1
# tau: restricted time
#' @noRd
m_cif <- function(x, time.cif, tau){
  index <- which(time.cif <= tau)
  t(diff(c(0, time.cif[index], tau))) %*% c(1, (1-x)[index])
}




##################################################
### calculate PA measure for competing risk model
##################################################
### input:
# ftime: vector of failure/censoring times
# fstatus: vector with a unique code for each failure type (1,2,3...)
#          and censored observations (0)
# tau: restricted time
# time.cif: time point of observation event of interest
# pred.cif: estimated CIF for all subject
#' @noRd
pam.censor.cr <- function(ftime, fstatus, tau, 
                          pred.cif, time.cif,
                          event.type = 1, X = NULL,
                          censoring = c("cox", "km")){
  

  ftime.new <- ifelse(ftime > tau, tau, ftime)
  ftime.new <- ifelse(ftime.new <= tau & !(fstatus %in% c(event.type, 0)), 
                      tau, ftime.new)
  

  pred <- apply(pred.cif, 2, m_cif, time.cif = time.cif, tau = tau)
  

  ipcw <- pa_cr_event_ipcw(ftime, fstatus, X = X, censoring = censoring)
  out <- pa_cr_weighted_r2(ftime.new, pred, ipcw$weight)

  out <- c(out,
           Pseudo.R = round(as.numeric(out$R.squared) *
                              as.numeric(out$L.squared), digits = 4))
  return(out)
}


####################################################################
### calculate competing risk PA measure for a given time point
####################################################################
#' @noRd
pam.censor.cr.point <- function(ftime, fstatus, tau, 
                                pred.cif,
                                event.type = 1, X = NULL,
                                censoring = c("cox", "km")){
  
  # obtain PA measure
  i.obs <- ifelse(ftime <= tau & fstatus == event.type, 1, 0)
  ipcw <- pa_cr_event_ipcw(ftime, fstatus, X = X, censoring = censoring)
  out <- pa_cr_weighted_r2(i.obs, pred.cif, ipcw$weight)
  out <- c(out,
           Pseudo.R = round(as.numeric(out$R.squared) *
                              as.numeric(out$L.squared), digits = 4))
  return(out)
}


####################################################################
### calculate PA measure for a given time point
####################################################################
#' @noRd
pam.censor.point <- function(i.obs, i.predict, y, delta){
  
  y.sorted<-sort(y)
  delta.sorted<-delta[order(y)]
  
  #KM estimate for censoring distribution
  fit.km.censoring <- survfit(Surv(y.sorted, 1-delta.sorted) ~ 1  )
  km.censoring.minus <- pa_censoring_left_limit(fit.km.censoring, y.sorted)
  ratio.km<- delta.sorted/km.censoring.minus
  
  ratio.km[is.nan(ratio.km)]<-0
  
  weight.km<-ratio.km/(sum(ratio.km))
  
  ### now change y (observe event time) to indicator version
  
  y.sorted<-i.obs[order(y)]
  y.predict.sorted<-i.predict[order(y)]
  
  wls.fitted<-tryCatch(lm(y.sorted~y.predict.sorted,weights=weight.km),error=function(e){return(c(NA,NA))})
  calibrate.fitted<-tryCatch(predict(wls.fitted),error=function(e){return(c(NA,NA))})
  
  num.rho2<-sum(weight.km*(calibrate.fitted-sum(weight.km*y.sorted))^2)
  denom.rho2<-sum(weight.km*(y.sorted-sum(weight.km*y.sorted))^2)
  
  R2 <- round(num.rho2 / denom.rho2, digits = 4)
  
  
  num.L2<- sum(weight.km*(y.sorted-calibrate.fitted)^2)
  denom.L2<- sum(weight.km*(y.sorted-y.predict.sorted)^2)
  L2 <- round(num.L2 / denom.L2, digits = 4)
  
  return(list(R.squared=R2,L.squared=L2))
}


#' @noRd
pam.censor<-function(y,y.predict,delta){
  
  
  
  y.sorted<-sort(y)
  delta.sorted<-delta[order(y)]
  y.predict.sorted<-y.predict[order(y)]
  
  fit.km.censoring <- survfit(Surv(y.sorted, 1-delta.sorted) ~ 1  )
  km.censoring.minus <- pa_censoring_left_limit(fit.km.censoring, y.sorted)
  ratio.km<- delta.sorted/km.censoring.minus
  
  ratio.km[is.nan(ratio.km)]<-0
  
  weight.km<-ratio.km/(sum(ratio.km))
  
  
  wls.fitted<-tryCatch(lm(y.sorted~y.predict.sorted,weights=weight.km),error=function(e){return(c(NA,NA))})
  calibrate.fitted<-tryCatch(predict(wls.fitted),error=function(e){return(c(NA,NA))})
  
  
  
  num.rho2<-sum(weight.km*(calibrate.fitted-sum(weight.km*y.sorted))^2)
  denom.rho2<-sum(weight.km*(y.sorted-sum(weight.km*y.sorted))^2)
  
  R2 <- round(num.rho2 / denom.rho2, digits = 4)
  
  num.L2<- sum(weight.km*(y.sorted-calibrate.fitted)^2)
  denom.L2<- sum(weight.km*(y.sorted-y.predict.sorted)^2)
  L2 <- round(num.L2 / denom.L2, digits = 4)
  
  return(list(R.squared=R2,L.squared=L2))
}

#' @noRd
C_cr <- function(time, status, predicted, tau = NULL, Cause_int = 1) {
  if (any(is.na(time))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(status))) {
    stop("The input vector cannot have NA")
  }
  if (any(!(status %in% c(0, 1, 2)))) {
    stop("The status must be 0 or 1 or 2")
  }
  if (any(is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input of Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }
  
  KM <- summary(survfit(Surv(time, status==0)~1))
  df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
  
  Time_survival <- time
  Censoring <- ifelse(status == 0, 0, 1)
  Cause <- ifelse(status == 2, 2, 1)
  Prediction <- -log(predicted)
  if(is.null(tau)) tau <- max(Time_survival) + 1
  n <- length(Prediction)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)
  Num <- 0
  Den <- 0
  
  w1 <- rep(0, n)
  w2 <- rep(0, n)
  for(i in 1:n){
    w1[i] <- min(df_KM$surv[df_KM$time < Time_survival[i]]) 
    w2[i] <- min(df_KM$surv[df_KM$time <= Time_survival[i]])
  }
  w3 <- w1 %*% t(w1)
  w1 <- w1 * w2
  
  for (i in 1:n) {
    # Only invert probabilities for pairs that enter the estimator. Censoring
    # survival may be zero after tau, where rows do not contribute.
    if (Time_survival[i] <= tau && Cause[i] == Cause_int &&
        Censoring[i] == 1) {
      later <- which(Time_survival[i] < Time_survival)
      competing <- which(Time_survival[i] >= Time_survival &
                           Cause != Cause_int & Censoring == 1)
      if ((length(later) && w1[i] <= 0) ||
          any(w3[i, competing] <= 0)) {
        stop("Censoring survival must be positive for comparable pairs.")
      }
      # Wolbers et al. (2014), equation (3.4): reciprocal products.
      A[i, later] <- 1 / w1[i]
      B[i, competing] <- 1 / w3[i, competing]
    }
    Q[i, which(Prediction[i] > Prediction)] <- 1
  }
  for (i in 1:n) {
    if (Time_survival[i] <= tau && Cause[i] == Cause_int && 
        Censoring[i] == 1) {
      N_t[i, ] <- 1
    }
  }
  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t
  Num <- sum(Num_mat)
  Den <- sum(Den_mat)
  return(Num/Den)
}

#' Summarize multiple competing-risks models into a wide comparison table
#'
#' @description
#' Calls \code{pa_survival_eval_cr()} for each model in a named list and
#' returns a wide table where rows are metrics and columns are model names.
#'
#' @param models A **named list**; each element is a list with components:
#'   \itemize{
#'     \item \code{pred_cif}   -- numeric matrix of predicted CIF values (rows: time grid; cols: subjects)
#'     \item \code{time.cif}   -- numeric vector of the time grid corresponding to \code{pred_cif} rows
#'     \item \code{event_time} -- numeric vector of observed times
#'     \item \code{status}     -- integer vector of event codes (0=censoring; \code{event_type}=target; others=competing)
#'   }
#' @param metrics Character vector of metrics to compute (passed through to
#'   \code{pa_survival_eval_cr()}); use \code{"all"} for all supported.
#' @param t_star (Optional) Evaluation time for Brier score, AUC, point R2,
#'   and calibration. Defaults to the median observed follow-up.
#' @param tau (Optional) The maximum prediction time horizon for R2 and C index. Default is maximum of observation time.
#' @param event_type Integer code for the cause of interest (passed through).
#' @param digits Number of decimal places to round the metric values (default = 2).
#' @param censoring Censoring estimator for generalized R-squared metrics,
#'   `"cox"` (default) or `"km"`; passed to [pa_survival_eval_cr()].
#' @param censor_covariates Optional numeric censoring-model covariates shared
#'   across models. Rows must match each model's subject order. Passed to
#'   [pa_survival_eval_cr()]; used only with `censoring = "cox"`.
#' @param calibration_args Optional named list passed to [pa_calibration_cr()].
#' @details See [pa_survival_eval_cr()] for the event-observation weighting
#'   rule, original-time support requirements and probability floor.
#'
#' @return A data.frame with:
#' \itemize{
#'   \item each row = metric
#'   \item each column = model name
#'   \item cell values = metric values
#' }
#' @export
#' @examples
#' # Competing risks: cause 1 is the event of interest, cause 2 the competing
#' # risk. Both cause-specific hazards are modelled.
#' cr <- pa_simulate_fine_gray(n = 150, v = 2, beta1 = c(0.5, -0.5),
#'                             beta2 = c(-0.3, 0.3), censor = 0.3, seed = 2025)
#' cr$time <- cr$obs.times
#' cr$status <- cr$obs.event
#' m1 <- survival::coxph(survival::Surv(time, status == 1) ~ X1 + X2,
#'                       data = cr, x = TRUE, y = TRUE)
#' m2 <- survival::coxph(survival::Surv(time, status == 2) ~ X1 + X2,
#'                       data = cr, x = TRUE, y = TRUE)
#' p <- pa_predict_cif(model1 = m1, model2 = m2, newdata = cr,
#'                     covs = c("X1", "X2"), event.type = 1, tau = max(cr$time))
#'
#' pa_summarize_cr(list(csh = p), event_type = 1)
pa_summarize_cr <- function(models,
                           metrics = NULL,
                           t_star = NULL,
                           tau = NULL,
                           event_type = 1,
                           digits = 2,
                           censoring = c("cox", "km"),
                           censor_covariates = NULL,
                           calibration_args = list()) {
  censoring <- match.arg(censoring)
  if (!is.list(models) || length(models) == 0)
    stop("'models' must be a non-empty named list.")
  if (is.null(names(models)) || any(names(models) == ""))
    names(models) <- paste0("Model_", seq_along(models))
  
  required <- c("cif_pred", "times", "status")
  
  eval_one <- function(mod, name) {
    if (!all(required %in% names(mod))) {
      stop(sprintf("Model '%s' must contain: %s",
                   name, paste(required, collapse = ", ")))
    }
    
    res <- pa_survival_eval_cr(
      pred_cif   = mod$cif_pred[, -1],
      event_time = mod$times,
      time.cif   = mod$cif_pred[, 1],
      status     = mod$status,
      metrics    = metrics,
      t_star     = t_star,
      tau        = tau,
      event_type = event_type,
      censoring = censoring,
      censor_covariates = censor_covariates,
      calibration_args = calibration_args
    )
    
    if (!all(c("Metric", "Value") %in% names(res)))
      stop(sprintf("Unexpected result from pa_survival_eval_cr() for '%s'", name))
    
    res$Model <- name
    res[, c("Model", "Metric", "Value")]
  }
  
  # evaluate each model
  res_list <- mapply(eval_one, models, names(models), SIMPLIFY = FALSE)
  res_long <- do.call(rbind, res_list)
  rownames(res_long) <- NULL
  
  # pivot to wide: rows = Metric, cols = Model, cells = Value
  res_wide <- reshape(
    res_long,
    idvar = "Metric",
    timevar = "Model",
    direction = "wide"
  )
  
  # clean column names "Value.Model" -> "Model"
  names(res_wide) <- sub("^Value\\.", "", names(res_wide))
  rownames(res_wide) <- NULL
  res_wide <- res_wide[, c("Metric", setdiff(names(res_wide), "Metric"))]
  
  # optional: enforce a pleasant metric order if present
  preferred_order <- c("pseudo_r2", "r_square", "l_square",
                       "pseudo_r2_point", "r2_point", "l2_point",
                       "c_index", "brier_score", "td_auc", "ici", "e50", "e90")
  present <- intersect(preferred_order, res_wide$Metric)
  others  <- setdiff(res_wide$Metric, preferred_order)
  res_wide$Metric <- factor(res_wide$Metric, levels = c(present, sort(others)))
  res_wide <- res_wide[order(res_wide$Metric), ]
  res_wide$Metric <- as.character(res_wide$Metric)
  
  # round numeric columns
  numeric_cols <- setdiff(names(res_wide), "Metric")
  res_wide[numeric_cols] <- lapply(res_wide[numeric_cols], function(x) round(as.numeric(x), digits))
  
  res_wide
}
