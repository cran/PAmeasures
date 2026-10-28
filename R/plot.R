#' Quick Diagnostic Plot for Survival Model Predicted Probabilities
#'
#' @description
#' Produces a quick diagnostic plot visualizing predicted survival or event
#' probabilities (e.g., from a Weibull AFT model) against the subject-specific
#' linear predictors (risk scores). The plot combines a smoothed line of
#' predicted values and points representing observed times, using shape to
#' indicate censoring status. This visualization provides a simple check of
#' how predicted risks correspond to observed outcomes.
#'
#' @param data A data frame (or list-like object) containing at least the
#'   following columns:
#'   \describe{
#'     \item{`linear.pred`}{Model linear predictor (risk score).}
#'     \item{`pred`}{Predicted probability or median survival value.}
#'     \item{`times`}{Observed event or censoring time.}
#'     \item{`status`}{Event indicator (0 = censored, 1 = event; optionally
#'       can include more levels if defined in \code{label_name}).}
#'   }
#' @param title Optional plot title (default = \code{NULL}).
#' @param xlab,ylab Character strings specifying the axis labels
#'   (defaults are \code{"Risk Score"} for x and \code{"Days"} for y).
#' @param label_name Character vector of labels corresponding to event
#'   status codes (e.g., \code{c("censored", "event")}); used in the legend.
#' @param levels Numeric or character vector giving the ordering of status
#'   levels; defaults inferred from the length of \code{label_name}.
#' @param shape_style Numeric vector of plotting symbols (pch values)
#'   corresponding to the levels of \code{status}. If \code{NULL}, defaults
#'   to \code{c(1, 19)} for two levels (censored vs event) or
#'   \code{c(2, 19, 1)} for three levels.
#' @param legend_name Character string for the legend title (default = "status").
#' @param invert_linear Logical; if \code{TRUE} (default), plots the negative
#'   of \code{linear.pred} so that higher values correspond to higher risk.
#'   Set to \code{FALSE} if you want higher values correspond to lower risk.
#' @param sample_index Optional integer vector selecting which subjects to
#'   plot, e.g. \code{sample(seq_len(n), 200)} to thin a crowded chart. If
#'   \code{NULL} (default), every subject is plotted.
#' @param restrict_time Optional numeric. If supplied, observed times above
#'   this value are capped at it and a dashed horizontal reference line is
#'   drawn, which keeps a few long follow-up times from compressing the plot.
#' @details
#' The function is primarily designed for quick visual diagnostics in survival
#' model development. The line (\code{pred} vs \code{linear.pred}) shows the
#' fitted or predicted trend, while the overlaid points (\code{times} vs
#' \code{linear.pred}) illustrate observed event times, distinguished by
#' censoring status. For more elaborate visualization (e.g., confidence bands,
#' multiple models), use faceting or combine plots with \pkg{patchwork}.
#'
#' @return A \code{ggplot} object that can be further customized.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_shape_manual
#'   theme_classic xlab ylab ggtitle theme element_text
#' @export
#' @examples
#' d <- pa_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' d <- d[, c("time", "status", "x1", "x2")]
#' fit <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
#'                        data = d, x = TRUE, y = TRUE)
#' pred <- pa_predict_coxph(model = fit, covs = c("x1", "x2"), new_data = d)
#'
#' # observed time against predicted risk score, shaped by event status
#' pa_plot_pred(pred, title = "Cox model")
pa_plot_pred <- function(data,
                      title = NULL,
                      xlab = "Risk Score",
                      ylab = "Days",
                      label_name = c("censored", "event"),
                      levels = NULL,
                      shape_style = NULL,
                      legend_name = "status",
                      invert_linear = TRUE, 
                      sample_index = NULL,
                      restrict_time = NULL) {
  if(is.null(shape_style) & length(label_name) == 2) shape_style = c(1, 19)
  if(is.null(shape_style) & length(label_name) == 3) shape_style = c(2, 19, 1)
  if(is.null(levels) & length(label_name) == 2) levels = c(0, 1)
  if(is.null(levels) & length(label_name) == 3) levels = c(0, 1, 2)
  # check required columns
  required_cols <- c("linear.pred", "pred", "times", "status")
  if (!all(required_cols %in% names(data))) {
    stop("Input data must contain columns: ",
         paste(required_cols, collapse = ", "))
  }
  
  # optionally invert the linear predictor
  x_var <- if (invert_linear) -data$linear.pred else data$linear.pred

  # NULL means "plot every subject". Without this, x_var[NULL] returns a
  # zero-length vector and the plot renders empty.
  if (is.null(sample_index)) sample_index <- seq_along(x_var)

  df <- data.frame(
    x_var = x_var[sample_index],
    pred  = data$pred[sample_index],
    times = data$times[sample_index],
    status = data$status[sample_index]
  )
  if(!is.null(restrict_time)){
    df$times <- ifelse(df$times <= restrict_time, df$times, restrict_time)
  }
  ggplot2::ggplot(
    data = df
  ) +
    ggplot2::geom_line(ggplot2::aes(x = x_var, y = pred)) +
    ggplot2::geom_point(ggplot2::aes(
      x = x_var, y = times,
      shape = factor(status, levels = levels,
                     labels = label_name)
    )) +
    ggplot2::scale_shape_manual(legend_name, values = shape_style) +
    ggplot2::theme_classic() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) -> p
  if(!is.null(restrict_time)){
    p + ggplot2::geom_hline(aes(yintercept = restrict_time),
                            linetype = "dashed", color = "black") -> p
  } 
  return(p)
}


#' Arrange multiple prediction plots with letter tags
#'
#' @description
#' Builds one panel per data frame using \code{plot_fun} (default: \code{pa_plot_pred})
#' and arranges them with letter tags. The function collects a single shared legend
#' and applies global styling (titles, axis labels, legend options) across panels.
#'
#' @param data_list List of data frames, each containing columns:
#'   \code{linear.pred}, \code{pred}, \code{times}, and \code{status}.
#' @param titles Optional character vector, same length as \code{data_list};
#'   if \code{NULL}, defaults to "Panel 1", "Panel 2", etc.
#' @param plot_fun Function that creates a single panel (default = \code{pa_plot_pred}).
#' @param ncol Number of columns in the panel layout (default = 2).
#' @param tag_levels Letter style for tags: \code{"a"}, \code{"A"}, \code{"1"}, \code{"i"}, or \code{"I"}.
#'
#' @param invert_linear Logical or logical vector; if length 1, recycled to all panels.
#' @param xlab,ylab Axis labels applied to all panels (defaults are "Risk Score" and "Days").
#' @param label_name Character vector of labels for status legend (e.g., \code{c("censored", "event")}).
#' @param levels Numeric vector for \code{status} level ordering.
#' @param shape_style Numeric vector of plotting symbols for censoring/event status.
#' @param legend_name Character string for the legend title (default = "status").
#' @param sample_size Optional integer. If supplied, that many subjects are
#'   sampled at random for each panel using \code{seed}; otherwise every
#'   subject is plotted.
#' @param restrict_time Optional numeric. If supplied, observed times above
#'   this value are capped at it and a dashed horizontal reference line is
#'   drawn, which keeps a few long follow-up times from compressing the plot.
#' @param seed Optional integer. If supplied, the random number generator is
#'   seeded with it before \code{sample_size} subjects are drawn, so the
#'   sampled points are reproducible. The default \code{NULL} uses the current random number stream
#'   without resetting it; random draws advance its state.
#'
#' @return A patchwork \code{ggplot} object combining all panels.
#' @importFrom patchwork plot_annotation wrap_plots plot_layout
#' @export
#'
#' @examples
#' # One panel per model. The previous example referenced undefined objects and
#' # so could not be executed; this one runs.
#' d <- pa_sim_cox_weibull(n = 150, pi_c = 0.3, v = 2,
#'                         beta = c(0.5, -0.5), seed = 2025)
#' d <- d[, c("time", "status", "x1", "x2")]
#' cox <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
#'                        data = d, x = TRUE, y = TRUE)
#' wei <- survival::survreg(survival::Surv(time, status) ~ x1 + x2,
#'                          data = d, dist = "weibull", x = TRUE, y = TRUE)
#' p_cox <- pa_predict_coxph(model = cox, covs = c("x1", "x2"), new_data = d)
#' p_wei <- pa_predict_survreg(model = wei, covs = c("x1", "x2"), new_data = d)
#'
#' pa_plot_summary(
#'   list(p_cox, p_wei),
#'   titles = c("Cox PH", "Weibull AFT"),
#'   ncol = 2,
#'   legend_name = "Status"
#' )
pa_plot_summary <- function(data_list,
                              titles = NULL,
                              plot_fun = pa_plot_pred,
                              ncol = 2,
                              tag_levels = "a",
                              invert_linear = TRUE,
                              xlab = "Risk Score",
                              ylab = "Days",
                              label_name = c("Censored", "Event"),
                              levels = NULL,
                              shape_style = NULL,
                              sample_size = NULL,
                              restrict_time = NULL,
                              legend_name = "Status", seed = NULL) {
  if(is.null(shape_style) & length(label_name) == 2) shape_style = c(1, 19)
  if(is.null(shape_style) & length(label_name) == 3) shape_style = c(2, 19, 1)
  if(is.null(levels)) levels = c(0:(length(label_name)-1))
  if(is.null(ncol)) ncol = length(label_name)
  # Basic checks
  if (!is.list(data_list) || length(data_list) == 0)
    stop("data_list must be a non-empty list of data frames.")
  
  k <- length(data_list)
  
  if (is.null(titles)) titles <- paste("Panel", seq_len(k))
  if (length(titles) != k)
    stop("titles must have the same length as data_list or be NULL.")
  
  # Recycle invert_linear if scalar
  if (length(invert_linear) == 1) invert_linear <- rep(invert_linear, k)
  if (length(invert_linear) != k)
    stop("invert_linear must be length 1 or length(data_list).")
  
  if(!is.null(sample_size)){
    if (!is.null(seed)) set.seed(seed)
    sample_index <- sample(1:length(data_list[[1]][["times"]]), size = sample_size)
  }else{
    sample_index <- 1:length(data_list[[1]][["times"]])
  } 
  # Create each plot
  plots <- Map(function(d, ttl, inv) {
    plot_fun(
      d,
      title = ttl,
      invert_linear = inv,
      xlab = xlab,
      ylab = ylab,
      label_name = label_name,
      levels = levels,
      shape_style = shape_style,
      legend_name = legend_name,
      sample_index = sample_index,
      restrict_time = restrict_time
    )
  }, data_list, titles, invert_linear)
  
  # Arrange and annotate with patchwork
  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      tag_levels = tag_levels,
      tag_prefix = "(",
      tag_suffix = ")"
    ) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

