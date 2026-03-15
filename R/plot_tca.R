#' @title TCA Plotting Functions
#' @description Visualisation of Transmission Channel Analysis results
#'   using \pkg{ggplot2}.
#' @name tca_plotting
NULL

#' Plot TCA Channel Decomposition
#'
#' Creates a stacked bar chart (for exhaustive modes) or a line chart
#' (for overlapping mode) showing channel contributions over horizons.
#'
#' @param x A \code{tca_result} object from \code{\link{tca_analyze}}.
#' @param target Integer index of the response variable to plot
#'   (default: first variable after shock).
#' @param type Plot type: \code{"bar"} for stacked bar chart (default
#'   for exhaustive modes) or \code{"line"} for overlaid lines (default
#'   for overlapping mode). If \code{NULL}, chosen automatically.
#' @param title Custom plot title (optional).
#' @param colors Named character vector of colours for channels (optional).
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' K <- 4
#' A1 <- matrix(c(0.7,-0.1,0.05,-0.05, -0.3,0.6,0.10,-0.10,
#'                  -0.2,0.1,0.70,0.05, -0.1,0.2,0.05,0.65), K, K, byrow=TRUE)
#' Sigma <- matrix(c(1,0.3,0.2,0.1, 0.3,1.5,0.25,0.15,
#'                    0.2,0.25,0.8,0.1, 0.1,0.15,0.1,0.6), K, K, byrow=TRUE)
#' Phi0 <- t(chol(Sigma))
#' sf <- tca_systems_form(Phi0, list(A1), h = 20)
#' res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
#'                     intermediates = c(2, 4), K = K, h = 20,
#'                     order = 1:K, mode = "exhaustive_4way",
#'                     var_names = c("IntRate","GDP","Inflation","Wages"))
#' plot_tca(res, target = 3)
plot_tca <- function(x, target = NULL, type = NULL, title = NULL,
                      colors = NULL) {
  if (!inherits(x, "tca_result")) {
    stop("x must be a 'tca_result' object from tca_analyze().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  if (is.null(target)) {
    # Default: first intermediate or first non-shock variable
    target <- if (x$from == 1) 2L else 1L
  }
  if (is.character(target)) {
    target <- match(target, x$var_names)
    if (is.na(target)) stop("target variable not found.")
  }

 # Auto-select plot type
  if (is.null(type)) {
    type <- if (x$mode == "overlapping") "line" else "bar"
  }

  horizons <- 0:x$h
  n_ch <- length(x$channel_names)

  # Build data frame
  df_list <- vector("list", n_ch)
  for (i in seq_len(n_ch)) {
    nm <- x$channel_names[i]
    df_list[[i]] <- data.frame(
      horizon = horizons,
      value   = x$irf_channels[[nm]][, target],
      channel = nm,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, df_list)
  df$channel <- factor(df$channel, levels = x$channel_names)

  # Default title
  if (is.null(title)) {
    title <- sprintf("TCA: %s shock -> %s (%s)",
                     x$var_names[x$from], x$var_names[target], x$mode)
  }

  # Default colors
  if (is.null(colors)) {
    default_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                     "#FF7F00", "#A65628", "#F781BF", "#999999")
    colors <- stats::setNames(default_pal[seq_len(n_ch)], x$channel_names)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$horizon,
                                         y = .data$value,
                                         fill = .data$channel,
                                         colour = .data$channel))

  if (type == "bar") {
    p <- p +
      ggplot2::geom_bar(stat = "identity", position = "stack",
                        width = 0.8) +
      ggplot2::scale_fill_manual(values = colors)
  } else {
    # Total IRF as background
    df_total <- data.frame(horizon = horizons,
                           value = x$irf_total[, target])
    p <- p +
      ggplot2::geom_line(data = df_total,
                         ggplot2::aes(x = .data$horizon, y = .data$value),
                         inherit.aes = FALSE,
                         linewidth = 1.2, colour = "grey40",
                         linetype = "dashed") +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::scale_colour_manual(values = colors)
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "grey50", linewidth = 0.3) +
    ggplot2::labs(title = title, x = "Horizon", y = "IRF",
                  fill = "Channel", colour = "Channel") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}
