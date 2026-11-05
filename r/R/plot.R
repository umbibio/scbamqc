BLUE <- "#2a78d6"
ORANGE <- "#eb6834"
GRAY <- "#9a9996"
INK <- "#52514e"

need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install the 'ggplot2' package to make plots.", call. = FALSE)
  }
}

theme_scbamqc <- function() {
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      legend.position = "top",
      plot.title = ggplot2::element_text(size = 10)
    )
}

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

#' Quality score against the alignment and count metrics
#'
#' Plots `qc_score` against each metric in `x`, one panel for each metric.
#' Selected barcodes are blue. Other scored barcodes are gray.
#'
#' @inheritParams gmm_labels
#' @param x Metadata columns to plot on the x-axis.
#' @param selected_key Logical metadata column with the selected barcodes, or
#'   `NULL`.
#' @param min_score,min_counts If set, dotted lines show these limits.
#' @param log_x Columns to plot as `log10`.
#' @param size,alpha Point size and opacity.
#' @return A ggplot object.
#' @export
qc_score_scatter <- function(object,
                             x = c("total_counts", "mean_clip_frac",
                                   "std_clip_frac", "mean_base_qual", "mean_nm"),
                             score_key = "qc_score", selected_key = "selected",
                             min_score = NULL, min_counts = NULL,
                             log_x = "total_counts", size = 0.3, alpha = 0.1) {
  need_ggplot2()
  check_input(object)
  meta <- get_meta(object)
  need_cols(meta, c(score_key, x))

  y <- meta[[score_key]]
  scored <- !is.na(y)
  has_sel <- !is.null(selected_key) && selected_key %in% colnames(meta)
  sel <- if (has_sel) meta[[selected_key]] %in% TRUE else rep(TRUE, nrow(meta))

  label <- ifelse(x %in% log_x, paste0("log10 ", x), x)
  long <- do.call(rbind, lapply(seq_along(x), function(i) {
    v <- as.numeric(meta[[x[i]]])
    if (x[i] %in% log_x) v <- log10(pmax(v, 1))
    data.frame(metric = label[i], value = v, score = y, selected = sel)[scored, ]
  }))
  long$metric <- factor(long$metric, levels = label)
  # Draw the unselected barcodes first, so that the selected stay visible.
  long <- long[order(long$selected), ]

  title <- paste0("N scored = ", fmt(sum(scored)))
  if (has_sel) title <- paste0(title, ", selected (blue) = ", fmt(sum(sel & scored)))

  p <- ggplot2::ggplot(long, ggplot2::aes(value, score,
                                          colour = selected)) +
    ggplot2::geom_point(size = size, alpha = alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = GRAY, `TRUE` = BLUE),
                                 guide = "none") +
    ggplot2::facet_wrap(~metric, nrow = 1, scales = "free_x",
                        strip.position = "bottom") +
    ggplot2::labs(x = NULL, y = score_key, title = title) +
    theme_scbamqc()
  if (!is.null(min_score)) {
    p <- p + ggplot2::geom_hline(yintercept = min_score, colour = INK,
                                 linetype = "dotted")
  }
  if (!is.null(min_counts) && "total_counts" %in% x) {
    xv <- if ("total_counts" %in% log_x) log10(max(min_counts, 1)) else min_counts
    vline <- data.frame(metric = factor(label[x == "total_counts"], levels = label),
                        value = xv)
    p <- p + ggplot2::geom_vline(data = vline,
                                 ggplot2::aes(xintercept = value),
                                 colour = INK, linetype = "dotted")
  }
  p
}

#' Histogram of the quality score by selection outcome
#'
#' A stacked histogram of the scored barcodes. Gray: not selected. Orange:
#' selected, with `qc_score <= threshold`. Blue: retained
#' (`qc_score > threshold`). A dashed line shows `threshold`.
#'
#' @inheritParams qc_score_scatter
#' @param threshold Score limit for the retained cells.
#' @param bins Number of bins from 0 to 1.
#' @param title Plot title. The default shows the number of retained cells.
#' @return A ggplot object.
#' @export
qc_score_hist <- function(object, threshold = 0.9, score_key = "qc_score",
                          selected_key = "selected", bins = 50, title = NULL) {
  need_ggplot2()
  check_input(object)
  meta <- get_meta(object)
  need_cols(meta, score_key)

  y <- meta[[score_key]]
  scored <- !is.na(y)
  sel <- if (selected_key %in% colnames(meta)) meta[[selected_key]] %in% TRUE else scored
  retained <- scored & sel & y > threshold
  candidate <- scored & sel & !retained
  background <- scored & !sel

  levels <- c(
    sprintf("not selected (%s)", fmt(sum(background))),
    sprintf("candidate, \u2264 %s (%s)", format(threshold), fmt(sum(candidate))),
    sprintf("retained, > %s (%s)", format(threshold), fmt(sum(retained)))
  )
  group <- ifelse(retained, levels[3], ifelse(candidate, levels[2], levels[1]))
  df <- data.frame(score = y[scored], group = factor(group[scored], levels = levels))
  if (is.null(title)) {
    title <- sprintf("%s retained of %s scored", fmt(sum(retained)), fmt(sum(scored)))
  }

  ggplot2::ggplot(df, ggplot2::aes(score, fill = group)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, length.out = bins + 1),
                            colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = stats::setNames(c(GRAY, ORANGE, BLUE), levels),
                               name = NULL, drop = FALSE) +
    ggplot2::geom_vline(xintercept = threshold, colour = INK, linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = score_key, y = "Barcodes", title = title) +
    theme_scbamqc() +
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(0.02, 0.98),
                   legend.justification = c(0, 1))
}

utils::globalVariables(c("value", "score", "selected", "group"))
