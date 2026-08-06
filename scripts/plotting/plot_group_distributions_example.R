#!/usr/bin/env Rscript

# Plot cohort-stratified group distributions with lightweight Wilcoxon annotations.
# Analysis context: cohort-stratified group comparison.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[[1]]
source(file.path(dirname(normalizePath(script_file)), "plotting_common.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing R package: ggplot2", call. = FALSE)
if (!requireNamespace("svglite", quietly = TRUE)) stop("Missing R package: svglite", call. = FALSE)

args <- parse_cli(
  commandArgs(trailingOnly = TRUE),
  required = c("input", "output_prefix"),
  defaults = list(title = "Group distributions across cohorts", formats = "pdf,svg,png", seed = "42")
)
set.seed(as.integer(args$seed))

data <- read_table_file(args$input)
require_columns(data, c("cohort", "metric", "group", "value"), "group distributions")
data <- require_nonempty_text(data, c("cohort", "metric", "group"), "group distributions")
data$value <- suppressWarnings(as.numeric(data$value))
if (any(!is.finite(data$value))) stop("group distributions value must be finite", call. = FALSE)

group_order <- unique(as.character(data$group))
if (length(group_order) != 2L) {
  stop("This example requires exactly two group values", call. = FALSE)
}
data$group <- factor(data$group, levels = group_order)
data$cohort <- factor(data$cohort, levels = unique(as.character(data$cohort)))
data$metric <- factor(data$metric, levels = unique(as.character(data$metric)))

split_data <- split(data, interaction(data$cohort, data$metric, drop = TRUE))
statistics <- lapply(split_data, function(frame) {
  observed_groups <- unique(as.character(frame$group))
  if (length(observed_groups) != 2L || any(table(frame$group) < 2L)) {
    stop("Every cohort/metric panel must contain both groups with at least two observations", call. = FALSE)
  }
  test <- stats::wilcox.test(value ~ group, data = frame, exact = FALSE)
  value_range <- range(frame$value)
  margin <- max(0.05, diff(value_range) * 0.12)
  data.frame(
    cohort = as.character(frame$cohort[[1]]),
    metric = as.character(frame$metric[[1]]),
    x = 1.5,
    y = max(frame$value) + margin,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})
statistics <- do.call(rbind, statistics)
statistics$q_value <- stats::p.adjust(statistics$p_value, method = "BH")
statistics$label <- vapply(statistics$q_value, format_q_value, character(1))

palette <- stats::setNames(c("#35637A", "#A23D52"), group_order)
plot <- ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group, color = group)) +
  ggplot2::geom_violin(width = 0.85, alpha = 0.18, linewidth = 0.5, trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.28, alpha = 0.72, outlier.shape = NA, linewidth = 0.45) +
  ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.2, alpha = 0.62, show.legend = FALSE) +
  ggplot2::geom_text(
    data = statistics,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    color = "#303030"
  ) +
  ggplot2::facet_grid(metric ~ cohort, scales = "free_y") +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::scale_color_manual(values = palette) +
  ggplot2::labs(title = args$title, x = NULL, y = "Value", fill = "Group", color = "Group") +
  ggplot2::theme_bw(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F2F2F2", color = "#D0D0D0"),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

save_ggplot_bundle(plot, args$output_prefix, args$formats, width = 10.0, height = 6.2)
