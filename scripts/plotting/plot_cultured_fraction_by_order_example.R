#!/usr/bin/env Rscript

# Plot precomputed culture-representation proportions and confidence intervals by bacterial order.
# Analysis context: phylogenomic culture-representation summary.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[[1]]
source(file.path(dirname(normalizePath(script_file)), "plotting_common.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing R package: ggplot2", call. = FALSE)
if (!requireNamespace("scales", quietly = TRUE)) stop("Missing R package: scales", call. = FALSE)
if (!requireNamespace("svglite", quietly = TRUE)) stop("Missing R package: svglite", call. = FALSE)

args <- parse_cli(
  commandArgs(trailingOnly = TRUE),
  required = c("input", "output_prefix"),
  defaults = list(
    title = "Culture representation by bacterial order",
    subtitle = "Points show the cultured proportion with Wilson 95% confidence intervals",
    formats = "pdf,svg,png"
  )
)

data <- read_table_file(args$input)
required_columns <- c(
  "dataset", "order", "cultured", "total", "ci_low", "ci_high",
  "q_value", "omnibus_p", "omnibus_q"
)
require_columns(data, required_columns, "culture fraction input")
data <- require_nonempty_text(data, c("dataset", "order"), "culture fraction input")

numeric_columns <- setdiff(required_columns, c("dataset", "order"))
for (column in numeric_columns) {
  data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
}
if (any(!is.finite(as.matrix(data[numeric_columns])))) {
  stop("culture fraction numeric columns must be finite", call. = FALSE)
}
if (any(data$total <= 0) || any(data$cultured < 0) || any(data$cultured > data$total)) {
  stop("cultured and total must satisfy 0 <= cultured <= total and total > 0", call. = FALSE)
}
if (any(abs(data$cultured - round(data$cultured)) > 1e-8) ||
    any(abs(data$total - round(data$total)) > 1e-8)) {
  stop("cultured and total must contain integer counts", call. = FALSE)
}
if (any(data$ci_low < 0 | data$ci_high > 1 | data$ci_low > data$ci_high)) {
  stop("confidence intervals must satisfy 0 <= ci_low <= ci_high <= 1", call. = FALSE)
}
if (any(data$q_value < 0 | data$q_value > 1) ||
    any(data$omnibus_p < 0 | data$omnibus_p > 1) ||
    any(data$omnibus_q < 0 | data$omnibus_q > 1)) {
  stop("p and q values must lie in [0, 1]", call. = FALSE)
}
if (anyDuplicated(data[c("dataset", "order")])) {
  stop("culture fraction input contains duplicate dataset/order rows", call. = FALSE)
}

data$proportion <- data$cultured / data$total
tolerance <- sqrt(.Machine$double.eps)
if (any(data$ci_low > data$proportion + tolerance | data$ci_high < data$proportion - tolerance)) {
  stop("each confidence interval must contain the cultured proportion", call. = FALSE)
}

dataset_order <- unique(data$dataset)
order_order <- unique(data$order)
for (dataset_name in dataset_order) {
  frame <- data[data$dataset == dataset_name, , drop = FALSE]
  if (length(unique(frame$omnibus_p)) != 1L || length(unique(frame$omnibus_q)) != 1L) {
    stop("omnibus_p and omnibus_q must be constant within each dataset", call. = FALSE)
  }
}

format_probability <- function(value) {
  if (value < 0.001) format(value, scientific = TRUE, digits = 2) else formatC(value, format = "f", digits = 3)
}

panel_rows <- do.call(rbind, lapply(dataset_order, function(dataset_name) {
  frame <- data[data$dataset == dataset_name, , drop = FALSE]
  data.frame(
    dataset = dataset_name,
    panel = paste0(
      dataset_name, "\nOverall order effect\nP=", format_probability(frame$omnibus_p[[1]]),
      "; BH q=", format_probability(frame$omnibus_q[[1]])
    ),
    stringsAsFactors = FALSE
  )
}))
panel_lookup <- stats::setNames(panel_rows$panel, panel_rows$dataset)

data$panel <- factor(unname(panel_lookup[data$dataset]), levels = panel_rows$panel)
data$order <- factor(data$order, levels = rev(order_order))
data$count_label <- paste0(round(data$cultured), "/", round(data$total))
data$significance <- ifelse(
  data$q_value < 0.001, "***",
  ifelse(data$q_value < 0.01, "**", ifelse(data$q_value < 0.05, "*", ""))
)

known_order_colors <- c(
  Actinomycetales = "#5B7DB1", Bacteroidales = "#6D8FBD",
  Burkholderiales = "#BFD7ED", Christensenellales = "#E9C08A",
  Coriobacteriales = "#B96B60", Enterobacterales = "#C95F48",
  Erysipelotrichales = "#D8BC78", Fusobacteriales = "#C98A63",
  Lachnospirales = "#D9CEE9", Lactobacillales = "#78BDB6",
  Oscillospirales = "#7EA66F", Peptostreptococcales = "#9BCB84",
  Veillonellales = "#B79A43", Other = "#9AA1A9"
)
unknown_orders <- setdiff(order_order, names(known_order_colors))
if (length(unknown_orders) > 0L) {
  fallback <- grDevices::hcl.colors(max(3L, length(unknown_orders)), "Dark 3")[seq_along(unknown_orders)]
  known_order_colors <- c(known_order_colors, stats::setNames(fallback, unknown_orders))
}
order_colors <- known_order_colors[order_order]

cultured_color <- "#8FBBD8"
uncultured_color <- "#C9E29A"
text_color <- "#20262E"

plot <- ggplot2::ggplot(data, ggplot2::aes(y = order)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = proportion, yend = order),
    color = cultured_color, linewidth = 3.2, alpha = 0.48, lineend = "butt"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = proportion, xend = 1, yend = order),
    color = uncultured_color, linewidth = 3.2, alpha = 0.48, lineend = "butt"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = ci_low, xend = ci_high, yend = order),
    color = cultured_color, linewidth = 1.0, lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = proportion, color = order),
    shape = 21, fill = cultured_color, size = 3.8, stroke = 1.1
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 1.03, label = count_label),
    hjust = 0, size = 3.5, color = text_color
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 1.27, label = significance),
    hjust = 0.5, size = 4.2, fontface = "bold", color = text_color
  ) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::scale_color_manual(values = order_colors, guide = "none", drop = FALSE) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.36), expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::labs(
    title = args$title, subtitle = args$subtitle,
    x = "Cultured proportion (95% CI)", y = NULL,
    caption = "n/N denotes cultured/total. Stars denote order-level q values: * <0.05, ** <0.01, *** <0.001."
  ) +
  ggplot2::theme_minimal(base_family = "sans", base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", color = text_color),
    plot.subtitle = ggplot2::element_text(color = "#56616C"),
    plot.caption = ggplot2::element_text(color = "#56616C", hjust = 0),
    strip.text = ggplot2::element_text(face = "bold", color = text_color, lineheight = 1.05),
    strip.background = ggplot2::element_rect(fill = "#F4F6F5", color = "#D9E2E1", linewidth = 0.4),
    axis.text = ggplot2::element_text(color = text_color),
    axis.title.x = ggplot2::element_text(face = "bold", color = text_color),
    panel.grid.major.y = ggplot2::element_line(color = "#EEF1F2", linewidth = 0.35),
    panel.grid.major.x = ggplot2::element_line(color = "#E1E6E8", linewidth = 0.4),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    plot.margin = ggplot2::margin(12, 18, 12, 12)
  )

height <- max(5.8, 2.8 + 0.42 * length(order_order))
save_ggplot_bundle(plot, args$output_prefix, args$formats, width = 11.8, height = height)
