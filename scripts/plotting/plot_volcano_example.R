#!/usr/bin/env Rscript

# Plot precomputed feature effects and FDR values as a labeled volcano plot.
# Representative source lineage: pipeline_Pan cross-cohort meta-analysis figures.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[[1]]
source(file.path(dirname(normalizePath(script_file)), "plotting_common.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing R package: ggplot2", call. = FALSE)
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Missing R package: ggrepel", call. = FALSE)
if (!requireNamespace("svglite", quietly = TRUE)) stop("Missing R package: svglite", call. = FALSE)

args <- parse_cli(
  commandArgs(trailingOnly = TRUE),
  required = c("input", "output_prefix"),
  defaults = list(
    title = "Feature association volcano plot",
    formats = "pdf,svg,png",
    effect_threshold = "0.5",
    q_threshold = "0.05",
    max_labels = "12",
    seed = "42"
  )
)
set.seed(as.integer(args$seed))
effect_threshold <- as.numeric(args$effect_threshold)
q_threshold <- as.numeric(args$q_threshold)
max_labels <- as.integer(args$max_labels)
if (!is.finite(effect_threshold) || effect_threshold < 0) stop("effect threshold must be non-negative", call. = FALSE)
if (!is.finite(q_threshold) || q_threshold <= 0 || q_threshold > 1) stop("q threshold must be in (0, 1]", call. = FALSE)
if (!is.finite(max_labels) || max_labels < 0) stop("max labels must be non-negative", call. = FALSE)

data <- read_table_file(args$input)
require_columns(data, c("feature", "effect", "q_value"), "volcano input")
data$effect <- suppressWarnings(as.numeric(data$effect))
data$q_value <- suppressWarnings(as.numeric(data$q_value))
if (any(!is.finite(data$effect)) || any(!is.finite(data$q_value))) {
  stop("effect and q_value must be finite", call. = FALSE)
}
if (any(data$q_value <= 0 | data$q_value > 1)) stop("q_value must be in (0, 1]", call. = FALSE)
if (anyDuplicated(data$feature)) stop("volcano input contains duplicate feature values", call. = FALSE)

data$minus_log10_q <- -log10(pmax(data$q_value, .Machine$double.xmin))
data$direction <- "Not significant"
data$direction[data$q_value <= q_threshold & data$effect >= effect_threshold] <- "Positive"
data$direction[data$q_value <= q_threshold & data$effect <= -effect_threshold] <- "Negative"

if (!"label" %in% names(data)) data$label <- ""
data$label <- as.character(data$label)
candidate_order <- order(data$q_value, -abs(data$effect))
selected <- head(candidate_order[data$direction[candidate_order] != "Not significant"], max_labels)
fallback_labels <- !nzchar(data$label)
data$label[fallback_labels] <- ""
data$label[selected] <- ifelse(nzchar(data$label[selected]), data$label[selected], data$feature[selected])
keep_label <- rep(FALSE, nrow(data))
keep_label[selected] <- TRUE
data$label[!keep_label] <- ""

palette <- c("Negative" = "#3B78A8", "Not significant" = "#B8B8B8", "Positive" = "#A85B50")
plot <- ggplot2::ggplot(data, ggplot2::aes(x = effect, y = minus_log10_q, color = direction)) +
  ggplot2::geom_hline(yintercept = -log10(q_threshold), linetype = "dashed", color = "#777777", linewidth = 0.45) +
  ggplot2::geom_vline(xintercept = c(-effect_threshold, effect_threshold), linetype = "dashed", color = "#777777", linewidth = 0.45) +
  ggplot2::geom_point(size = 2.0, alpha = 0.78) +
  ggrepel::geom_text_repel(
    data = data[nzchar(data$label), , drop = FALSE],
    ggplot2::aes(label = label),
    size = 2.8,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.35,
    seed = as.integer(args$seed),
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = palette, breaks = names(palette)) +
  ggplot2::labs(title = args$title, x = "Effect size", y = expression(-log[10](q)), color = NULL) +
  ggplot2::theme_bw(base_size = 11, base_family = "sans") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

save_ggplot_bundle(plot, args$output_prefix, args$formats, width = 7.2, height = 5.8)
