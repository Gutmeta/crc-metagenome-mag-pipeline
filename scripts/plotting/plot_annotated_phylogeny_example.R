#!/usr/bin/env Rscript

# Plot a fan phylogeny with culture-status sectors and portable annotation rings.
# Analysis context: phylogenomic annotation.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[[1]]
source(file.path(dirname(normalizePath(script_file)), "plotting_common.R"))

if (!requireNamespace("ape", quietly = TRUE)) stop("Missing R package: ape", call. = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing R package: ggplot2", call. = FALSE)
if (!requireNamespace("ggnewscale", quietly = TRUE)) stop("Missing R package: ggnewscale", call. = FALSE)
if (!requireNamespace("scales", quietly = TRUE)) stop("Missing R package: scales", call. = FALSE)
if (!requireNamespace("svglite", quietly = TRUE)) stop("Missing R package: svglite", call. = FALSE)

args <- parse_cli(
  commandArgs(trailingOnly = TRUE),
  required = c("tree", "annotation", "output_prefix"),
  defaults = list(title = "Annotated phylogeny", formats = "pdf,svg,png", open_angle = "70")
)
open_angle <- suppressWarnings(as.numeric(args$open_angle))
if (!is.finite(open_angle) || open_angle < 0 || open_angle >= 180) {
  stop("open angle must be in [0, 180)", call. = FALSE)
}

tree <- ape::read.tree(args$tree)
if (is.null(tree) || length(tree$tip.label) < 3L) {
  stop("Tree must contain at least three tips", call. = FALSE)
}
if (anyDuplicated(tree$tip.label)) stop("Tree contains duplicate tip labels", call. = FALSE)
tree <- ape::ladderize(tree, right = FALSE)

annotation <- read_table_file(args$annotation)
require_columns(annotation, c("tip_id", "group", "order", "status", "score"), "tree annotation")
annotation <- require_nonempty_text(annotation, c("tip_id", "group", "order", "status"), "tree annotation")
if (anyDuplicated(annotation$tip_id)) stop("tree annotation contains duplicate tip_id values", call. = FALSE)
missing_tips <- setdiff(tree$tip.label, annotation$tip_id)
extra_tips <- setdiff(annotation$tip_id, tree$tip.label)
if (length(missing_tips) > 0L || length(extra_tips) > 0L) {
  stop(
    "tree and annotation labels do not match; missing=", length(missing_tips),
    ", extra=", length(extra_tips), call. = FALSE
  )
}
annotation <- annotation[match(tree$tip.label, annotation$tip_id), , drop = FALSE]
annotation$score <- suppressWarnings(as.numeric(annotation$score))
if (any(!is.finite(annotation$score)) || any(annotation$score < 0 | annotation$score > 1)) {
  stop("tree annotation score must lie in [0, 1]", call. = FALSE)
}
if (!all(annotation$status %in% c("Cultured", "Uncultured"))) {
  stop("tree annotation status must contain only Cultured or Uncultured", call. = FALSE)
}

make_sector_data <- function(frame, inner_radius, outer_radius, value_column, prefix, angle_width) {
  pieces <- vector("list", nrow(frame))
  for (index in seq_len(nrow(frame))) {
    theta <- seq(
      frame$angle[[index]] - angle_width / 2,
      frame$angle[[index]] + angle_width / 2,
      length.out = 7
    )
    pieces[[index]] <- data.frame(
      polygon = paste0(prefix, "_", index),
      x = c(outer_radius * cos(theta), rev(inner_radius * cos(theta))),
      y = c(outer_radius * sin(theta), rev(inner_radius * sin(theta))),
      value = frame[[value_column]][[index]],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pieces)
}

make_tree_layout <- function(tree, open_angle, inner_radius, band_width) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ape::plot.phylo(
    tree, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
    open.angle = open_angle, no.margin = TRUE
  )
  plot_environment <- get(".PlotPhyloEnv", envir = asNamespace("ape"))
  layout <- get("last_plot.phylo", envir = plot_environment)
  raw_radius <- sqrt(layout$xx^2 + layout$yy^2)
  raw_angle <- atan2(layout$yy, layout$xx)
  radius_scale <- max(raw_radius)
  radius <- inner_radius + raw_radius / radius_scale * band_width
  coordinates <- data.frame(
    node = seq_along(layout$xx), angle = raw_angle, radius = radius,
    x = radius * cos(raw_angle), y = radius * sin(raw_angle)
  )
  edges <- data.frame(
    parent = layout$edge[, 1], child = layout$edge[, 2],
    parent_radius = coordinates$radius[layout$edge[, 1]],
    child_radius = coordinates$radius[layout$edge[, 2]],
    child_angle = coordinates$angle[layout$edge[, 2]]
  )
  edges$x <- edges$parent_radius * cos(edges$child_angle)
  edges$y <- edges$parent_radius * sin(edges$child_angle)
  edges$xend <- edges$child_radius * cos(edges$child_angle)
  edges$yend <- edges$child_radius * sin(edges$child_angle)

  arc_pieces <- vector("list", nrow(edges))
  for (index in seq_len(nrow(edges))) {
    parent_angle <- coordinates$angle[edges$parent[[index]]]
    child_angle <- coordinates$angle[edges$child[[index]]]
    delta <- child_angle - parent_angle
    if (delta > pi) delta <- delta - 2 * pi
    if (delta < -pi) delta <- delta + 2 * pi
    theta <- seq(parent_angle, parent_angle + delta, length.out = 10)
    arc_pieces[[index]] <- data.frame(
      arc = paste0("arc_", index),
      x = edges$parent_radius[[index]] * cos(theta),
      y = edges$parent_radius[[index]] * sin(theta)
    )
  }
  tips <- coordinates[seq_along(tree$tip.label), , drop = FALSE]
  tips$tip_id <- tree$tip.label
  angle_differences <- diff(sort(tips$angle))
  angle_differences <- angle_differences[angle_differences > 1e-4]
  angle_width <- stats::median(angle_differences, na.rm = TRUE) * 0.90
  if (!is.finite(angle_width) || angle_width <= 0) {
    stop("Could not determine angular spacing between tree tips", call. = FALSE)
  }
  list(edges = edges, arcs = do.call(rbind, arc_pieces), tips = tips, angle_width = angle_width)
}

tree_inner_radius <- 0.28
tree_band_width <- 1.15
layout <- make_tree_layout(tree, open_angle, tree_inner_radius, tree_band_width)
tip_radius <- tree_inner_radius + tree_band_width
tip_data <- merge(layout$tips, annotation, by = "tip_id", sort = FALSE)
tip_data <- tip_data[match(tree$tip.label, tip_data$tip_id), , drop = FALSE]

ring_width <- 0.13
ring_gap <- 0.018
status_data <- make_sector_data(tip_data, 0, tip_radius + 0.02, "status", "status", layout$angle_width)
group_inner <- tip_radius + 0.05
group_data <- make_sector_data(tip_data, group_inner, group_inner + ring_width, "group", "group", layout$angle_width)
order_inner <- group_inner + ring_width + ring_gap
order_data <- make_sector_data(tip_data, order_inner, order_inner + ring_width, "order", "order", layout$angle_width)
score_inner <- order_inner + ring_width + 0.045
score_data <- make_sector_data(tip_data, score_inner, score_inner + 0.18, "score", "score", layout$angle_width)

status_palette <- c(Cultured = "#8FBBD8", Uncultured = "#C9E29A")
groups <- unique(annotation$group)
group_palette <- stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
orders <- unique(annotation$order)
known_order_palette <- c(
  Actinomycetales = "#5B7DB1", Bacteroidales = "#6D8FBD",
  Burkholderiales = "#BFD7ED", Christensenellales = "#E9C08A",
  Coriobacteriales = "#B96B60", Enterobacterales = "#C95F48",
  Erysipelotrichales = "#D8BC78", Fusobacteriales = "#C98A63",
  Lachnospirales = "#D9CEE9", Lactobacillales = "#78BDB6",
  Oscillospirales = "#7EA66F", Peptostreptococcales = "#9BCB84",
  Veillonellales = "#B79A43", Other = "#9AA1A9"
)
unknown_orders <- setdiff(orders, names(known_order_palette))
if (length(unknown_orders) > 0L) {
  colors <- grDevices::hcl.colors(max(3L, length(unknown_orders)), "Dark 3")[seq_along(unknown_orders)]
  known_order_palette <- c(known_order_palette, stats::setNames(colors, unknown_orders))
}
order_palette <- known_order_palette[orders]

plot <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = status_data,
    ggplot2::aes(x = x, y = y, group = polygon, fill = value),
    color = "white", linewidth = 0.01, alpha = 0.48
  ) +
  ggplot2::scale_fill_manual(
    name = "Status", values = status_palette, drop = FALSE,
    guide = ggplot2::guide_legend(order = 3, override.aes = list(alpha = 1))
  ) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_path(
    data = layout$arcs,
    ggplot2::aes(x = x, y = y, group = arc),
    color = "#4D565F", linewidth = 0.34, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = layout$edges,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    color = "#4D565F", linewidth = 0.34, lineend = "round"
  ) +
  ggplot2::geom_polygon(
    data = group_data,
    ggplot2::aes(x = x, y = y, group = polygon, fill = value),
    color = "#F1F3F5", linewidth = 0.03
  ) +
  ggplot2::scale_fill_manual(
    name = "Group", values = group_palette, drop = FALSE,
    guide = ggplot2::guide_legend(order = 1)
  ) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_polygon(
    data = order_data,
    ggplot2::aes(x = x, y = y, group = polygon, fill = value),
    color = "#F1F3F5", linewidth = 0.03
  ) +
  ggplot2::scale_fill_manual(
    name = "Order", values = order_palette, drop = FALSE,
    guide = ggplot2::guide_legend(order = 2)
  ) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_polygon(
    data = score_data,
    ggplot2::aes(x = x, y = y, group = polygon, fill = as.numeric(value)),
    color = "#F1F3F5", linewidth = 0.03
  ) +
  ggplot2::scale_fill_gradient(
    name = "Evidence score", low = "#F7F5FB", high = "#8C6BC4", limits = c(0, 1),
    guide = ggplot2::guide_colorbar(order = 4)
  ) +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::labs(title = args$title) +
  ggplot2::theme_void(base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "#20262E"),
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.key.size = grid::unit(0.42, "cm"),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

save_ggplot_bundle(plot, args$output_prefix, args$formats, width = 10.2, height = 8.4)
