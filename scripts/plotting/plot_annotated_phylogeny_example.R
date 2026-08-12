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
  defaults = list(
    title = "Annotated phylogeny", formats = "pdf,svg,png", open_angle = "70",
    cohort_evidence = "", module_input = ""
  )
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
require_columns(annotation, c("tip_id", "group", "order", "status"), "tree annotation")
if (!"global_rank_score" %in% names(annotation)) {
  if (!"score" %in% names(annotation)) {
    stop("tree annotation is missing required column: global_rank_score", call. = FALSE)
  }
  annotation$global_rank_score <- annotation$score
}
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
annotation$global_rank_score <- suppressWarnings(as.numeric(annotation$global_rank_score))
if (any(!is.finite(annotation$global_rank_score)) ||
    any(annotation$global_rank_score < 0 | annotation$global_rank_score > 1)) {
  stop("tree annotation global_rank_score must lie in [0, 1]", call. = FALSE)
}
if (!all(annotation$status %in% c("Cultured", "Uncultured"))) {
  stop("tree annotation status must contain only Cultured or Uncultured", call. = FALSE)
}

validate_exact_tip_labels <- function(frame, table_name) {
  missing <- setdiff(tree$tip.label, frame$tip_id)
  extra <- setdiff(frame$tip_id, tree$tip.label)
  if (length(missing) > 0L || length(extra) > 0L) {
    stop(
      table_name, " and tree labels do not match; missing=", length(missing),
      ", extra=", length(extra), call. = FALSE
    )
  }
}

cohort_evidence <- NULL
if (nzchar(args$cohort_evidence)) {
  cohort_evidence <- read_table_file(args$cohort_evidence)
  require_columns(cohort_evidence, c("tip_id", "cohort_evidence_score"), "cohort evidence")
  cohort_evidence <- require_nonempty_text(cohort_evidence, "tip_id", "cohort evidence")
  if (anyDuplicated(cohort_evidence$tip_id)) {
    stop("cohort evidence contains duplicate tip_id values", call. = FALSE)
  }
  validate_exact_tip_labels(cohort_evidence, "cohort evidence")
  cohort_evidence$cohort_evidence_score <- suppressWarnings(as.numeric(cohort_evidence$cohort_evidence_score))
  if (any(!is.finite(cohort_evidence$cohort_evidence_score)) ||
      any(cohort_evidence$cohort_evidence_score < 0 | cohort_evidence$cohort_evidence_score > 1)) {
    stop("cohort evidence score must lie in [0, 1]", call. = FALSE)
  }
  cohort_evidence <- cohort_evidence[match(tree$tip.label, cohort_evidence$tip_id), , drop = FALSE]
}

modules <- NULL
module_names <- character()
if (nzchar(args$module_input)) {
  modules <- read_table_file(args$module_input)
  require_columns(modules, c("tip_id", "module", "score"), "module input")
  modules <- require_nonempty_text(modules, c("tip_id", "module"), "module input")
  if (anyDuplicated(modules[c("tip_id", "module")])) {
    stop("module input contains duplicate tip_id/module keys", call. = FALSE)
  }
  modules$score <- suppressWarnings(as.numeric(modules$score))
  if (any(!is.finite(modules$score)) || any(modules$score < 0 | modules$score > 1)) {
    stop("module input score must lie in [0, 1]", call. = FALSE)
  }
  module_names <- unique(modules$module)
  validate_exact_tip_labels(unique(modules[c("tip_id")]), "module input")
  expected_rows <- length(tree$tip.label) * length(module_names)
  if (nrow(modules) != expected_rows ||
      any(table(modules$tip_id) != length(module_names)) ||
      any(table(modules$module) != length(tree$tip.label))) {
    stop("module input must contain a complete tip-by-module matrix", call. = FALSE)
  }
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
numeric_ring_width <- 0.16
ring_gap <- 0.018
status_data <- make_sector_data(tip_data, 0, tip_radius + 0.02, "status", "status", layout$angle_width)
next_ring_inner <- tip_radius + 0.05

group_data <- make_sector_data(
  tip_data, next_ring_inner, next_ring_inner + ring_width,
  "group", "group", layout$angle_width
)
next_ring_inner <- next_ring_inner + ring_width + ring_gap
order_data <- make_sector_data(
  tip_data, next_ring_inner, next_ring_inner + ring_width,
  "order", "order", layout$angle_width
)
next_ring_inner <- next_ring_inner + ring_width + ring_gap
global_rank_data <- make_sector_data(
  tip_data, next_ring_inner, next_ring_inner + numeric_ring_width,
  "global_rank_score", "global_rank", layout$angle_width
)
next_ring_inner <- next_ring_inner + numeric_ring_width + ring_gap

cohort_evidence_data <- NULL
if (!is.null(cohort_evidence)) {
  evidence_tip_data <- merge(
    layout$tips, cohort_evidence[c("tip_id", "cohort_evidence_score")],
    by = "tip_id", sort = FALSE
  )
  evidence_tip_data <- evidence_tip_data[match(tree$tip.label, evidence_tip_data$tip_id), , drop = FALSE]
  cohort_evidence_data <- make_sector_data(
    evidence_tip_data, next_ring_inner, next_ring_inner + numeric_ring_width,
    "cohort_evidence_score", "cohort_evidence", layout$angle_width
  )
  next_ring_inner <- next_ring_inner + numeric_ring_width + ring_gap
}

module_data <- NULL
if (!is.null(modules)) {
  module_pieces <- vector("list", length(module_names))
  for (module_index in seq_along(module_names)) {
    module_name <- module_names[[module_index]]
    module_frame <- modules[modules$module == module_name, c("tip_id", "score"), drop = FALSE]
    module_frame <- merge(layout$tips, module_frame, by = "tip_id", sort = FALSE)
    module_frame <- module_frame[match(tree$tip.label, module_frame$tip_id), , drop = FALSE]
    module_pieces[[module_index]] <- make_sector_data(
      module_frame, next_ring_inner, next_ring_inner + ring_width,
      "score", paste0("module_", module_index), layout$angle_width
    )
    module_pieces[[module_index]]$module <- module_name
    next_ring_inner <- next_ring_inner + ring_width + ring_gap
  }
  module_data <- do.call(rbind, module_pieces)
}

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
    data = global_rank_data,
    ggplot2::aes(x = x, y = y, group = polygon, fill = as.numeric(value)),
    color = "#F1F3F5", linewidth = 0.03
  ) +
  ggplot2::scale_fill_gradient(
    name = "Global ML rank", low = "#F7F5FB", high = "#7655A3", limits = c(0, 1),
    guide = ggplot2::guide_colorbar(order = 4)
  )

if (!is.null(cohort_evidence_data)) {
  plot <- plot +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_polygon(
      data = cohort_evidence_data,
      ggplot2::aes(x = x, y = y, group = polygon, fill = as.numeric(value)),
      color = "#F1F3F5", linewidth = 0.03
    ) +
    ggplot2::scale_fill_gradient(
      name = "Cohort ML evidence", low = "#FFF4E6", high = "#D55E00", limits = c(0, 1),
      guide = ggplot2::guide_colorbar(order = 5)
    )
}

if (!is.null(module_data)) {
  plot <- plot +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_polygon(
      data = module_data,
      ggplot2::aes(x = x, y = y, group = polygon, fill = as.numeric(value)),
      color = "#F1F3F5", linewidth = 0.03
    ) +
    ggplot2::scale_fill_gradient(
      name = "Functional module score", low = "#EFF8F2", high = "#238B45", limits = c(0, 1),
      guide = ggplot2::guide_colorbar(order = 6)
    )
}

ring_labels <- c("Group", "Order", "Global ML rank")
if (!is.null(cohort_evidence_data)) ring_labels <- c(ring_labels, "Cohort ML evidence")
if (length(module_names) > 0L) ring_labels <- c(ring_labels, module_names)

plot <- plot +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::labs(
    title = args$title,
    subtitle = paste("Rings (inner to outer):", paste(ring_labels, collapse = "  |  "))
  ) +
  ggplot2::theme_void(base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "#20262E"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#4D565F", size = 9),
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.key.size = grid::unit(0.42, "cm"),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

save_ggplot_bundle(plot, args$output_prefix, args$formats, width = 11.2, height = 8.8)
