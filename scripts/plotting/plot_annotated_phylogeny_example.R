#!/usr/bin/env Rscript

# Plot a fan phylogeny with portable group, culture-status, and score annotations.
# Representative source lineage: pipeline_Pan annotated phylogeny figures.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[[1]]
source(file.path(dirname(normalizePath(script_file)), "plotting_common.R"))

if (!requireNamespace("ape", quietly = TRUE)) stop("Missing R package: ape", call. = FALSE)
if (!requireNamespace("svglite", quietly = TRUE)) stop("Missing R package: svglite", call. = FALSE)

args <- parse_cli(
  commandArgs(trailingOnly = TRUE),
  required = c("tree", "annotation", "output_prefix"),
  defaults = list(title = "Annotated phylogeny", formats = "pdf,svg,png")
)

tree <- ape::read.tree(args$tree)
if (is.null(tree) || length(tree$tip.label) < 3L) {
  stop("Tree must contain at least three tips", call. = FALSE)
}
annotation <- read_table_file(args$annotation)
require_columns(annotation, c("tip_id", "group", "status", "score"), "tree annotation")
annotation <- require_nonempty_text(annotation, c("tip_id", "group", "status"), "tree annotation")
if (anyDuplicated(annotation$tip_id)) stop("tree annotation contains duplicate tip_id values", call. = FALSE)
missing_tips <- setdiff(tree$tip.label, annotation$tip_id)
if (length(missing_tips) > 0L) {
  stop("tree annotation is missing tips: ", paste(head(missing_tips, 8L), collapse = ", "), call. = FALSE)
}
annotation <- annotation[match(tree$tip.label, annotation$tip_id), , drop = FALSE]
annotation$score <- suppressWarnings(as.numeric(annotation$score))
if (any(!is.finite(annotation$score)) || any(annotation$score < 0)) {
  stop("tree annotation score must contain finite, non-negative values", call. = FALSE)
}

groups <- unique(as.character(annotation$group))
statuses <- unique(as.character(annotation$status))
group_colors <- stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
status_colors <- stats::setNames(grDevices::hcl.colors(length(statuses), "Set 2"), statuses)
score_range <- range(annotation$score)
if (diff(score_range) == 0) {
  score_cex <- rep(1.15, nrow(annotation))
} else {
  score_cex <- 0.75 + 1.15 * (annotation$score - score_range[[1]]) / diff(score_range)
}

draw_tree <- function() {
  graphics::par(mar = c(1.5, 1.5, 3.2, 1.5), family = "sans")
  ape::plot.phylo(
    tree,
    type = "fan",
    tip.color = unname(group_colors[annotation$group]),
    edge.color = "#59636D",
    edge.width = 0.8,
    cex = 0.62,
    label.offset = 0.02,
    no.margin = FALSE
  )
  ape::tiplabels(
    pch = 21,
    bg = unname(status_colors[annotation$status]),
    col = "white",
    cex = score_cex,
    lwd = 0.5
  )
  graphics::title(main = args$title, font.main = 2, cex.main = 1.15)
  graphics::legend(
    "topleft",
    legend = names(group_colors),
    text.col = unname(group_colors),
    bty = "n",
    title = "Group",
    cex = 0.72
  )
  graphics::legend(
    "topright",
    legend = names(status_colors),
    pch = 21,
    pt.bg = unname(status_colors),
    col = "white",
    bty = "n",
    title = "Status (point size = score)",
    cex = 0.72
  )
}

save_base_bundle(draw_tree, args$output_prefix, args$formats, width = 8.2, height = 8.2)
