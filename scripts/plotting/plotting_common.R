#!/usr/bin/env Rscript

# Shared, privacy-safe helpers for the portable R plotting examples.

parse_cli <- function(args, required = character(), defaults = list()) {
  values <- defaults
  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, call. = FALSE)
    }
    if (grepl("=", token, fixed = TRUE)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      value <- sub("^--[^=]+=", "", token)
    } else {
      key <- sub("^--", "", token)
      index <- index + 1L
      if (index > length(args)) {
        stop("Missing value for --", key, call. = FALSE)
      }
      value <- args[[index]]
    }
    values[[gsub("-", "_", key, fixed = TRUE)]] <- value
    index <- index + 1L
  }
  missing <- required[!vapply(required, function(key) {
    !is.null(values[[key]]) && nzchar(values[[key]])
  }, logical(1))]
  if (length(missing) > 0L) {
    stop("Missing required arguments: ", paste0("--", gsub("_", "-", missing), collapse = ", "), call. = FALSE)
  }
  values
}

read_table_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input table does not exist: ", path, call. = FALSE)
  }
  if (tolower(tools::file_ext(path)) == "csv") {
    frame <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    frame <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (nrow(frame) == 0L) {
    stop("Input table is empty: ", path, call. = FALSE)
  }
  frame
}

require_columns <- function(frame, columns, table_name = "input") {
  missing <- setdiff(columns, names(frame))
  if (length(missing) > 0L) {
    stop(table_name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

parse_formats <- function(value) {
  formats <- unique(tolower(trimws(strsplit(value, ",", fixed = TRUE)[[1]])))
  formats <- formats[nzchar(formats)]
  allowed <- c("pdf", "svg", "png")
  unknown <- setdiff(formats, allowed)
  if (length(formats) == 0L || length(unknown) > 0L) {
    stop("Formats must be one or more of: pdf, svg, png", call. = FALSE)
  }
  formats
}

save_ggplot_bundle <- function(plot, output_prefix, formats, width, height) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  for (format in parse_formats(formats)) {
    path <- paste0(output_prefix, ".", format)
    if (format == "svg") {
      ggplot2::ggsave(path, plot = plot, width = width, height = height, device = svglite::svglite)
    } else if (format == "pdf") {
      ggplot2::ggsave(path, plot = plot, width = width, height = height, device = grDevices::cairo_pdf)
    } else {
      ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 300, bg = "white")
    }
    if (!file.exists(path) || file.info(path)$size <= 0) {
      stop("Figure was not written correctly: ", path, call. = FALSE)
    }
  }
}

save_base_bundle <- function(draw_function, output_prefix, formats, width, height) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  for (format in parse_formats(formats)) {
    path <- paste0(output_prefix, ".", format)
    if (format == "pdf") {
      grDevices::cairo_pdf(path, width = width, height = height)
    } else if (format == "svg") {
      svglite::svglite(path, width = width, height = height)
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = 300, type = "cairo")
    }
    draw_function()
    grDevices::dev.off()
    if (!file.exists(path) || file.info(path)$size <= 0) {
      stop("Figure was not written correctly: ", path, call. = FALSE)
    }
  }
}

format_q_value <- function(value) {
  if (!is.finite(value)) return("q = NA")
  if (value < 0.001) return("q < 0.001")
  sprintf("q = %.3f", value)
}
