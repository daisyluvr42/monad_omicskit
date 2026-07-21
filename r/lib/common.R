#!/usr/bin/env Rscript
# Shared helpers for Scholar R analysis scripts.
#
# Every Scholar R script is invoked as:
#   Rscript <script>.R <args.json> <result.json>
# and must write a JSON object to <result.json>. Figures follow the same
# {png, svg} contract as the Python side so Skills can treat both alike.

suppressPackageStartupMessages({
  library(jsonlite)
})

# base R gained %||% in 4.4; define it for older versions.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

omics_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) < 2L) {
    stop("Usage: Rscript <script>.R <args.json> <result.json>", call. = FALSE)
  }
  list(
    params = jsonlite::fromJSON(argv[[1]], simplifyVector = TRUE),
    out = argv[[2]]
  )
}

# Fail with an actionable message instead of an R namespace error.
omics_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      sprintf(
        "Missing R packages: %s. Install them with: Rscript r/bootstrap.R %s",
        paste(missing, collapse = ", "),
        paste(missing, collapse = " ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

omics_output_dir <- function(subdir = NULL) {
  base <- Sys.getenv("OMICS_OUTPUT_DIR", unset = "")
  if (!nzchar(base)) {
    base <- file.path(path.expand("~"), ".workbuddy", "workspace", "omics")
  }
  path <- if (is.null(subdir)) base else file.path(base, subdir)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Strip characters that would break a filename; keep it recognisable.
omics_safe_name <- function(value, default = "output") {
  value <- gsub("[^A-Za-z0-9._-]+", "-", as.character(value)[1])
  value <- gsub("^-+|-+$", "", value)
  if (!nzchar(value)) default else value
}

# Save a ggplot (or grid object) as 300 dpi PNG plus editable SVG.
omics_save_plot <- function(plot, subdir, name, width = 7, height = 5) {
  omics_require("ggplot2")
  dir <- omics_output_dir(subdir)
  name <- omics_safe_name(name)
  png_path <- file.path(dir, paste0(name, ".png"))
  svg_path <- file.path(dir, paste0(name, ".svg"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(svg_path, plot, width = width, height = height, bg = "white")
  list(png = png_path, svg = svg_path)
}

# For base-graphics output (maftools, WGCNA, forestplot, nomogram, ...).
omics_save_base_plot <- function(draw, subdir, name, width = 7, height = 5) {
  dir <- omics_output_dir(subdir)
  name <- omics_safe_name(name)
  png_path <- file.path(dir, paste0(name, ".png"))
  svg_path <- file.path(dir, paste0(name, ".svg"))
  grDevices::png(png_path, width = width, height = height, units = "in", res = 300, bg = "white")
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::svg(svg_path, width = width, height = height, bg = "white")
  draw()
  grDevices::dev.off()
  list(png = png_path, svg = svg_path)
}

# Write a result table next to the figures so users get the numbers, not just a picture.
omics_save_table <- function(df, subdir, name) {
  dir <- omics_output_dir(subdir)
  path <- file.path(dir, paste0(omics_safe_name(name), ".csv"))
  utils::write.csv(df, path, row.names = FALSE, na = "")
  path
}

# Read a matrix from an explicit path or from inline records.
# Bioinformatics matrices are gene-by-sample with gene IDs in the first column.
omics_read_matrix <- function(params, path_key = "matrix_path", records_key = "matrix") {
  path <- params[[path_key]]
  if (!is.null(path) && nzchar(path)) {
    path <- path.expand(path)
    if (!file.exists(path)) stop(sprintf("File not found: %s", path), call. = FALSE)
    sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
    df <- utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (!is.null(params[[records_key]])) {
    df <- as.data.frame(params[[records_key]], check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop(sprintf("Provide %s or %s", path_key, records_key), call. = FALSE)
  }
  if (ncol(df) < 2L) stop("Matrix needs an ID column plus at least one sample column", call. = FALSE)
  ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  # Collapse duplicated gene symbols by highest mean expression, the usual convention.
  if (anyDuplicated(ids)) {
    keep <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
    mat <- mat[keep, , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  mat
}

omics_read_table <- function(params, path_key, records_key) {
  path <- params[[path_key]]
  if (!is.null(path) && nzchar(path)) {
    path <- path.expand(path)
    if (!file.exists(path)) stop(sprintf("File not found: %s", path), call. = FALSE)
    sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
    return(utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (!is.null(params[[records_key]])) {
    return(as.data.frame(params[[records_key]], check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop(sprintf("Provide %s or %s", path_key, records_key), call. = FALSE)
}

# Guard against the single most common bioinformatics input error:
# feeding already-normalised or log-transformed values to a count-based model.
omics_assert_counts <- function(mat, method) {
  if (any(mat < 0, na.rm = TRUE)) {
    stop(sprintf("%s requires raw integer counts, but the matrix contains negative values (looks log-transformed or centred). Use limma instead.", method), call. = FALSE)
  }
  finite <- mat[is.finite(mat)]
  if (length(finite) && max(finite) < 50) {
    stop(sprintf("%s requires raw integer counts, but the maximum value is %.3f (looks log-transformed or TPM-scaled). Use limma instead.", method, max(finite)), call. = FALSE)
  }
  fractional <- abs(mat - round(mat)) > 1e-8
  if (mean(fractional, na.rm = TRUE) > 0.01) {
    stop(sprintf("%s requires integer counts, but the matrix is mostly non-integer (looks like TPM/FPKM/normalised data). Use limma instead.", method), call. = FALSE)
  }
  invisible(TRUE)
}

omics_species <- function(params, default = NULL) {
  species <- params$species
  if (is.null(species) || !nzchar(as.character(species)[1])) {
    if (!is.null(default)) return(default)
    stop("species is required (human or mouse); gene identifiers differ and must not be guessed.", call. = FALSE)
  }
  species <- tolower(as.character(species)[1])
  known <- list(
    human = list(species = "human", orgdb = "org.Hs.eg.db", kegg = "hsa", msigdb = "Homo sapiens", reactome = "human"),
    mouse = list(species = "mouse", orgdb = "org.Mm.eg.db", kegg = "mmu", msigdb = "Mus musculus", reactome = "mouse")
  )
  if (!species %in% names(known)) {
    stop(sprintf("Unsupported species '%s'. Use human or mouse.", species), call. = FALSE)
  }
  known[[species]]
}

# Fields that are conceptually lists must stay JSON arrays even when they hold a
# single element. auto_unbox would otherwise collapse them to a bare scalar, so
# the result's shape would depend on the data and every caller would need to
# handle both forms.
omics_arr <- function(x) {
  if (is.null(x)) return(NULL)
  I(unname(x))
}

omics_write <- function(result, out) {
  json <- jsonlite::toJSON(
    result,
    auto_unbox = TRUE,
    digits = NA,
    na = "null",
    null = "null",
    force = TRUE
  )
  writeLines(json, out, useBytes = TRUE)
  invisible(TRUE)
}

# Wrap a handler so R errors become a clean JSON payload the MCP layer can surface.
omics_main <- function(handler) {
  ctx <- omics_args()
  result <- tryCatch(
    handler(ctx$params),
    error = function(e) list(error = conditionMessage(e))
  )
  omics_write(result, ctx$out)
  if (!is.null(result$error)) quit(status = 3, save = "no")
  invisible(TRUE)
}
