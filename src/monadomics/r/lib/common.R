#!/usr/bin/env Rscript
# Shared helpers for MonadOmics R analysis scripts.
#
# The CLI sets OMICS_R_DIR and invokes each script as:
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
        "Missing R packages: %s. Run monadomics doctor for installation instructions.",
        paste(missing, collapse = ", ")
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
  omics_require("svglite")
  dir <- omics_output_dir(subdir)
  name <- omics_safe_name(name)
  png_path <- file.path(dir, paste0(name, ".png"))
  svg_path <- file.path(dir, paste0(name, ".svg"))
  grDevices::png(png_path, width = width, height = height, units = "in", res = 300, bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  svglite::svglite(svg_path, width = width, height = height, bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
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
omics_matrix_type <- function(params) {
  matrix_type <- tolower(as.character(params$matrix_type %||% "")[1])
  if (!matrix_type %in% c("counts", "normalized")) {
    stop("matrix_type is required and must be 'counts' or 'normalized'.", call. = FALSE)
  }
  matrix_type
}

omics_read_matrix <- function(params, path_key = "matrix_path", records_key = "matrix", matrix_type = NULL) {
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
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("Gene IDs must be non-empty and unique; resolve duplicates explicitly before analysis.", call. = FALSE)
  }
  sample_names <- names(df)[-1]
  if (any(!nzchar(trimws(sample_names))) || anyDuplicated(sample_names)) {
    stop("Matrix sample IDs must be non-empty and unique.", call. = FALSE)
  }
  values <- df[, -1, drop = FALSE]
  converted <- lapply(values, function(column) suppressWarnings(as.numeric(column)))
  bad_columns <- names(values)[vapply(
    seq_along(values),
    function(index) any(!is.na(values[[index]]) & is.na(converted[[index]])),
    logical(1)
  )]
  if (length(bad_columns)) {
    stop(sprintf("Matrix contains non-numeric values in column(s): %s", paste(bad_columns, collapse = ", ")), call. = FALSE)
  }
  mat <- as.matrix(as.data.frame(converted, check.names = FALSE))
  if (anyNA(mat)) {
    stop("Matrix contains missing values; impute or remove them before analysis.", call. = FALSE)
  }
  if (any(!is.finite(mat))) stop("Matrix contains non-finite values.", call. = FALSE)
  rownames(mat) <- ids
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
  fractional <- abs(mat - round(mat)) > 1e-8
  if (any(fractional, na.rm = TRUE)) {
    stop(sprintf("%s requires raw integer counts, but the matrix contains non-integer values. Use limma with matrix_type='normalized' instead.", method), call. = FALSE)
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

# Wrap a handler so R errors become a clean JSON payload the CLI can surface.
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
