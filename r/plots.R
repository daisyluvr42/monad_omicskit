#!/usr/bin/env Rscript
# Publication figures for expression analyses: volcano, heatmap, venn, PCA.
#
# params:
#   type          volcano | heatmap | venn | pca
#   output_name   basename for the saved figure
#
#   volcano: deg_path/deg records with gene, log2FoldChange, padj (or pvalue)
#            log2fc, padj thresholds; label_top; label_genes
#   heatmap: matrix_path/matrix; genes (subset); coldata_path/coldata + annotation_columns
#            scale_rows (default TRUE); show_rownames
#   venn:    sets = named list of character vectors (2-4 sets)
#   pca:     matrix_path/matrix; coldata_path/coldata; colour_column; label_samples

OMICS_R_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(OMICS_R_DIR, "lib", "common.R"))

plot_volcano <- function(params) {
  omics_require(c("ggplot2"))
  df <- omics_read_table(params, "deg_path", "deg")
  gene_col <- as.character(params$gene_column %||% "gene")[1]
  if (!gene_col %in% colnames(df)) gene_col <- colnames(df)[1]
  p_col <- if ("padj" %in% colnames(df)) "padj" else "pvalue"
  if (!"log2FoldChange" %in% colnames(df) || !p_col %in% colnames(df)) {
    stop("DEG table needs log2FoldChange and padj (or pvalue) columns.", call. = FALSE)
  }

  lfc_cut <- as.numeric(params$log2fc %||% 1)
  p_cut <- as.numeric(params$padj %||% 0.05)
  df$gene <- as.character(df[[gene_col]])
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  df$pval <- as.numeric(df[[p_col]])
  df <- df[is.finite(df$log2FoldChange) & is.finite(df$pval), ]
  df$pval[df$pval <= 0] <- .Machine$double.xmin
  df$direction <- ifelse(df$pval < p_cut & df$log2FoldChange >= lfc_cut, "Up",
                  ifelse(df$pval < p_cut & df$log2FoldChange <= -lfc_cut, "Down", "NS"))
  df$direction <- factor(df$direction, levels = c("Down", "NS", "Up"))

  label_genes <- as.character(params$label_genes %||% character())
  label_top <- as.integer(params$label_top %||% 10)
  sig <- df[df$direction != "NS", ]
  if (length(label_genes)) {
    labelled <- df[df$gene %in% label_genes, ]
  } else if (label_top > 0 && nrow(sig)) {
    ordered <- sig[order(sig$pval, -abs(sig$log2FoldChange)), ]
    labelled <- utils::head(ordered, label_top)
  } else {
    labelled <- df[0, ]
  }

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = -log10(pval), colour = direction)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::scale_colour_manual(values = c(Down = "#2f5f8f", NS = "#bdbdbd", Up = "#b5483a")) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "#666666", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", colour = "#666666", linewidth = 0.4) +
    ggplot2::labs(
      x = expression(log[2]~fold~change),
      y = bquote(-log[10]~.(p_col)),
      colour = NULL,
      title = as.character(params$title %||% "")[1]
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (nrow(labelled)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot + ggrepel::geom_text_repel(
        data = labelled, ggplot2::aes(label = gene),
        size = 3, max.overlaps = 30, show.legend = FALSE, colour = "#222222"
      )
    } else {
      plot <- plot + ggplot2::geom_text(
        data = labelled, ggplot2::aes(label = gene),
        size = 3, vjust = -0.6, show.legend = FALSE, colour = "#222222"
      )
    }
  }

  name <- as.character(params$output_name %||% "volcano")[1]
  figure <- omics_save_plot(plot, "figures", name, width = as.numeric(params$width %||% 6.5), height = as.numeric(params$height %||% 5.5))
  list(
    type = "volcano",
    genes_plotted = nrow(df),
    up = sum(df$direction == "Up"),
    down = sum(df$direction == "Down"),
    thresholds = list(log2fc = lfc_cut, p = p_cut, p_column = p_col),
    labelled_genes = labelled$gene,
    figure = figure
  )
}

plot_heatmap <- function(params) {
  omics_require(c("pheatmap"))
  mat <- omics_read_matrix(params)
  genes <- params$genes
  if (!is.null(genes)) {
    genes <- unique(as.character(genes))
    present <- intersect(genes, rownames(mat))
    if (!length(present)) {
      stop("None of the requested genes are present in the matrix. Check identifier type and species.", call. = FALSE)
    }
    mat <- mat[present, , drop = FALSE]
  }
  if (nrow(mat) < 2L) stop("Heatmap needs at least 2 genes.", call. = FALSE)

  # Drop zero-variance rows: they break row scaling and render as blank bands.
  variances <- apply(mat, 1, stats::var, na.rm = TRUE)
  dropped <- rownames(mat)[!is.finite(variances) | variances == 0]
  mat <- mat[is.finite(variances) & variances > 0, , drop = FALSE]
  if (nrow(mat) < 2L) stop("Fewer than 2 genes have non-zero variance across samples.", call. = FALSE)

  annotation <- NA
  ann_cols <- params$annotation_columns
  if (!is.null(params$coldata) || !is.null(params$coldata_path)) {
    coldata <- omics_read_table(params, "coldata_path", "coldata")
    rownames(coldata) <- as.character(coldata[[1]])
    coldata <- coldata[intersect(colnames(mat), rownames(coldata)), , drop = FALSE]
    if (!is.null(ann_cols)) {
      ann_cols <- intersect(as.character(ann_cols), colnames(coldata))
      if (length(ann_cols)) annotation <- coldata[, ann_cols, drop = FALSE]
    }
  }

  scale_rows <- if (is.null(params$scale_rows)) TRUE else isTRUE(params$scale_rows)
  name <- as.character(params$output_name %||% "heatmap")[1]
  show_rownames <- if (is.null(params$show_rownames)) nrow(mat) <= 60 else isTRUE(params$show_rownames)

  draw <- function() {
    pheatmap::pheatmap(
      mat,
      scale = if (scale_rows) "row" else "none",
      annotation_col = annotation,
      show_rownames = show_rownames,
      show_colnames = ncol(mat) <= 40,
      cluster_rows = nrow(mat) > 2,
      cluster_cols = isTRUE(params$cluster_columns %||% TRUE) && ncol(mat) > 2,
      color = grDevices::colorRampPalette(c("#2f5f8f", "#f7f7f7", "#b5483a"))(100),
      border_color = NA,
      main = as.character(params$title %||% "")[1],
      silent = TRUE
    )
  }
  height <- max(4, min(14, 0.16 * nrow(mat) + 2))
  figure <- omics_save_base_plot(
    function() grid::grid.draw(draw()$gtable),
    "figures", name,
    width = as.numeric(params$width %||% 7.5), height = as.numeric(params$height %||% height)
  )

  list(
    type = "heatmap",
    genes = nrow(mat),
    samples = ncol(mat),
    scaled = if (scale_rows) "row z-score" else "none",
    dropped_zero_variance = dropped,
    figure = figure,
    notes = if (scale_rows) "Colours are row z-scores, so they show relative pattern across samples, not absolute expression." else NULL
  )
}

plot_venn <- function(params) {
  sets <- params$sets
  if (is.null(sets) || length(sets) < 2L) stop("Provide `sets`: a named list of 2-4 gene vectors.", call. = FALSE)
  if (length(sets) > 4L) stop("Venn diagrams support at most 4 sets; use an UpSet plot for more.", call. = FALSE)
  sets <- lapply(sets, function(x) unique(as.character(x)))

  name <- as.character(params$output_name %||% "venn")[1]
  if (requireNamespace("ggvenn", quietly = TRUE)) {
    omics_require("ggplot2")
    plot <- ggvenn::ggvenn(
      sets, fill_color = c("#2f5f8f", "#b5483a", "#5f8f2f", "#8f6f2f")[seq_along(sets)],
      stroke_size = 0.4, set_name_size = 4, text_size = 3.5, show_percentage = FALSE
    ) + ggplot2::labs(title = as.character(params$title %||% "")[1])
    figure <- omics_save_plot(plot, "figures", name, width = as.numeric(params$width %||% 6), height = as.numeric(params$height %||% 5))
  } else {
    omics_require("ggvenn")
    figure <- NULL
  }

  # Every region, so the numbers in the figure can be traced back to gene lists.
  members <- unique(unlist(sets))
  membership <- vapply(sets, function(s) members %in% s, logical(length(members)))
  membership <- matrix(membership, nrow = length(members), dimnames = list(members, names(sets)))
  key <- apply(membership, 1, function(r) paste(names(sets)[r], collapse = " & "))
  regions <- split(members, key)

  intersection <- Reduce(intersect, sets)
  table_path <- omics_save_table(
    data.frame(gene = members, region = key, stringsAsFactors = FALSE),
    "figures", paste0(name, "_membership")
  )

  list(
    type = "venn",
    sets = lapply(sets, length),
    intersection_all = intersection,
    intersection_size = length(intersection),
    regions = lapply(regions, length),
    membership_table = table_path,
    figure = figure
  )
}

plot_pca <- function(params) {
  omics_require(c("ggplot2"))
  mat <- omics_read_matrix(params)
  variances <- apply(mat, 1, stats::var, na.rm = TRUE)
  mat <- mat[is.finite(variances) & variances > 0, , drop = FALSE]
  if (nrow(mat) < 2L) stop("Need at least 2 variable genes for PCA.", call. = FALSE)
  top_n <- as.integer(params$top_variable %||% 2000)
  if (nrow(mat) > top_n) {
    keep <- order(apply(mat, 1, stats::var, na.rm = TRUE), decreasing = TRUE)[seq_len(top_n)]
    mat <- mat[keep, , drop = FALSE]
  }

  pca <- stats::prcomp(t(mat), scale. = TRUE)
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  scores <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)

  colour_column <- as.character(params$colour_column %||% "")[1]
  if (nzchar(colour_column) && (!is.null(params$coldata) || !is.null(params$coldata_path))) {
    coldata <- omics_read_table(params, "coldata_path", "coldata")
    coldata$.sample <- as.character(coldata[[1]])
    if (colour_column %in% colnames(coldata)) {
      scores <- merge(scores, coldata[, c(".sample", colour_column)], by.x = "sample", by.y = ".sample", all.x = TRUE)
    } else {
      colour_column <- ""
    }
  } else {
    colour_column <- ""
  }

  plot <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2)) +
    (if (nzchar(colour_column)) ggplot2::geom_point(ggplot2::aes(colour = .data[[colour_column]]), size = 3, alpha = 0.85)
     else ggplot2::geom_point(size = 3, alpha = 0.85, colour = "#2f5f8f")) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", percent[1]),
      y = sprintf("PC2 (%.1f%%)", percent[2]),
      colour = if (nzchar(colour_column)) colour_column else NULL,
      title = as.character(params$title %||% "")[1]
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (isTRUE(params$label_samples)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3, max.overlaps = 30)
    } else {
      plot <- plot + ggplot2::geom_text(ggplot2::aes(label = sample), size = 3, vjust = -0.8)
    }
  }

  name <- as.character(params$output_name %||% "pca")[1]
  figure <- omics_save_plot(plot, "figures", name, width = as.numeric(params$width %||% 6.5), height = as.numeric(params$height %||% 5))

  list(
    type = "pca",
    samples = ncol(mat),
    genes_used = nrow(mat),
    variance_explained = list(PC1 = percent[1], PC2 = percent[2], PC3 = if (length(percent) >= 3) percent[3] else NULL),
    scores = scores,
    figure = figure,
    notes = "Inspect PCA for batch structure before interpreting group separation as biology."
  )
}

handler <- function(params) {
  type <- tolower(as.character(params$type %||% "volcano")[1])
  switch(
    type,
    volcano = plot_volcano(params),
    heatmap = plot_heatmap(params),
    venn = plot_venn(params),
    pca = plot_pca(params),
    stop(sprintf("Unknown type '%s'. Use volcano, heatmap, venn, or pca.", type), call. = FALSE)
  )
}

omics_main(handler)
