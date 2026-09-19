#!/usr/bin/env Rscript
# Differential expression: DESeq2, edgeR, or limma.
#
# params:
#   method        deseq2 | edger | limma
#   matrix_path   gene-by-sample matrix, gene IDs in column 1 (or `matrix` records)
#   coldata_path  sample metadata, sample IDs in column 1 (or `coldata` records)
#   group_column  column in coldata holding the comparison groups
#   treat, control   group levels; log2FC is treat vs control
#   covariates    optional character vector of additional coldata columns
#   log2fc, padj  significance thresholds (default 1, 0.05)
#   padj_method   BH (default), bonferroni, ...
#   voom          limma only: TRUE when the input is raw counts
#   output_name   basename for saved tables

OMICS_R_DIR <- Sys.getenv("OMICS_R_DIR")
source(file.path(OMICS_R_DIR, "lib", "common.R"))

build_design <- function(group_column, covariates) {
  stats::reformulate(c(covariates, group_column))
}

align_inputs <- function(params) {
  matrix_type <- omics_matrix_type(params)
  mat <- omics_read_matrix(params, matrix_type = matrix_type)
  coldata <- omics_read_table(params, "coldata_path", "coldata")
  rownames(coldata) <- as.character(coldata[[1]])

  group_column <- params$group_column
  if (is.null(group_column) || !nzchar(group_column)) {
    stop("group_column is required", call. = FALSE)
  }
  if (!group_column %in% colnames(coldata)) {
    stop(sprintf("group_column '%s' not found in coldata columns: %s",
                 group_column, paste(colnames(coldata), collapse = ", ")), call. = FALSE)
  }

  shared <- intersect(colnames(mat), rownames(coldata))
  if (length(shared) < 4L) {
    stop(sprintf("Only %d sample(s) shared between the matrix columns and coldata rows. Check that sample IDs match exactly.", length(shared)), call. = FALSE)
  }
  mat <- mat[, shared, drop = FALSE]
  coldata <- coldata[shared, , drop = FALSE]

  treat <- as.character(params$treat)
  control <- as.character(params$control)
  levels_present <- unique(as.character(coldata[[group_column]]))
  if (is.null(params$treat) || is.null(params$control)) {
    stop(sprintf("treat and control are required. Levels present: %s", paste(levels_present, collapse = ", ")), call. = FALSE)
  }
  for (lvl in c(treat, control)) {
    if (!lvl %in% levels_present) {
      stop(sprintf("Group '%s' not found. Levels present: %s", lvl, paste(levels_present, collapse = ", ")), call. = FALSE)
    }
  }

  keep <- as.character(coldata[[group_column]]) %in% c(treat, control)
  mat <- mat[, keep, drop = FALSE]
  coldata <- coldata[keep, , drop = FALSE]

  covariates <- params$covariates
  if (!is.null(covariates)) {
    covariates <- as.character(covariates)
    missing <- setdiff(covariates, colnames(coldata))
    if (length(missing)) {
      stop(sprintf("covariates not found in coldata: %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
  }

  n_treat <- sum(as.character(coldata[[group_column]]) == treat)
  n_control <- sum(as.character(coldata[[group_column]]) == control)
  if (n_treat < 2L || n_control < 2L) {
    stop(sprintf("Each group needs at least 2 samples (%s=%d, %s=%d). Differential expression on a single replicate is not interpretable.",
                 treat, n_treat, control, n_control), call. = FALSE)
  }

  safe_coldata <- data.frame(row.names = rownames(coldata))
  safe_coldata$.group <- factor(as.character(coldata[[group_column]]), levels = c(control, treat))
  safe_covariates <- character()
  if (length(covariates)) {
    safe_covariates <- sprintf(".covariate_%d", seq_along(covariates))
    for (index in seq_along(covariates)) {
      value <- coldata[[covariates[[index]]]]
      if (is.character(value)) value <- factor(value)
      safe_coldata[[safe_covariates[[index]]]] <- value
    }
  }

  list(mat = mat, coldata = safe_coldata, group_column = ".group",
       original_group_column = group_column, treat = treat, control = control,
       covariates = safe_covariates, original_covariates = covariates,
       matrix_type = matrix_type,
       n_treat = n_treat, n_control = n_control)
}

run_deseq2 <- function(io, padj_method, alpha) {
  omics_require(c("DESeq2"))
  counts <- round(io$mat)
  storage.mode(counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = io$coldata,
    design = build_design(io$group_column, io$covariates)
  )
  # Standard low-count prefilter; keeps dispersion estimation stable.
  dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(
    dds,
    contrast = c(io$group_column, io$treat, io$control),
    pAdjustMethod = padj_method,
    alpha = alpha
  )
  df <- as.data.frame(res)
  data.frame(
    gene = rownames(df),
    log2FoldChange = df$log2FoldChange,
    baseMean = df$baseMean,
    stat = df$stat,
    pvalue = df$pvalue,
    padj = df$padj,
    stringsAsFactors = FALSE
  )
}

run_edger <- function(io, padj_method) {
  omics_require(c("edgeR"))
  counts <- round(io$mat)
  design <- stats::model.matrix(build_design(io$group_column, io$covariates), data = io$coldata)
  dge <- edgeR::DGEList(counts = counts, group = io$coldata[[io$group_column]])
  keep <- edgeR::filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  # Last coefficient is the treat-vs-control contrast given the control-first factor.
  test <- edgeR::glmQLFTest(fit, coef = ncol(design))
  tab <- edgeR::topTags(test, n = Inf, adjust.method = padj_method, sort.by = "none")$table
  data.frame(
    gene = rownames(tab),
    log2FoldChange = tab$logFC,
    baseMean = 2^tab$logCPM,
    stat = tab$F,
    pvalue = tab$PValue,
    padj = tab$FDR,
    stringsAsFactors = FALSE
  )
}

run_limma <- function(io, padj_method, use_voom) {
  omics_require(c("limma"))
  design <- stats::model.matrix(build_design(io$group_column, io$covariates), data = io$coldata)
  if (isTRUE(use_voom)) {
    omics_require(c("edgeR"))
    dge <- edgeR::DGEList(counts = round(io$mat))
    keep <- edgeR::filterByExpr(dge, design = design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)
    object <- limma::voom(dge, design)
  } else {
    object <- io$mat
  }
  fit <- limma::eBayes(limma::lmFit(object, design))
  tab <- limma::topTable(fit, coef = ncol(design), number = Inf, adjust.method = padj_method, sort.by = "none")
  data.frame(
    gene = rownames(tab),
    log2FoldChange = tab$logFC,
    baseMean = tab$AveExpr,
    stat = tab$t,
    pvalue = tab$P.Value,
    padj = tab$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

handler <- function(params) {
  method <- tolower(as.character(params$method %||% "deseq2")[1])
  padj_method <- as.character(params$padj_method %||% "BH")[1]
  io <- align_inputs(params)

  use_voom <- isTRUE(params$voom)
  if (method %in% c("deseq2", "edger") && io$matrix_type != "counts") {
    stop(sprintf("%s requires matrix_type='counts'.", method), call. = FALSE)
  }
  if (method == "limma" && use_voom && io$matrix_type != "counts") {
    stop("limma-voom requires matrix_type='counts'.", call. = FALSE)
  }
  if (method == "limma" && !use_voom && io$matrix_type != "normalized") {
    stop("limma without voom requires matrix_type='normalized'.", call. = FALSE)
  }
  if (method %in% c("deseq2", "edger")) {
    omics_assert_counts(io$mat, toupper(method))
  } else if (method == "limma" && use_voom) {
    omics_assert_counts(io$mat, "limma-voom")
  }

  lfc_cut <- as.numeric(params$log2fc %||% 1)
  padj_cut <- as.numeric(params$padj %||% 0.05)
  if (!is.finite(padj_cut) || padj_cut <= 0 || padj_cut >= 1 ||
      !is.finite(lfc_cut) || lfc_cut < 0) {
    stop("padj must be between 0 and 1, and log2fc must be non-negative.", call. = FALSE)
  }
  table <- switch(
    method,
    deseq2 = run_deseq2(io, padj_method, padj_cut),
    edger = run_edger(io, padj_method),
    limma = run_limma(io, padj_method, use_voom),
    stop(sprintf("Unknown method '%s'. Use deseq2, edger, or limma.", method), call. = FALSE)
  )

  available <- is.finite(table$padj) & is.finite(table$log2FoldChange)
  table$direction <- ifelse(available, "ns", "unavailable")
  table$direction[available & table$padj < padj_cut & table$log2FoldChange > 0 & table$log2FoldChange >= lfc_cut] <- "up"
  table$direction[available & table$padj < padj_cut & table$log2FoldChange < 0 & table$log2FoldChange <= -lfc_cut] <- "down"
  table <- table[order(table$padj, -abs(table$log2FoldChange)), ]

  name <- as.character(params$output_name %||% sprintf("deg_%s_%s_vs_%s", method, io$treat, io$control))[1]
  all_path <- omics_save_table(table, "deg", name)
  sig <- table[table$direction %in% c("up", "down"), ]
  sig_path <- omics_save_table(sig, "deg", paste0(name, "_significant"))

  list(
    method = method,
    matrix_type = io$matrix_type,
    comparison = sprintf("%s vs %s", io$treat, io$control),
    design = paste(c(io$original_covariates, io$original_group_column), collapse = " + "),
    samples = list(treat = io$n_treat, control = io$n_control),
    covariates = omics_arr(io$original_covariates),
    thresholds = list(log2fc = lfc_cut, padj = padj_cut, padj_method = padj_method),
    genes_tested = nrow(table),
    filtering = list(
      input_genes = nrow(io$mat),
      prefilter_removed = nrow(io$mat) - nrow(table),
      pvalue_unavailable = sum(!is.finite(table$pvalue)),
      padj_unavailable = sum(!is.finite(table$padj)),
      independent_filtering_alpha = if (method == "deseq2") padj_cut else NULL
    ),
    significant = list(
      total = nrow(sig),
      up = sum(table$direction == "up"),
      down = sum(table$direction == "down")
    ),
    top_genes = utils::head(sig[, c("gene", "log2FoldChange", "pvalue", "padj", "direction")], 25),
    result_table = all_path,
    significant_table = sig_path,
    notes = omics_arr(c(
      "log2FoldChange is treat vs control; positive means higher in treat.",
      sprintf("Multiple testing correction: %s.", padj_method),
      "Missing statistical values remain missing (empty CSV cells); unavailable is distinct from non-significant.",
      if (method == "deseq2") "Prefilter: total count >= 10; independent filtering inside results() is separate and may leave padj unavailable."
    ))
  )
}

omics_main(handler)
