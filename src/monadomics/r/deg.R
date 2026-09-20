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
source(file.path(OMICS_R_DIR, "lib", "annotation.R"))

build_design <- function(group_column, covariates) {
  stats::reformulate(c(covariates, group_column))
}

align_inputs <- function(params) {
  matrix_type <- omics_matrix_type(params)
  mat <- omics_read_matrix(params, matrix_type = matrix_type)
  coldata <- omics_read_table(params, "coldata_path", "coldata")
  sample_ids <- as.character(coldata[[1]])
  if (anyNA(sample_ids) || any(!nzchar(trimws(sample_ids))) || anyDuplicated(sample_ids)) {
    stop("Metadata sample IDs must be non-empty and unique.", call. = FALSE)
  }
  rownames(coldata) <- sample_ids

  group_column <- params$group_column
  if (is.null(group_column) || !nzchar(group_column)) {
    stop("group_column is required", call. = FALSE)
  }
  if (!group_column %in% colnames(coldata)) {
    stop(sprintf("group_column '%s' not found in coldata columns: %s",
                 group_column, paste(colnames(coldata), collapse = ", ")), call. = FALSE)
  }

  input_samples <- colnames(mat)
  missing_samples <- setdiff(input_samples, sample_ids)
  if (length(missing_samples)) {
    stop(sprintf("Missing metadata for matrix samples: %s", paste(missing_samples, collapse = ", ")), call. = FALSE)
  }
  coldata <- coldata[input_samples, , drop = FALSE]
  if (anyNA(coldata[[group_column]]) || any(!nzchar(trimws(as.character(coldata[[group_column]]))))) {
    stop("Group labels are missing; resolve sample metadata before analysis.", call. = FALSE)
  }

  treat <- as.character(params$treat)
  control <- as.character(params$control)
  levels_present <- unique(as.character(coldata[[group_column]]))
  if (is.null(params$treat) || is.null(params$control)) {
    stop(sprintf("treat and control are required. Levels present: %s", paste(levels_present, collapse = ", ")), call. = FALSE)
  }
  if (identical(treat, control)) stop("treat and control must be different groups.", call. = FALSE)
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
      if (anyNA(value) || (is.numeric(value) && any(!is.finite(value))) ||
          (is.character(value) && any(!nzchar(trimws(value))))) {
        stop(sprintf("Covariate '%s' contains missing or non-finite values.", covariates[[index]]), call. = FALSE)
      }
      if (is.character(value)) value <- factor(value)
      safe_coldata[[safe_covariates[[index]]]] <- value
    }
  }

  design <- stats::model.matrix(build_design(".group", safe_covariates), data = safe_coldata)
  design_rank <- qr(design)$rank
  if (design_rank < ncol(design)) {
    stop("Design matrix is not full rank; groups and covariates are confounded or redundant.", call. = FALSE)
  }
  if (nrow(design) <= design_rank) {
    stop("Design has no residual degrees of freedom; reduce the model or provide more replicates.", call. = FALSE)
  }

  list(mat = mat, coldata = safe_coldata, group_column = ".group",
       original_group_column = group_column, treat = treat, control = control,
       covariates = safe_covariates, original_covariates = covariates,
       matrix_type = matrix_type, design_matrix = design,
       sample_selection = list(input = omics_arr(input_samples), used = omics_arr(colnames(mat)),
                               excluded_other_groups = omics_arr(setdiff(input_samples, colnames(mat))),
                               metadata_only = omics_arr(setdiff(sample_ids, input_samples))),
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
  design <- io$design_matrix
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
  design <- io$design_matrix
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
  mapping <- NULL
  if (!is.null(params$species) || !is.null(params$id_type)) {
    if (is.null(params$species) || is.null(params$id_type)) {
      stop("Provide both species and id_type for gene annotation, or omit both and report original IDs.", call. = FALSE)
    }
    org <- omics_species(params)
    species <- org$species
    mapping <- omics_map_ids(rownames(io$mat), org$orgdb, params$id_type)
  }

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
  annotation <- NULL
  if (!is.null(mapping)) {
    unambiguous <- mapping$mapped[!mapping$mapped$INPUT %in% mapping$ambiguous, , drop = FALSE]
    table$symbol <- unambiguous$SYMBOL[match(table$gene, unambiguous$INPUT)]
    annotation <- omics_save_annotation(mapping, species, "deg", name)
  }
  all_path <- omics_save_table(table, "deg", name)
  sig <- table[table$direction %in% c("up", "down"), ]
  sig_path <- omics_save_table(sig, "deg", paste0(name, "_significant"))

  list(
    method = method,
    matrix_type = io$matrix_type,
    comparison = sprintf("%s vs %s", io$treat, io$control),
    design = paste(c(io$original_covariates, io$original_group_column), collapse = " + "),
    fit_scope = "two_group_subset",
    group_levels = omics_arr(c(io$control, io$treat)),
    sample_selection = io$sample_selection,
    design_rank = ncol(io$design_matrix),
    residual_df = nrow(io$design_matrix) - ncol(io$design_matrix),
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
    annotation = annotation,
    top_genes = utils::head(sig[, intersect(c("gene", "symbol", "log2FoldChange", "pvalue", "padj", "direction"), names(sig))], 25),
    result_table = all_path,
    significant_table = sig_path,
    notes = omics_arr(c(
      "log2FoldChange is treat vs control; positive means higher in treat.",
      "Only the two requested groups are fitted; this is not a joint fit of every group in the input.",
      if (is.null(annotation)) "Gene symbols were not annotated; report original IDs or obtain a sourced mapping before naming genes.",
      sprintf("Multiple testing correction: %s.", padj_method),
      "Missing statistical values remain missing (empty CSV cells); unavailable is distinct from non-significant.",
      if (method == "deseq2") "Prefilter: total count >= 10; independent filtering inside results() is separate and may leave padj unavailable."
    ))
  )
}

omics_main(handler)
