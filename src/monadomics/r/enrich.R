#!/usr/bin/env Rscript
# Functional enrichment: GO / KEGG / Reactome over-representation, GSEA, GSVA, ssGSEA.
#
# params:
#   method       go | kegg | reactome | gsea | gsva | ssgsea
#   species      human | mouse  (required; identifiers differ and must not be guessed)
#   genes        character vector of gene symbols (ORA methods)
#   ranked       named list/data.frame of gene + metric, or gene_column/metric_column
#                with a table (GSEA); ranking metric is usually log2FoldChange
#   gene_column, metric_column, ranked_path / ranked
#   ontology     GO only: BP (default) | CC | MF | ALL
#   universe     optional background gene symbols for ORA
#   pvalue, qvalue  cutoffs (default 0.05, 0.2)
#   gene_sets    gsva/ssgsea: required named list of gene sets
#   matrix_path  gsva/ssgsea: expression matrix
#   top_n        how many terms to plot (default 10)
#   output_name  basename for saved outputs

OMICS_R_DIR <- Sys.getenv("OMICS_R_DIR")
source(file.path(OMICS_R_DIR, "lib", "common.R"))

input_id_type <- function(params) {
  id_type <- toupper(as.character(params$id_type %||% "")[1])
  if (!id_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID")) {
    stop("id_type is required and must be SYMBOL, ENSEMBL, or ENTREZID.", call. = FALSE)
  }
  id_type
}

map_ids <- function(ids, orgdb, id_type) {
  omics_require(c("clusterProfiler", orgdb))
  ids <- unique(stats::na.omit(as.character(ids)))
  ids <- ids[nzchar(ids)]
  if (!length(ids)) stop("No gene identifiers supplied.", call. = FALSE)
  if (id_type == "ENTREZID") {
    mapped <- data.frame(INPUT = ids, ENTREZID = ids, stringsAsFactors = FALSE)
  } else {
    mapped <- suppressWarnings(suppressMessages(
      clusterProfiler::bitr(ids, fromType = id_type, toType = "ENTREZID", OrgDb = orgdb)
    ))
    if (nrow(mapped)) colnames(mapped)[colnames(mapped) == id_type] <- "INPUT"
  }
  if (!nrow(mapped)) {
    stop(sprintf("None of the supplied %s identifiers mapped to ENTREZID. Check species and identifier type.", id_type), call. = FALSE)
  }
  list(
    mapped = mapped,
    entrez = unique(mapped$ENTREZID),
    input_n = length(ids),
    mapped_n = length(unique(mapped$INPUT)),
    unmapped = setdiff(ids, mapped$INPUT)
  )
}

enrich_table <- function(result) {
  if (is.null(result)) return(data.frame())
  df <- as.data.frame(result)
  if (!nrow(df)) return(df)
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID", "NES", "setSize", "enrichmentScore", "core_enrichment"), colnames(df))
  df[, keep, drop = FALSE]
}

dot_plot <- function(df, title, top_n) {
  omics_require(c("ggplot2"))
  df <- utils::head(df[order(df$p.adjust), ], top_n)
  if (!nrow(df)) return(NULL)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  count_col <- if ("Count" %in% colnames(df)) "Count" else "setSize"
  ggplot2::ggplot(df, ggplot2::aes(x = -log10(p.adjust), y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = .data[[count_col]], colour = p.adjust)) +
    ggplot2::scale_colour_gradient(low = "#b5483a", high = "#2f5f8f") +
    ggplot2::labs(x = expression(-log[10](adjusted~p)), y = NULL, title = title,
                  size = "Genes", colour = "p.adjust") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

run_ora <- function(params, org, method) {
  genes <- params$genes
  if (is.null(genes)) stop("genes is required for over-representation analysis.", call. = FALSE)
  id_type <- input_id_type(params)
  mapping <- map_ids(genes, org$orgdb, id_type)

  universe_entrez <- NULL
  if (!is.null(params$universe)) {
    universe_entrez <- map_ids(params$universe, org$orgdb, id_type)$entrez
  }

  pvalue <- as.numeric(params$pvalue %||% 0.05)
  qvalue <- as.numeric(params$qvalue %||% 0.2)

  result <- switch(
    method,
    go = {
      ontology <- toupper(as.character(params$ontology %||% "BP")[1])
      clusterProfiler::enrichGO(
        gene = mapping$entrez, OrgDb = org$orgdb, keyType = "ENTREZID",
        ont = ontology, pvalueCutoff = pvalue, qvalueCutoff = qvalue,
        universe = universe_entrez, readable = TRUE
      )
    },
    kegg = {
      res <- clusterProfiler::enrichKEGG(
        gene = mapping$entrez, organism = org$kegg,
        pvalueCutoff = pvalue, qvalueCutoff = qvalue, universe = universe_entrez
      )
      if (!is.null(res)) res <- clusterProfiler::setReadable(res, OrgDb = org$orgdb, keyType = "ENTREZID")
      res
    },
    reactome = {
      omics_require("ReactomePA")
      res <- ReactomePA::enrichPathway(
        gene = mapping$entrez, organism = org$reactome,
        pvalueCutoff = pvalue, qvalueCutoff = qvalue, universe = universe_entrez, readable = TRUE
      )
      res
    }
  )

  df <- enrich_table(result)
  name <- as.character(params$output_name %||% paste0("enrich_", method))[1]
  table_path <- if (nrow(df)) omics_save_table(df, "enrich", name) else NULL

  figure <- NULL
  if (nrow(df)) {
    top_n <- as.integer(params$top_n %||% 10)
    label <- switch(method, go = paste0("GO ", toupper(as.character(params$ontology %||% "BP")[1])), kegg = "KEGG", reactome = "Reactome")
    plot <- dot_plot(df, paste0(label, " enrichment"), top_n)
    if (!is.null(plot)) figure <- omics_save_plot(plot, "enrich", paste0(name, "_dotplot"), width = 7.5, height = 5.5)
  }

  list(
    method = method,
    species = org$species,
    id_type = id_type,
    ontology = if (method == "go") toupper(as.character(params$ontology %||% "BP")[1]) else NULL,
    genes_supplied = mapping$input_n,
    genes_mapped = mapping$mapped_n,
    genes_unmapped = omics_arr(mapping$unmapped),
    background = if (is.null(universe_entrez)) "genome-wide default" else sprintf("user-supplied (%d genes)", length(universe_entrez)),
    thresholds = list(pvalue = pvalue, qvalue = qvalue, padj_method = "BH"),
    terms_significant = nrow(df),
    top_terms = utils::head(df, 15),
    result_table = table_path,
    figure = figure,
    notes = omics_arr(c(
      "Enrichment shows association between a gene list and annotated terms; it is not evidence of mechanism or causation.",
      "Unmapped symbols were excluded; check them for outdated aliases before reporting gene counts."
    ))
  )
}

run_gsea <- function(params, org) {
  omics_require(c("clusterProfiler", org$orgdb))
  id_type <- input_id_type(params)
  tab <- omics_read_table(params, "ranked_path", "ranked")
  gene_column <- as.character(params$gene_column %||% colnames(tab)[1])[1]
  metric_column <- as.character(params$metric_column %||% "log2FoldChange")[1]
  for (col in c(gene_column, metric_column)) {
    if (!col %in% colnames(tab)) {
      stop(sprintf("Column '%s' not found. Available: %s", col, paste(colnames(tab), collapse = ", ")), call. = FALSE)
    }
  }
  metric <- suppressWarnings(as.numeric(tab[[metric_column]]))
  invalid_metric <- !is.na(tab[[metric_column]]) & is.na(metric)
  if (any(invalid_metric)) {
    stop(sprintf("Ranking column '%s' contains non-numeric values.", metric_column), call. = FALSE)
  }
  tab[[metric_column]] <- metric
  mapping <- map_ids(tab[[gene_column]], org$orgdb, id_type)
  merged <- merge(tab, mapping$mapped, by.x = gene_column, by.y = "INPUT")
  merged <- merged[is.finite(merged[[metric_column]]), , drop = FALSE]
  grouped <- split(merged[[metric_column]], merged$ENTREZID)
  ranks <- vapply(grouped, function(values) values[[which.max(abs(values))]], numeric(1))
  ranks <- sort(ranks, decreasing = TRUE)
  if (length(ranks) < 50) {
    stop(sprintf("Only %d ranked genes after mapping. GSEA needs the full ranked list, not a filtered DEG subset.", length(ranks)), call. = FALSE)
  }

  pvalue <- as.numeric(params$pvalue %||% 0.05)
  result <- clusterProfiler::gseGO(
    geneList = ranks, OrgDb = org$orgdb, keyType = "ENTREZID",
    ont = toupper(as.character(params$ontology %||% "BP")[1]),
    pvalueCutoff = pvalue, verbose = FALSE
  )
  if (!is.null(result)) result <- clusterProfiler::setReadable(result, OrgDb = org$orgdb, keyType = "ENTREZID")
  df <- enrich_table(result)

  name <- as.character(params$output_name %||% "gsea")[1]
  table_path <- if (nrow(df)) omics_save_table(df, "enrich", name) else NULL
  figure <- NULL
  if (nrow(df)) {
    plot <- dot_plot(df, "GSEA", as.integer(params$top_n %||% 10))
    if (!is.null(plot)) figure <- omics_save_plot(plot, "enrich", paste0(name, "_dotplot"), width = 7.5, height = 5.5)
  }

  list(
    method = "gsea",
    species = org$species,
    id_type = id_type,
    genes_ranked = length(ranks),
    ranking_metric = metric_column,
    thresholds = list(pvalue = pvalue, padj_method = "BH"),
    terms_significant = nrow(df),
    top_terms = utils::head(df[order(df$p.adjust), c("ID", "Description", "NES", "pvalue", "p.adjust", "setSize")], 15),
    result_table = table_path,
    figure = figure,
    notes = omics_arr(c(
      "GSEA requires the complete ranked gene list, not a pre-filtered DEG set.",
      "NES sign follows the ranking metric; positive means enriched at the top of the ranking."
    ))
  )
}

run_gsva <- function(params, org, method) {
  omics_require(c("GSVA"))
  id_type <- input_id_type(params)
  matrix_type <- omics_matrix_type(params)
  mat <- omics_read_matrix(params, matrix_type = matrix_type)
  gene_sets <- params$gene_sets
  if (is.null(gene_sets)) {
    stop("gene_sets is required: a named list of gene symbol vectors (e.g. Hallmark sets).", call. = FALSE)
  }
  gene_sets <- lapply(gene_sets, function(x) unique(as.character(x)))
  covered <- vapply(gene_sets, function(g) sum(g %in% rownames(mat)), integer(1))
  required <- vapply(
    gene_sets,
    function(g) as.integer(min(length(g), max(5L, ceiling(0.2 * length(g))))),
    integer(1)
  )
  usable <- covered >= required
  dropped <- names(gene_sets)[!usable]
  gene_sets <- gene_sets[usable]
  covered <- covered[usable]
  if (!length(gene_sets)) {
    stop("No gene-set members are present in the expression matrix. Check that both use the same identifier type.", call. = FALSE)
  }

  kcdf <- if (matrix_type == "counts") "Poisson" else "Gaussian"
  scores <- if (utils::packageVersion("GSVA") >= "1.52.0") {
    par <- if (method == "gsva") {
      GSVA::gsvaParam(exprData = mat, geneSets = gene_sets, kcdf = kcdf)
    } else {
      GSVA::ssgseaParam(exprData = mat, geneSets = gene_sets)
    }
    GSVA::gsva(par, verbose = FALSE)
  } else {
    GSVA::gsva(mat, gene_sets, method = method, kcdf = kcdf, verbose = FALSE)
  }

  out <- data.frame(gene_set = rownames(scores), as.data.frame(scores), check.names = FALSE)
  name <- as.character(params$output_name %||% method)[1]
  table_path <- omics_save_table(out, "enrich", name)

  list(
    method = method,
    species = org$species,
    id_type = id_type,
    matrix_type = matrix_type,
    gene_sets = length(gene_sets),
    genes_per_set_found = as.list(covered),
    gene_sets_dropped_low_coverage = omics_arr(dropped),
    samples = ncol(mat),
    kcdf = if (method == "gsva") kcdf else NULL,
    score_table = table_path,
    notes = omics_arr(c(
      sprintf("%s returns per-sample pathway scores; compare them between groups with a downstream test.", toupper(method)),
      "Scores are relative within this dataset and are not comparable across datasets."
    ))
  )
}

handler <- function(params) {
  method <- tolower(as.character(params$method %||% "go")[1])
  org <- omics_species(params)
  switch(
    method,
    go = ,
    kegg = ,
    reactome = run_ora(params, org, method),
    gsea = run_gsea(params, org),
    gsva = ,
    ssgsea = run_gsva(params, org, method),
    stop(sprintf("Unknown method '%s'. Use go, kegg, reactome, gsea, gsva, or ssgsea.", method), call. = FALSE)
  )
}

omics_main(handler)
