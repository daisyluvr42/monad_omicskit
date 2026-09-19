root <- normalizePath(commandArgs(trailingOnly = TRUE)[1])
Sys.setenv(OMICS_R_DIR = file.path(root, "src/monadomics/r"),
           OMICS_OUTPUT_DIR = file.path(root, "test-output/statistical-regression"))
load_functions <- function(name) {
  expressions <- parse(file.path(Sys.getenv("OMICS_R_DIR"), name))
  env <- new.env(parent = globalenv())
  for (expression in expressions[-length(expressions)]) eval(expression, env)
  env
}
deg <- load_functions("deg.R")
plots <- load_functions("plots.R")
enrich <- load_functions("enrich.R")

# Two non-default cutoffs on heterogeneous synthetic counts, including low-count genes.
set.seed(731)
ng <- 2400L
means <- c(rep(1.5, 1600), runif(800, 40, 400))
cts <- matrix(rnbinom(ng * 8, mu = rep(means, 8), size = 5), ng, 8)
cts[1601:1800, 5:8] <- cts[1601:1800, 5:8] * 5L
cts[1801:2000, 1:4] <- cts[1801:2000, 1:4] * 5L
dimnames(cts) <- list(paste0("g", seq_len(ng)), paste0("s", 1:8))
meta <- data.frame(sample = colnames(cts), group = rep(c("control", "treated"), each = 4),
                   donor = rep(paste0("donor", 1:4), 2))
params <- list(method = "deseq2", matrix_type = "counts",
               matrix = data.frame(gene = rownames(cts), cts), coldata = meta,
               group_column = "group", treat = "treated", control = "control", covariates = "donor")
direct_meta <- meta
rownames(direct_meta) <- meta$sample
direct_meta$group <- factor(meta$group, levels = c("control", "treated"))
dds <- DESeq2::DESeqDataSetFromMatrix(cts, direct_meta, ~ donor + group)
dds <- DESeq2::DESeq(dds[rowSums(cts) >= 10, ], quiet = TRUE)
for (alpha in c(0.01, 0.1)) {
  result <- deg$handler(c(params, list(padj = alpha, output_name = paste0("deg_", alpha))))
  actual <- read.csv(result$result_table, row.names = 1)
  expected <- as.data.frame(DESeq2::results(dds, contrast = c("group", "treated", "control"), alpha = alpha))
  actual <- actual[rownames(expected), ]
  stopifnot(identical(is.na(actual$padj), is.na(expected$padj)),
            isTRUE(all.equal(actual$padj, expected$padj, tolerance = 1e-8)),
            isTRUE(all.equal(actual$log2FoldChange, expected$log2FoldChange, tolerance = 1e-8)),
            result$filtering$independent_filtering_alpha == alpha,
            result$filtering$prefilter_removed == sum(rowSums(cts) < 10))
  expected_sig <- is.finite(expected$padj) & expected$padj < alpha & abs(expected$log2FoldChange) >= 1
  significant <- read.csv(result$significant_table)
  stopifnot(setequal(significant$gene, rownames(expected)[expected_sig]),
            nrow(significant) == result$significant$total)
}
cat("PASS direct DESeq2 agreement at two cutoffs with paired synthetic data\n")

# Exercise missing values and zero effect without depending on a particular fitted dataset.
deg$run_deseq2 <- function(io, padj_method, alpha) data.frame(
  gene = c("up", "down", "low", "outlier", "zero"), log2FoldChange = c(2, -2, 3, NA, 0),
  baseMean = 10, stat = c(3, -3, 1, NA, 0), pvalue = c(.001, .002, .1, NA, .001),
  padj = c(.01, .02, NA, NA, .01))
result <- deg$handler(c(params, list(log2fc = 0, output_name = "missing")))
tab <- read.csv(result$result_table)
stopifnot(sum(is.na(tab$padj)) == 2L, sum(is.na(tab$pvalue)) == 1L,
          all(tab$direction[is.na(tab$padj)] == "unavailable"),
          result$significant$total == 2L, tab$direction[tab$gene == "zero"] == "ns")
volcano <- plots$plot_volcano(list(deg_path = result$result_table, log2fc = 0, output_name = "missing_volcano"))
stopifnot(volcano$genes_omitted_unavailable == 2L, volcano$up == 1L, volcano$down == 1L)
cat("PASS missing values, zero effect and volcano counts\n")

prepared <- plots$expression_for_plot(params)
reference <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts = cts)), log = TRUE, prior.count = 2)
stopifnot(isTRUE(all.equal(prepared$matrix, reference)))
heatmap <- plots$plot_heatmap(c(params, list(genes = c("g1601", "g1801", "g2100"), output_name = "normalized_heatmap")))
stopifnot(identical(heatmap$normalization, prepared$normalization), heatmap$genes == 3L)
normalized <- list(matrix_type = "normalized", matrix = data.frame(gene = rownames(reference), reference))
stopifnot(identical(plots$expression_for_plot(normalized)$matrix, reference))
pca <- plots$plot_pca(c(params, list(top_variable = 100, output_name = "normalized_pca")))
keep <- order(apply(reference, 1, var), decreasing = TRUE)[1:100]
expected_pc <- prcomp(t(reference[keep, ]), center = TRUE, scale. = FALSE)
stopifnot(isTRUE(all.equal(abs(as.matrix(pca$scores[, c("PC1", "PC2")])),
                          abs(expected_pc$x[, 1:2]), check.attributes = FALSE)))
scaled <- cts[, 1:4]
scaled[, 2] <- 10L * scaled[, 1]
same <- plots$expression_for_plot(list(matrix_type = "counts", matrix = data.frame(gene = rownames(scaled), scaled)))
raw_depth_difference <- median(abs(log2(scaled[, 1] + 1) - log2(scaled[, 2] + 1)))
stopifnot(max(abs(same$matrix[, 1] - same$matrix[, 2])) < raw_depth_difference / 10)
cat("PASS TMM reference, full-matrix heatmap normalization, normalized passthrough and PCA\n")

# A controlled upstream response tests filtering independently of GO membership and stochastic P values.
namespace <- asNamespace("clusterProfiler")
original_gsea <- get("gseGO", namespace)
original_table <- enrich$enrich_table
fake <- data.frame(ID = c("positive", "negative", "nominal", "q_only", "missing_q", NA),
                   Description = c("positive", "negative", "nominal", "q_only", "missing_q", NA),
                   NES = c(2, -2, 1, 1, 1, NA), pvalue = c(.001, .002, .01, .003, .004, NA),
                   p.adjust = c(.01, .02, .4, .03, .04, NA), qvalue = c(.01, .02, .3, .4, NA, NA),
                   setSize = c(20, 30, 40, 20, 20, NA))
assignInNamespace("gseGO", function(...) NULL, ns = "clusterProfiler")
enrich$enrich_table <- function(result) if (is.null(fake)) data.frame() else fake
gsea_params <- list(ranked = data.frame(gene = as.character(1:100), log2FoldChange = seq(-3, 3, length.out = 100)),
                    id_type = "ENTREZID", pvalue = .05, padj = .05, qvalue = .2, seed = 123,
                    output_name = "gsea_filters")
org <- list(species = "human", orgdb = "org.Hs.eg.db")
result <- enrich$run_gsea(gsea_params, org)
sig <- read.csv(result$significant_table)
stopifnot(result$terms_tested == 5L, result$terms_invalid == 1L, result$terms_significant == 2L,
          setequal(sig$ID, c("positive", "negative")),
          nrow(result$top_terms) == 2L, nrow(read.csv(result$result_table)) == 5L,
          identical(enrich$dot_plot(sig, "test", 10)$data$.x, sig$NES))
gsea_params$qvalue <- NULL
result <- enrich$run_gsea(gsea_params, org)
stopifnot(result$terms_significant == 4L)
old_figure <- result$figure
gsea_params$padj <- .001
result <- enrich$run_gsea(gsea_params, org)
stopifnot(result$terms_significant == 0L, is.null(result$figure), nrow(read.csv(result$significant_table)) == 0L,
          !file.exists(old_figure$png), !file.exists(old_figure$svg))
fake <- NULL
result <- enrich$run_gsea(gsea_params, org)
stopifnot(result$terms_tested == 0L, result$terms_significant == 0L,
          nrow(read.csv(result$result_table)) == 0L,
          nrow(read.csv(result$significant_table)) == 0L)
assignInNamespace("gseGO", original_gsea, ns = "clusterProfiler")
enrich$enrich_table <- original_table
cat("PASS GSEA nominal-only, q-value, missing and empty result contracts\n")

# Real GO execution, species-independent seed contract, no expected pathway count.
keys <- AnnotationDbi::mappedkeys(org.Mm.eg.db::org.Mm.egGO)
set.seed(814)
gsea_params <- list(id_type = "ENTREZID", ranked = data.frame(gene = sample(keys, 700),
                    log2FoldChange = rnorm(700)), seed = 19, pvalue = 1, padj = 1,
                    output_name = "mouse_seed")
mouse <- list(species = "mouse", orgdb = "org.Mm.eg.db")
first <- enrich$run_gsea(gsea_params, mouse)
stopifnot(first$terms_tested > 0L)
first_table <- read.csv(first$result_table)
second <- enrich$run_gsea(gsea_params, mouse)
stopifnot(identical(first_table, read.csv(second$result_table)), first$seed == 19L)
cat("PASS real mouse GO GSEA seed reproducibility\n")
