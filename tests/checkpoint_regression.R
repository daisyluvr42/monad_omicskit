root <- normalizePath(commandArgs(trailingOnly = TRUE)[1])
Sys.setenv(OMICS_R_DIR = file.path(root, "src/monadomics/r"),
           OMICS_OUTPUT_DIR = file.path(root, "test-output/checkpoint-regression"))
load_functions <- function(name) {
  expressions <- parse(file.path(Sys.getenv("OMICS_R_DIR"), name))
  env <- new.env(parent = globalenv())
  for (expression in expressions[-length(expressions)]) eval(expression, env)
  env
}
deg <- load_functions("deg.R")
enrich <- load_functions("enrich.R")
reject <- function(expression, pattern) {
  error <- tryCatch({force(expression); NULL}, error = conditionMessage)
  stopifnot(!is.null(error), grepl(pattern, error, ignore.case = TRUE))
}

set.seed(218)
db <- org.Hs.eg.db::org.Hs.eg.db
ids <- head(AnnotationDbi::keys(db, keytype = "ENTREZID"), 120)
counts <- matrix(rnbinom(120 * 12, mu = 150, size = 10), 120, 12,
                 dimnames = list(ids, paste0("sample", 1:12)))
metadata <- data.frame(sample = colnames(counts), group = rep(c("left", "right", "unused"), each = 4))
p <- list(method = "edger", matrix_type = "counts", matrix = data.frame(gene = ids, counts),
          coldata = metadata, group_column = "group", control = "left", treat = "right", output_name = "pairwise")
bad <- p
bad$matrix$gene[2] <- bad$matrix$gene[1]
reject(deg$align_inputs(bad), "unique")
bad$matrix$gene[2] <- ""
reject(deg$align_inputs(bad), "non-empty")
bad <- p
bad$matrix[1,2] <- Inf
reject(deg$align_inputs(bad), "non-finite")
bad <- p
names(bad$matrix)[3] <- names(bad$matrix)[2]
reject(deg$align_inputs(bad), "sample IDs")
bad <- p
bad$coldata <- metadata[-1,]
reject(deg$align_inputs(bad), "Missing metadata")
bad <- p
bad$coldata$sample[2] <- bad$coldata$sample[1]
reject(deg$align_inputs(bad), "sample IDs")
bad <- p
bad$coldata$group[1] <- NA
reject(deg$align_inputs(bad), "labels are missing")
bad <- p
bad$treat <- bad$control
reject(deg$align_inputs(bad), "different groups")
bad <- p
bad$covariates <- "batch"
bad$coldata$batch <- rep(1:4, 3)
bad$coldata$batch[1] <- NA
reject(deg$align_inputs(bad), "Covariate.*missing")
bad$coldata$batch <- metadata$group
reject(deg$align_inputs(bad), "full rank")
# Six independent sample effects plus group/intercept saturate the eight selected observations.
bad$covariates <- c("x1", "x2", "x3", "x4", "x5", "x6")
for (j in 1:6) bad$coldata[[paste0("x", j)]] <- as.integer(seq_len(12) == c(1,2,3,5,6,7)[j])
reject(deg$align_inputs(bad), "residual degrees")
cat("PASS matrix IDs, finite values, sample coverage, groups and model estimability\n")

p$species <- "human"
p$id_type <- "ENTREZID"
first <- deg$handler(p)
stopifnot(first$fit_scope == "two_group_subset", identical(as.character(first$sample_selection$used), colnames(counts)[1:8]),
          identical(as.character(first$sample_selection$excluded_other_groups), colnames(counts)[9:12]),
          first$design_rank == 2L, first$residual_df == 6L, first$annotation$database == "org.Hs.eg.db")
tab <- read.csv(first$result_table, colClasses = c(gene = "character"))
reference <- suppressMessages(AnnotationDbi::select(db, keys = tab$gene, keytype = "ENTREZID", columns = "SYMBOL"))
stopifnot(identical(tab$symbol, reference$SYMBOL[match(tab$gene, reference$ENTREZID)]))
shuffled <- p
shuffled$coldata <- metadata[sample(12),]
shuffled$matrix <- p$matrix[,c(1, sample(2:13))]
shuffled$output_name <- "reordered"
second <- deg$handler(shuffled)
tab2 <- read.csv(second$result_table, colClasses = c(gene = "character"))
tab2 <- tab2[match(tab$gene,tab2$gene),]
stopifnot(isTRUE(all.equal(tab$log2FoldChange, tab2$log2FoldChange, tolerance = 1e-8)),
          isTRUE(all.equal(tab$padj, tab2$padj, tolerance = 1e-8)))
p$log2fc <- 1e9
p$output_name <- "no_degs"
empty <- deg$handler(p)
stopifnot(empty$significant$total == 0L, nrow(read.csv(empty$significant_table)) == 0L)
bad <- p
bad$species <- NULL
reject(deg$handler(bad), "both species and id_type")
cat("PASS actual fitting scope, sourced symbols, order independence and zero significant genes\n")

mapped <- omics_map_ids(c(ids[1], "NOT_A_REAL_ENTREZ_ID"), "org.Hs.eg.db", "ENTREZID")
stopifnot(mapped$input_n == 2L, mapped$mapped_n == 1L, identical(mapped$unmapped, "NOT_A_REAL_ENTREZ_ID"))
none <- omics_map_ids("NOT_A_REAL_ENTREZ_ID", "org.Hs.eg.db", "ENTREZID")
stopifnot(none$mapped_n == 0L, nrow(none$table) == 1L, is.na(none$table$SYMBOL))
ens <- head(AnnotationDbi::keys(db, keytype = "ENSEMBL"), 1000)
mapping <- omics_map_ids(ens, "org.Hs.eg.db", "ENSEMBL")
expected <- suppressMessages(AnnotationDbi::select(db, keys = ens, keytype = "ENSEMBL", columns = c("ENTREZID","SYMBOL")))
stopifnot(nrow(mapping$mapped) == nrow(unique(expected[!is.na(expected$ENTREZID),])),
          setequal(mapping$ambiguous, names(which(table(mapping$mapped$INPUT)>1))))
cat("PASS validated ENTREZ IDs, unmapped IDs and complete one-to-many annotation\n")

go_ids <- AnnotationDbi::mappedkeys(org.Hs.eg.db::org.Hs.egGO)
without_go <- head(setdiff(AnnotationDbi::keys(db, keytype = "ENTREZID"), go_ids), 1)
background <- c(head(go_ids, 500), without_go)
ora <- list(genes = c(head(go_ids, 80), without_go), universe = background, id_type = "ENTREZID",
            ontology = "BP", pvalue = 1, qvalue = 1, output_name = "coverage")
org <- list(species = "human", orgdb = "org.Hs.eg.db")
result <- enrich$run_ora(ora, org, "go")
terms <- read.csv(result$result_table)
stopifnot(nrow(terms)>0L, result$directionality == "unsigned_gene_set",
          result$effective_counts$foreground < result$genes_mapped,
          all(as.integer(sub(".*/", "", terms$GeneRatio)) == result$effective_counts$foreground),
          all(as.integer(sub(".*/", "", terms$BgRatio)) == result$effective_counts$background))
old_figure <- result$figure
ora$pvalue <- .Machine$double.xmin
empty <- enrich$run_ora(ora, org, "go")
stopifnot(empty$terms_significant == 0L, nrow(read.csv(empty$result_table)) == 0L,
          is.null(empty$figure), !file.exists(old_figure$png), !file.exists(old_figure$svg),
          identical(empty$effective_counts,result$effective_counts))
ora$genes <- "NOT_A_REAL_ENTREZ_ID"
reject(enrich$run_ora(ora, org, "go"), "No foreground IDs mapped")
ora$genes <- go_ids[501]
reject(enrich$run_ora(ora, org, "go"), "outside.*background")
cat("PASS effective ORA coverage, unsigned interpretation, empty results and background boundary\n")
