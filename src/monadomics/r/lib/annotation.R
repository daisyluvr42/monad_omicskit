omics_map_ids <- function(ids, orgdb, id_type) {
  omics_require(c("AnnotationDbi", orgdb))
  if (!id_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID")) {
    stop("id_type must be SYMBOL, ENSEMBL, or ENTREZID.", call. = FALSE)
  }
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  db <- getExportedValue(orgdb, orgdb)
  valid <- intersect(ids, AnnotationDbi::keys(db, keytype = id_type))
  mapped <- data.frame(INPUT = character(), ENTREZID = character(), SYMBOL = character())
  if (length(valid)) {
    selected <- suppressMessages(AnnotationDbi::select(
      db, keys = valid, keytype = id_type, columns = c("ENTREZID", "SYMBOL")
    ))
    mapped <- unique(data.frame(INPUT = selected[[id_type]], ENTREZID = selected$ENTREZID,
                                SYMBOL = selected$SYMBOL, stringsAsFactors = FALSE))
    mapped <- mapped[!is.na(mapped$ENTREZID), , drop = FALSE]
  }
  mapping <- merge(data.frame(INPUT = ids), mapped, by = "INPUT", all.x = TRUE, sort = FALSE)
  mapping <- mapping[order(match(mapping$INPUT, ids)), , drop = FALSE]
  ambiguous <- names(which(table(mapped$INPUT) > 1L))
  list(mapped = mapped, table = mapping, entrez = unique(mapped$ENTREZID), input_n = length(ids),
       mapped_n = length(unique(mapped$INPUT)), unmapped = setdiff(ids, mapped$INPUT),
       ambiguous = ambiguous, orgdb = orgdb, version = as.character(utils::packageVersion(orgdb)),
       id_type = id_type)
}

omics_save_annotation <- function(mapping, species, subdir, name) {
  list(species = species, id_type = mapping$id_type, database = mapping$orgdb,
       database_version = mapping$version, input_ids = mapping$input_n,
       mapped_ids = mapping$mapped_n, unmapped_ids = omics_arr(mapping$unmapped),
       ambiguous_ids = omics_arr(mapping$ambiguous),
       mapping_table = omics_save_table(mapping$table, subdir, paste0(name, "_annotation")))
}
