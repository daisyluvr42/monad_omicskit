#!/usr/bin/env Rscript
# Install the R packages Scholar's bioinformatics tools depend on.
#
#   Rscript r/bootstrap.R              # core + deg + enrich + plot + survival
#   Rscript r/bootstrap.R enrich       # one group
#   Rscript r/bootstrap.R DESeq2 limma # explicit packages
#   Rscript r/bootstrap.R --check      # report status, install nothing
#
# Kept out of the Scholar installer on purpose: a first Bioconductor build can
# take 20+ minutes and must not be able to fail the MCP install.

GROUPS <- list(
  core = c("jsonlite", "ggplot2"),
  deg = c("DESeq2", "edgeR", "limma"),
  enrich = c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db", "ReactomePA", "GSVA", "msigdbr"),
  # svglite backs ggsave's SVG device; without it every figure call fails at save time.
  plot = c("pheatmap", "ggrepel", "ggvenn", "RColorBrewer", "svglite"),
  survival = c("survival", "survminer", "glmnet", "timeROC", "rms")
)

DEFAULT_GROUPS <- c("core", "deg", "enrich", "plot", "survival")

# Bioconductor packages need BiocManager; everything else comes from CRAN.
BIOC <- c(
  "DESeq2", "edgeR", "limma", "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", "org.Mm.eg.db", "ReactomePA", "GSVA"
)

resolve <- function(argv) {
  argv <- argv[argv != "--check"]
  if (!length(argv)) return(unique(unlist(GROUPS[DEFAULT_GROUPS])))
  out <- character()
  for (item in argv) {
    out <- c(out, if (item %in% names(GROUPS)) GROUPS[[item]] else item)
  }
  unique(out)
}

status <- function(pkgs) {
  vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
}

report <- function(pkgs) {
  have <- status(pkgs)
  for (pkg in pkgs) {
    cat(sprintf("%-20s %s\n", pkg, if (have[[pkg]]) "OK" else "MISSING"))
  }
  cat(sprintf("\n%d/%d installed\n", sum(have), length(have)))
  invisible(have)
}

# A crashed install leaves 00LOCK-<pkg> behind, and every later attempt at that
# package fails until it is removed.
clear_stale_locks <- function() {
  cleared <- character()
  for (lib in .libPaths()) {
    locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
    for (lock in locks) {
      if (unlink(lock, recursive = TRUE, force = TRUE) == 0) {
        cleared <- c(cleared, basename(lock))
      }
    }
  }
  if (length(cleared)) {
    cat(sprintf("Cleared %d stale install lock(s): %s\n", length(cleared), paste(cleared, collapse = ", ")))
  }
  invisible(cleared)
}

warn_if_source_only <- function() {
  if (!identical(getOption("pkgType"), "source")) return(invisible(NULL))
  cat("\n! This R build installs every package from source.\n")
  cat("  Platform:", R.version$platform, "\n")
  cat("  CRAN and Bioconductor ship macOS binaries for the official CRAN build only,\n")
  cat("  so a Homebrew R compiles the whole dependency tree and needs system libraries\n")
  cat("  (cmake, imagemagick, hdf5, ...). Installing the CRAN build of R is far more\n")
  cat("  reliable for Bioconductor: https://cran.r-project.org/bin/macosx/\n\n")
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  pkgs <- resolve(argv)

  if ("--check" %in% argv) {
    report(pkgs)
    return(invisible(NULL))
  }

  if (getRversion() < "4.2.0") {
    stop("R 4.2 or newer is required.", call. = FALSE)
  }

  options(repos = c(CRAN = "https://cloud.r-project.org"))
  # Annotation packages run to hundreds of MB; the 60s default aborts them
  # mid-download on slower or unstable connections.
  options(timeout = max(getOption("timeout"), 1800))
  warn_if_source_only()
  clear_stale_locks()
  missing <- pkgs[!status(pkgs)]
  if (!length(missing)) {
    cat("All requested R packages are already installed.\n")
    return(invisible(NULL))
  }

  cat(sprintf("Installing %d package(s): %s\n", length(missing), paste(missing, collapse = ", ")))

  bioc_missing <- intersect(missing, BIOC)
  cran_missing <- setdiff(missing, BIOC)

  if (length(cran_missing)) {
    utils::install.packages(cran_missing, Ncpus = max(1L, parallel::detectCores() - 1L))
  }
  if (length(bioc_missing)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
  }

  cat("\nFinal status:\n")
  have <- report(pkgs)
  if (!all(have)) {
    cat("\nSome packages failed to install. Read the log above for the first error;\n")
    cat("system libraries (e.g. gfortran, libxml2) are the usual cause on macOS.\n")
    quit(status = 1, save = "no")
  }
}

main()
