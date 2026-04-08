ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  if (!bioc) {
    install.packages(pkg)
    return(invisible(TRUE))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[1]))
root <- normalizePath(file.path(dirname(script_file), "..", ".."))

r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
for (f in r_files) {
  sys.source(f, envir = globalenv())
}

ensure_pkg("bladderbatch", bioc = TRUE)
ensure_pkg("Biobase", bioc = TRUE)

suppressPackageStartupMessages({
  library(bladderbatch)
  library(Biobase)
})

data("bladderdata")
expr <- t(Biobase::exprs(bladderEset))
expr <- as.data.frame(expr, check.names = FALSE)
meta <- Biobase::pData(bladderEset)[, c("cancer", "outcome", "batch"), drop = FALSE]
meta <- as.data.frame(meta, check.names = FALSE)

outdir <- file.path(root, "inst", "extdata", "case_studies", "bladderbatch")

source_note <- sprintf(
  "Bioconductor bladderbatch package v%s (dataset: bladderEset)",
  as.character(utils::packageVersion("bladderbatch"))
)

run_case_study(
  expr = expr,
  dataset_id = "bladderbatch",
  source = source_note,
  output_dir = outdir,
  metadata = meta,
  mode = "paper",
  cgm_min_cluster_size = 12,
  fgm_min_cluster_size = 6,
  top_n_genes = 1200,
  seed = 42
)
