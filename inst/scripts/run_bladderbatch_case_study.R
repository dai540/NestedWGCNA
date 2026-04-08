ensure_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

ensure_bioc("bladderbatch")
ensure_bioc("Biobase")

suppressPackageStartupMessages({
  library(bladderbatch)
  library(Biobase)
})

data("bladderdata")

expr <- t(exprs(bladderEset))
expr <- as.data.frame(expr, check.names = FALSE)

meta <- pData(bladderEset)[, c("cancer", "outcome", "batch"), drop = FALSE]
meta <- as.data.frame(meta, check.names = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[1]))
root <- normalizePath(file.path(dirname(script_file), "..", ".."))

script <- file.path(root, "inst", "scripts", "run_matrix_case_study.py")
outdir <- file.path(root, "inst", "extdata", "case_studies", "bladderbatch")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tmp_expr <- tempfile(fileext = ".tsv")
tmp_meta <- tempfile(fileext = ".tsv")
write.table(expr, file = tmp_expr, sep = "\t", quote = FALSE, col.names = NA)
write.table(meta, file = tmp_meta, sep = "\t", quote = FALSE, col.names = NA)

source_note <- sprintf(
  "Bioconductor bladderbatch package v%s (dataset: bladderEset)",
  as.character(packageVersion("bladderbatch"))
)

args <- c(
  shQuote(script),
  "--input", shQuote(tmp_expr),
  "--output-dir", shQuote(outdir),
  "--dataset-id", "bladderbatch",
  "--source", shQuote(source_note),
  "--metadata", shQuote(tmp_meta),
  "--cgm-min-size", "12",
  "--fgm-min-size", "6",
  "--top-var-genes", "1200"
)

status <- system2("python", args = args)
if (status != 0) {
  stop("Python case-study pipeline failed for bladderbatch dataset.")
}
