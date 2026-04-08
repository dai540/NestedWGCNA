ensure_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

ensure_bioc("ALL")
ensure_bioc("Biobase")

suppressPackageStartupMessages({
  library(ALL)
  library(Biobase)
})

data(ALL)

expr <- t(exprs(ALL))
expr <- as.data.frame(expr, check.names = FALSE)

meta <- pData(ALL)[, c("BT", "mol.biol", "sex", "age"), drop = FALSE]
meta <- as.data.frame(meta, check.names = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[1]))
root <- normalizePath(file.path(dirname(script_file), "..", ".."))
script <- file.path(root, "inst", "scripts", "run_matrix_case_study.py")
outdir <- file.path(root, "inst", "extdata", "case_studies", "all_leukemia")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tmp_expr <- tempfile(fileext = ".tsv")
tmp_meta <- tempfile(fileext = ".tsv")
write.table(expr, file = tmp_expr, sep = "\t", quote = FALSE, col.names = NA)
write.table(meta, file = tmp_meta, sep = "\t", quote = FALSE, col.names = NA)

source_note <- sprintf(
  "Bioconductor ALL package v%s (dataset: ALL ExpressionSet)",
  as.character(packageVersion("ALL"))
)

args <- c(
  shQuote(script),
  "--input", shQuote(tmp_expr),
  "--output-dir", shQuote(outdir),
  "--dataset-id", "all_leukemia",
  "--source", shQuote(source_note),
  "--metadata", shQuote(tmp_meta),
  "--cgm-min-size", "30",
  "--fgm-min-size", "8",
  "--top-var-genes", "1200"
)

status <- system2("python", args = args)
if (status != 0) {
  stop("Python case-study pipeline failed for ALL dataset.")
}
