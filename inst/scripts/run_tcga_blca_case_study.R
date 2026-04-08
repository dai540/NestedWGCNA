tcga_blca_url <- "https://raw.githubusercontent.com/ilyada/NestedWGCNA/main/data/TCGA_BLCA.tsv"

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[1]))
root <- normalizePath(file.path(dirname(script_file), "..", ".."))

r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
for (f in r_files) {
  sys.source(f, envir = globalenv())
}

tmp <- tempfile(fileext = ".tsv")
utils::download.file(tcga_blca_url, tmp, mode = "wb")
expr <- utils::read.delim(tmp, check.names = FALSE, row.names = 1)
expr <- as.data.frame(expr, check.names = FALSE)

outdir <- file.path(root, "inst", "extdata", "case_studies", "tcga_blca")

run_case_study(
  expr = expr,
  dataset_id = "tcga_blca",
  source = tcga_blca_url,
  output_dir = outdir,
  metadata = NULL,
  mode = "paper",
  cgm_min_cluster_size = 50,
  fgm_min_cluster_size = 10,
  top_n_genes = 1200,
  seed = 42
)
