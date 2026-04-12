#' Coerce input to expression matrix
#'
#' @param x Numeric matrix-like object with samples in rows and genes in columns.
#' @return Numeric matrix.
#' @export
as_expression_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("`x` must be a matrix or data.frame.")
  }
  storage.mode(x) <- "double"
  if (nrow(x) < 2L || ncol(x) < 2L) {
    stop("`x` must have at least 2 rows and 2 columns.")
  }
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("sample_", seq_len(nrow(x)))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("gene_", seq_len(ncol(x)))
  }
  x
}

#' Clean expression matrix
#'
#' Removes genes with non-finite values or zero variance.
#'
#' @param x Numeric matrix-like object.
#' @param min_samples Minimum sample count after cleaning.
#' @return Cleaned numeric matrix.
#' @export
clear_data <- function(x, min_samples = 10L) {
  x <- as_expression_matrix(x)
  finite_gene <- apply(x, 2L, function(v) all(is.finite(v)))
  x <- x[, finite_gene, drop = FALSE]
  if (ncol(x) < 2L) {
    stop("Too few genes after removing non-finite genes.")
  }
  vars <- apply(x, 2L, stats::var)
  keep <- is.finite(vars) & vars > 0
  x <- x[, keep, drop = FALSE]
  if (ncol(x) < 2L) {
    stop("Too few genes after removing zero-variance genes.")
  }
  if (nrow(x) < as.integer(min_samples)) {
    stop("Too few samples for analysis.")
  }
  x
}

#' Select top variable genes
#'
#' @param x Numeric matrix-like object.
#' @param n Number of genes to keep.
#' @return Matrix with selected genes.
#' @export
top_variable_genes <- function(x, n = 2000L) {
  x <- clear_data(x, min_samples = 2L)
  n <- min(as.integer(n), ncol(x))
  vars <- apply(x, 2L, stats::var)
  keep <- order(vars, decreasing = TRUE)[seq_len(n)]
  x[, keep, drop = FALSE]
}
