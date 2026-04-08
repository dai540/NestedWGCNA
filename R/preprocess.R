#' Coerce input to numeric expression matrix
#'
#' @param x Matrix-like object (`matrix` or `data.frame`) with rows as samples
#'   and columns as genes.
#' @param require_names If `TRUE`, requires row and column names.
#'
#' @return Numeric matrix.
#' @export
as_expression_matrix <- function(x, require_names = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("`x` must be a matrix or data.frame.", call. = FALSE)
  }
  storage.mode(x) <- "double"

  if (require_names) {
    if (is.null(rownames(x))) {
      stop("`x` must have rownames (sample ids).", call. = FALSE)
    }
    if (is.null(colnames(x))) {
      stop("`x` must have colnames (gene ids).", call. = FALSE)
    }
  }

  if (nrow(x) < 10) {
    stop("At least 10 samples are required.", call. = FALSE)
  }
  if (ncol(x) < 20) {
    stop("At least 20 genes are required.", call. = FALSE)
  }
  x
}

#' Basic cleaning (Python `clear_data` equivalent)
#'
#' @param x Numeric expression matrix (`samples x genes`).
#' @param fillna If `NULL`, drop genes containing `NA`; otherwise fill `NA`
#'   with this value.
#'
#' @return Cleaned matrix.
#' @export
clear_data <- function(x, fillna = NULL) {
  x <- as_expression_matrix(x, require_names = FALSE)

  if (is.null(fillna)) {
    keep <- colSums(is.na(x)) == 0
    x <- x[, keep, drop = FALSE]
  } else {
    x[is.na(x)] <- fillna
  }

  # Remove constant genes.
  v <- apply(x, 2, stats::var)
  keep_var <- is.finite(v) & v > 0
  x[, keep_var, drop = FALSE]
}

#' Classification error of a vector
#'
#' @param x Vector.
#'
#' @return Classification error in the range `0` to `1`.
#' @export
classification_err <- function(x) {
  tab <- table(x, useNA = "ifany")
  1 - max(tab) / length(x)
}

#' Bootstrap-style filter (Python `bootstrap_filter` equivalent)
#'
#' Keeps genes with classification error larger than `threshold`.
#'
#' @param x Expression matrix (`samples x genes`).
#' @param threshold Threshold (default `1 / exp(1)`).
#'
#' @return Filtered matrix.
#' @export
bootstrap_filter <- function(x, threshold = 1 / exp(1)) {
  x <- as_expression_matrix(x, require_names = FALSE)
  errs <- apply(x, 2, classification_err)
  x[, errs > threshold, drop = FALSE]
}

#' Select top variable genes
#'
#' @param x Expression matrix (`samples x genes`).
#' @param n Number of genes to keep.
#'
#' @return Matrix with top variable genes.
#' @export
top_variable_genes <- function(x, n = 1200) {
  x <- as_expression_matrix(x, require_names = FALSE)
  if (n >= ncol(x)) {
    return(x)
  }
  vars <- apply(x, 2, stats::var)
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(n)]
  x[, keep, drop = FALSE]
}
