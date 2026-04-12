#' Compute adjacency matrix
#'
#' @param x Numeric matrix-like object with samples in rows and genes in columns.
#' @param mode `paper` for `r^2`, `python_compat` for `|r|`.
#' @param method Correlation method.
#' @return Gene-by-gene adjacency matrix.
#' @export
compute_adjacency <- function(
  x,
  mode = c("paper", "python_compat"),
  method = c("pearson", "spearman")
) {
  mode <- match.arg(mode)
  method <- match.arg(method)
  x <- clear_data(x, min_samples = 2L)
  cor_mat <- suppressWarnings(
    stats::cor(x, method = method, use = "pairwise.complete.obs")
  )
  cor_mat[!is.finite(cor_mat)] <- 0
  if (mode == "paper") {
    adj <- cor_mat^2
  } else {
    adj <- abs(cor_mat)
  }
  adj[adj < 0] <- 0
  adj[adj > 1] <- 1
  diag(adj) <- 1
  dimnames(adj) <- list(colnames(x), colnames(x))
  adj
}
