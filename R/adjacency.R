#' Compute adjacency matrix
#'
#' @param x Expression matrix (`samples x genes`).
#' @param method Correlation method (`"pearson"` or `"spearman"`).
#' @param mode Similarity mode: `"paper"`, `"python_compat"`, or `"custom"`.
#' @param signed Whether to use signed adjacency.
#' @param power Optional power parameter. If `NULL`, defaults by mode.
#'
#' @return Gene-gene adjacency matrix.
#' @export
compute_adjacency <- function(
    x,
    method = c("pearson", "spearman"),
    mode = c("paper", "python_compat", "custom"),
    signed = FALSE,
    power = NULL) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  x <- as_expression_matrix(x, require_names = FALSE)

  cor_mat <- stats::cor(x, method = method, use = "pairwise.complete.obs")
  cor_mat[!is.finite(cor_mat)] <- 0

  if (mode == "paper") {
    # Paper-style default: r^2, with optional signed variant.
    if (is.null(power)) {
      power <- 2
    }
  } else if (mode == "python_compat") {
    # Python get_clust() behavior uses alpha = 1 in get_adjacency().
    if (is.null(power)) {
      power <- 1
    }
  } else if (is.null(power)) {
    power <- 1
  }

  if (signed) {
    adj <- ((1 + cor_mat) / 2) ^ power
  } else {
    adj <- abs(cor_mat) ^ power
  }
  diag(adj) <- 1
  adj
}
