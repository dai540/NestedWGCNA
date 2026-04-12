#' Compute dissimilarity matrix
#'
#' @param adjacency Square adjacency matrix.
#' @return Gene-by-gene dissimilarity matrix.
#' @export
compute_dissimilarity <- function(adjacency) {
  if (!is.matrix(adjacency) || nrow(adjacency) != ncol(adjacency)) {
    stop("`adjacency` must be a square matrix.")
  }
  adj <- adjacency
  adj[!is.finite(adj)] <- 0
  adj[adj < 0] <- 0
  adj[adj > 1] <- 1
  d <- 1 - adj
  d[d < 0] <- 0
  d <- sqrt(d)
  d <- (d + t(d)) / 2
  diag(d) <- 0
  d
}
