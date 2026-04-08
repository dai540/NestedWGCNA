#' Weighted core decomposition
#'
#' @param x Weighted adjacency-like square matrix.
#'
#' @return List with integer `order` and numeric `degeneracy`.
#' @export
core_decompose_weighted <- function(x) {
  a <- as.matrix(x)
  if (!is.numeric(a) || nrow(a) != ncol(a)) {
    stop("`x` must be a numeric square matrix.", call. = FALSE)
  }
  n <- nrow(a)
  conn <- colSums(a)
  idx <- integer(n)
  deg <- numeric(n)

  for (i in seq_len(n)) {
    j <- which.min(conn)
    conn <- conn - a[j, ]
    idx[i] <- j
    deg[i] <- conn[j]
    conn[j] <- Inf
  }
  list(order = idx, degeneracy = deg)
}

#' Erode clusters to core genes
#'
#' @param clusters Integer cluster labels (`-1` for noise).
#' @param dissimilarity Gene-gene dissimilarity matrix.
#'
#' @return Updated core cluster labels.
#' @export
erode_to_core <- function(clusters, dissimilarity) {
  cl <- as.integer(clusters)
  d <- as.matrix(dissimilarity)
  if (nrow(d) != length(cl) || ncol(d) != length(cl)) {
    stop("`dissimilarity` dimension must match `clusters` length.", call. = FALSE)
  }

  out <- cl
  if (all(cl < 0)) {
    return(out)
  }

  for (k in 0:max(cl)) {
    members <- which(cl == k)
    if (length(members) == 0) {
      next
    }
    if (length(members) == 1) {
      out[members] <- k
      next
    }
    sub_d <- d[members, members, drop = FALSE]
    cd <- core_decompose_weighted(1 - sub_d)
    cut_idx <- which.max(cd$degeneracy)
    keep_local <- cd$order[cut_idx:length(cd$order)]
    out[members] <- -1L
    out[members[keep_local]] <- k
  }
  out
}
