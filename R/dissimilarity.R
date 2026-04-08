symmetrize_matrix <- function(x) {
  (x + t(x)) / 2
}

tom_numerator <- function(a) {
  a %*% a + a
}

tom_denominator <- function(a) {
  s <- colSums(a)
  mm <- outer(s, s, pmin)
  mm + 1 - a
}

mot_denominator <- function(a) {
  s <- colSums(a)
  mm <- outer(s, s, pmax)
  mm + 1 - a
}

#' Compute dissimilarity matrix
#'
#' @param adjacency Adjacency matrix.
#' @param method One of `"paper_metric"`, `"python_compat"`, `"tom"`, `"mot"`,
#'   or `"custom"`.
#' @param custom_fun Optional function for custom dissimilarity.
#'
#' @return Gene-gene dissimilarity matrix.
#' @export
compute_dissimilarity <- function(
    adjacency,
    method = c("paper_metric", "python_compat", "tom", "mot", "custom"),
    custom_fun = NULL) {
  method <- match.arg(method)
  a <- as.matrix(adjacency)
  if (!is.numeric(a) || nrow(a) != ncol(a)) {
    stop("`adjacency` must be a numeric square matrix.", call. = FALSE)
  }

  if (method %in% c("paper_metric", "python_compat")) {
    d <- symmetrize_matrix(1 - a)
    d <- round(d, 12)
    d[d < 0] <- 0
    d <- sqrt(d)
    d <- as.matrix(d)
    diag(d) <- 0
    return(d)
  }

  if (method == "tom") {
    b <- a
    diag(b) <- 0
    d <- 1 - (tom_numerator(b) / tom_denominator(b))
    diag(d) <- 0
    return(d)
  }

  if (method == "mot") {
    b <- a
    diag(b) <- 0
    d <- 1 - (tom_numerator(b) / mot_denominator(b))
    diag(d) <- 0
    return(d)
  }

  if (!is.function(custom_fun)) {
    stop("`custom_fun` must be a function when method = 'custom'.", call. = FALSE)
  }
  d <- custom_fun(a)
  d <- as.matrix(d)
  if (!is.numeric(d) || any(dim(d) != dim(a))) {
    stop("`custom_fun` must return a numeric matrix with same dimension.", call. = FALSE)
  }
  d
}
