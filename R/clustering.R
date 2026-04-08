#' Run UMAP on genes using precomputed dissimilarity
#'
#' @param dissimilarity Gene-gene dissimilarity matrix.
#' @param n_components UMAP embedding dimensions.
#' @param n_neighbors UMAP neighbors.
#' @param seed Random seed.
#'
#' @return Numeric embedding matrix.
#' @export
run_gene_umap <- function(
    dissimilarity,
    n_components = 10,
    n_neighbors = 30,
    seed = 42) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required.", call. = FALSE)
  }
  d <- as.matrix(dissimilarity)
  if (!is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`dissimilarity` must be a numeric square matrix.", call. = FALSE)
  }
  n_genes <- nrow(d)
  if (n_genes < 3) {
    stop("At least 3 genes are required for UMAP embedding.", call. = FALSE)
  }
  n_components <- max(2L, min(as.integer(n_components), n_genes - 1L))
  n_neighbors <- max(2L, min(as.integer(n_neighbors), n_genes - 1L))

  set.seed(seed)
  umap_args <- list(
    n_components = n_components,
    n_neighbors = n_neighbors,
    min_dist = 0,
    n_threads = 1,
    verbose = FALSE,
    ret_model = FALSE
  )
  if ("input" %in% names(formals(uwot::umap))) {
    emb <- do.call(
      uwot::umap,
      c(list(X = stats::as.dist(d), input = "dist"), umap_args)
    )
    return(emb)
  }

  # Fallback for newer uwot releases where direct distance input is unavailable:
  # project the dissimilarity matrix with classical MDS and run UMAP on that space.
  mds_k <- max(2L, min(50L, n_genes - 1L))
  mds <- stats::cmdscale(stats::as.dist(d), k = mds_k)
  emb <- do.call(
    uwot::umap,
    c(list(X = mds, metric = "euclidean"), umap_args)
  )
  emb
}

#' Run HDBSCAN on gene embedding
#'
#' @param embedding UMAP embedding.
#' @param min_cluster_size Minimum cluster size.
#' @param min_samples Ignored in current implementation; kept for API parity.
#'
#' @return Integer labels (`-1` for noise, `0..k-1` for clusters).
#' @export
run_gene_hdbscan <- function(
    embedding,
    min_cluster_size = 30,
    min_samples = NULL) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required.", call. = FALSE)
  }
  emb <- as.matrix(embedding)
  if (!is.numeric(emb) || nrow(emb) < 2) {
    stop("`embedding` must be a numeric matrix with at least 2 rows.", call. = FALSE)
  }
  min_pts <- if (is.null(min_samples)) {
    max(2L, as.integer(min_cluster_size))
  } else {
    max(2L, as.integer(min_samples))
  }
  fit <- dbscan::hdbscan(emb, minPts = min_pts)
  # dbscan::hdbscan uses 0 as noise and 1..k as cluster labels.
  ifelse(fit$cluster == 0L, -1L, fit$cluster - 1L)
}

#' Discover coarse-grained modules (CGMs)
#'
#' @param x Expression matrix (`samples x genes`).
#' @param mode Similarity mode.
#' @param corr_method Correlation method.
#' @param min_cluster_size Minimum cluster size.
#' @param n_components UMAP dimensions; if `NULL`, computed as in Python code.
#' @param seed Random seed.
#'
#' @return List with CGM assignments, core assignments, embedding, and matrices.
#' @export
find_cgm <- function(
    x,
    mode = c("paper", "python_compat"),
    corr_method = c("pearson", "spearman"),
    min_cluster_size = 30,
    n_components = NULL,
    seed = 42) {
  mode <- match.arg(mode)
  corr_method <- match.arg(corr_method)
  x <- as_expression_matrix(x, require_names = FALSE)

  if (is.null(n_components)) {
    n_components <- as.integer(round(log2(nrow(x)), 0) + 1)
  }
  n_neighbors <- min(2L * as.integer(min_cluster_size), ncol(x) - 1L)

  adj <- compute_adjacency(
    x = x,
    method = corr_method,
    mode = mode,
    signed = FALSE,
    power = if (mode == "paper") 2 else 1
  )
  dis <- compute_dissimilarity(
    adjacency = adj,
    method = if (mode == "paper") "paper_metric" else "python_compat"
  )
  emb <- run_gene_umap(
    dissimilarity = dis,
    n_components = n_components,
    n_neighbors = n_neighbors,
    seed = seed
  )
  labels <- run_gene_hdbscan(
    embedding = emb,
    min_cluster_size = min_cluster_size,
    min_samples = min_cluster_size %/% 2
  )
  core <- erode_to_core(labels, dis)

  list(
    clusters = labels,
    core_clusters = core,
    embedding = emb,
    adjacency = adj,
    dissimilarity = dis
  )
}
