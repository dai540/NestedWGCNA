#' Calculate per-sample module scores
#'
#' @param x Expression matrix (`samples x genes`).
#' @param assignment Integer module assignment vector for genes (`-1` = noise).
#' @param prefix Prefix for score column names.
#'
#' @return Data frame (`samples x modules`) of module scores.
#' @export
module_score <- function(x, assignment, prefix = "M_") {
  x <- as_expression_matrix(x, require_names = FALSE)
  a <- as.integer(assignment)
  if (length(a) != ncol(x)) {
    stop("`assignment` length must match number of genes in `x`.", call. = FALSE)
  }
  out <- list()
  for (k in sort(unique(a[a >= 0]))) {
    genes <- which(a == k)
    if (length(genes) == 0) {
      next
    }
    out[[paste0(prefix, k)]] <- row_means(x[, genes, drop = FALSE])
  }
  if (length(out) == 0) {
    return(data.frame(row.names = rownames(x)))
  }
  as.data.frame(out, check.names = FALSE)
}

#' Discover FGMs inside a CGM
#'
#' @param x Expression matrix (`samples x genes`) for a single target CGM.
#' @param mode Similarity mode.
#' @param min_cluster_size Minimum FGM cluster size.
#' @param genfocus_corr_method Correlation method for GenFocus.
#' @param genfocus_corr_thr Correlation threshold for INGS.
#' @param genfocus_CVR_thr CVR clipping threshold.
#' @param seed Random seed.
#'
#' @return List with FGM assignments and normalization output.
#' @export
find_fgm <- function(
    x,
    mode = c("paper", "python_compat"),
    min_cluster_size = 10,
    genfocus_corr_method = c("spearman", "pearson"),
    genfocus_corr_thr = 0.9,
    genfocus_CVR_thr = 0.6,
    seed = 42) {
  mode <- match.arg(mode)
  genfocus_corr_method <- match.arg(genfocus_corr_method)
  x <- as_expression_matrix(x, require_names = FALSE)

  thresholds <- unique(c(genfocus_corr_thr, 0.85, 0.8, 0.75, 0.7))
  gf <- NULL
  used_thr <- NA_real_
  for (thr in thresholds) {
    gf <- genfocus_normalize(
      x = x,
      focus = "eigengene",
      corr_method = genfocus_corr_method,
      input_scale = "tpm",
      corr_thr = thr,
      CVR_thr = genfocus_CVR_thr
    )
    if (!is.null(gf)) {
      used_thr <- thr
      break
    }
  }
  if (is.null(gf)) {
    # Fallback: if GenFocus cannot establish INGS, run FGM clustering on the
    # original CGM expression matrix to keep the second-stage decomposition usable.
    fallback <- find_cgm(
      x = x,
      mode = mode,
      corr_method = "pearson",
      min_cluster_size = min_cluster_size,
      n_components = as.integer(round(log2(nrow(x)), 0) + 1),
      seed = seed
    )
    fallback <- c(
      fallback,
      list(
        genfocus = NULL,
        genfocus_corr_thr_used = NA_real_
      )
    )
    return(fallback)
  }

  cgm2 <- find_cgm(
    x = gf$normalized_matrix,
    mode = mode,
    corr_method = "pearson",
    min_cluster_size = min_cluster_size,
    n_components = as.integer(round(log2(nrow(x)), 0) + 1),
    seed = seed
  )
  cgm2$genfocus <- gf
  cgm2$genfocus_corr_thr_used <- used_thr
  cgm2
}

#' Run NestedWGCNA end-to-end
#'
#' @param x Expression matrix (`samples x genes`).
#' @param mode Similarity mode.
#' @param min_cgm_size Minimum CGM cluster size.
#' @param min_fgm_size Minimum FGM cluster size.
#' @param corr_method Correlation method for CGM stage.
#' @param genfocus_corr_method Correlation method for GenFocus.
#' @param target_cgm Optional target CGM id for FGM stage. If `NULL`, uses
#'   largest non-noise CGM.
#' @param top_n_genes Optional variable-gene filter size.
#' @param seed Random seed.
#'
#' @return Object of class `NestedWGCNAResult`.
#' @export
run_nested_wgcna <- function(
    x,
    mode = c("paper", "python_compat"),
    min_cgm_size = 100,
    min_fgm_size = 10,
    corr_method = c("pearson", "spearman"),
    genfocus_corr_method = c("spearman", "pearson"),
    target_cgm = NULL,
    top_n_genes = NULL,
    seed = 42) {
  mode <- match.arg(mode)
  corr_method <- match.arg(corr_method)
  genfocus_corr_method <- match.arg(genfocus_corr_method)

  x0 <- as_expression_matrix(x, require_names = FALSE)
  x1 <- clear_data(x0)
  if (!is.null(top_n_genes)) {
    x1 <- top_variable_genes(x1, n = as.integer(top_n_genes))
  }

  cgm <- find_cgm(
    x = x1,
    mode = mode,
    corr_method = corr_method,
    min_cluster_size = min_cgm_size,
    seed = seed
  )
  cgm_assign <- cgm$clusters

  if (is.null(target_cgm)) {
    non_noise <- cgm_assign[cgm_assign >= 0]
    if (length(non_noise) == 0) {
      target_cgm <- -1L
    } else {
      target_cgm <- as.integer(names(sort(table(non_noise), decreasing = TRUE)[1]))
    }
  }

  if (target_cgm < 0) {
    fgm <- list(
      clusters = rep(-1L, ncol(x1)),
      core_clusters = rep(-1L, ncol(x1)),
      embedding = matrix(NA_real_, nrow = ncol(x1), ncol = 2),
      genfocus = NULL,
      genfocus_corr_thr_used = NA_real_
    )
    fgm_assign_full <- rep(-1L, ncol(x1))
    fgm_core_full <- rep(-1L, ncol(x1))
    names(fgm_assign_full) <- colnames(x1)
    names(fgm_core_full) <- colnames(x1)
  } else {
    target_idx <- which(cgm_assign == target_cgm)
    target_expr <- x1[, target_idx, drop = FALSE]
    fgm <- find_fgm(
      x = target_expr,
      mode = mode,
      min_cluster_size = min_fgm_size,
      genfocus_corr_method = genfocus_corr_method,
      seed = seed
    )
    fgm_assign_full <- rep(-1L, ncol(x1))
    fgm_core_full <- rep(-1L, ncol(x1))
    fgm_assign_full[target_idx] <- fgm$clusters
    fgm_core_full[target_idx] <- fgm$core_clusters
    names(fgm_assign_full) <- colnames(x1)
    names(fgm_core_full) <- colnames(x1)
  }

  cgm_scores <- module_score(x1, cgm_assign, prefix = "CGM_")
  fgm_scores <- module_score(x1, fgm_assign_full, prefix = "FGM_")

  res <- list(
    input_info = list(n_samples = nrow(x1), n_genes = ncol(x1)),
    params = list(
      mode = mode,
      min_cgm_size = min_cgm_size,
      min_fgm_size = min_fgm_size,
      corr_method = corr_method,
      genfocus_corr_method = genfocus_corr_method,
      genfocus_corr_thr_used = fgm[["genfocus_corr_thr_used", exact = TRUE]],
      target_cgm = target_cgm,
      top_n_genes = top_n_genes
    ),
    cgm = list(
      assignment = stats::setNames(cgm_assign, colnames(x1)),
      core_assignment = stats::setNames(cgm$core_clusters, colnames(x1)),
      embedding = cgm$embedding
    ),
    fgm = list(
      assignment = fgm_assign_full,
      core_assignment = fgm_core_full,
      embedding = fgm$embedding
    ),
    genfocus = fgm[["genfocus", exact = TRUE]],
    module_scores = list(cgm = cgm_scores, fgm = fgm_scores),
    call = match.call(),
    seed = seed
  )
  class(res) <- "NestedWGCNAResult"
  res
}
