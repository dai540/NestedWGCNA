.cluster_by_hclust <- function(dissimilarity, min_module_size) {
  n_genes <- ncol(dissimilarity)
  if (n_genes < 2L) {
    out <- rep(1L, n_genes)
    names(out) <- colnames(dissimilarity)
    return(out)
  }

  tree <- stats::hclust(stats::as.dist(dissimilarity), method = "average")
  k <- max(2L, floor(n_genes / max(2L, as.integer(min_module_size))))
  k <- min(k, n_genes)
  modules <- stats::cutree(tree, k = k)
  modules <- as.integer(modules)
  names(modules) <- colnames(dissimilarity)

  size_tbl <- table(modules)
  small_ids <- as.integer(names(size_tbl[size_tbl < as.integer(min_module_size)]))
  if (length(small_ids) > 0L) {
    modules[modules %in% small_ids] <- 0L
  }

  positive <- sort(unique(modules[modules > 0L]))
  if (length(positive) == 0L) {
    modules[] <- 1L
    return(modules)
  }

  map <- stats::setNames(seq_along(positive), positive)
  relabeled <- integer(length(modules))
  relabeled[modules > 0L] <- unname(map[as.character(modules[modules > 0L])])
  names(relabeled) <- names(modules)
  relabeled
}

.select_core_genes <- function(adjacency, modules, core_fraction = 0.25, min_core = 2L) {
  module_ids <- sort(unique(modules[modules > 0L]))
  core_list <- vector("list", length(module_ids))
  names(core_list) <- as.character(module_ids)

  for (id in module_ids) {
    genes <- names(modules)[modules == id]
    if (length(genes) == 0L) {
      core_list[[as.character(id)]] <- character(0)
      next
    }
    if (length(genes) == 1L) {
      core_list[[as.character(id)]] <- genes
      next
    }
    sub_adj <- adjacency[genes, genes, drop = FALSE]
    connectivity <- rowSums(sub_adj) - 1
    n_core <- min(length(genes), max(as.integer(min_core), floor(length(genes) * core_fraction)))
    core_list[[as.character(id)]] <- names(sort(connectivity, decreasing = TRUE))[seq_len(n_core)]
  }

  core_df <- do.call(
    rbind,
    lapply(names(core_list), function(id) {
      genes <- core_list[[id]]
      if (length(genes) == 0L) {
        return(NULL)
      }
      data.frame(
        gene = genes,
        module = as.integer(id),
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(core_df)) {
    core_df <- data.frame(gene = character(0), module = integer(0))
  }

  list(core_genes_by_module = core_list, assignment = core_df)
}

.normalize_by_core <- function(x, modules, core_genes_by_module) {
  x_norm <- x
  module_ids <- sort(unique(modules[modules > 0L]))

  for (id in module_ids) {
    genes <- names(modules)[modules == id]
    cores <- core_genes_by_module[[as.character(id)]]
    cores <- intersect(cores, colnames(x_norm))
    genes <- intersect(genes, colnames(x_norm))
    if (length(cores) == 0L || length(genes) == 0L) {
      next
    }
    baseline <- rowMeans(x_norm[, cores, drop = FALSE])
    x_norm[, genes] <- sweep(x_norm[, genes, drop = FALSE], 1L, baseline, "-")
  }

  x_norm
}

#' Find coarse-grained modules (CGM)
#'
#' @param x Numeric matrix-like object.
#' @param mode Similarity mode.
#' @param min_cgm_size Minimum module size.
#' @param top_n_genes Number of high-variance genes used for CGM.
#' @param corr_method Correlation method.
#' @param keep_matrices If `TRUE`, include adjacency and dissimilarity in output.
#' @return List with CGM assignments.
#' @export
find_cgm <- function(
  x,
  mode = c("paper", "python_compat"),
  min_cgm_size = 100L,
  top_n_genes = 2000L,
  corr_method = c("pearson", "spearman"),
  keep_matrices = FALSE
) {
  mode <- match.arg(mode)
  corr_method <- match.arg(corr_method)

  x_use <- top_variable_genes(x, n = top_n_genes)
  adj <- compute_adjacency(x_use, mode = mode, method = corr_method)
  dis <- compute_dissimilarity(adj)
  modules <- .cluster_by_hclust(dis, min_module_size = min_cgm_size)

  out <- list(
    assignment = data.frame(
      gene = names(modules),
      module = as.integer(modules),
      stringsAsFactors = FALSE
    ),
    module_vector = modules,
    matrix = x_use
  )
  if (isTRUE(keep_matrices)) {
    out$adjacency <- adj
    out$dissimilarity <- dis
  }
  out
}

#' Find fine-grained modules (FGM)
#'
#' @param x Numeric matrix-like object.
#' @param mode Similarity mode.
#' @param min_fgm_size Minimum module size.
#' @param corr_method Correlation method.
#' @param keep_matrices If `TRUE`, include adjacency and dissimilarity in output.
#' @return List with FGM assignments.
#' @export
find_fgm <- function(
  x,
  mode = c("paper", "python_compat"),
  min_fgm_size = 20L,
  corr_method = c("pearson", "spearman"),
  keep_matrices = FALSE
) {
  mode <- match.arg(mode)
  corr_method <- match.arg(corr_method)

  x_use <- clear_data(x, min_samples = 2L)
  adj <- compute_adjacency(x_use, mode = mode, method = corr_method)
  dis <- compute_dissimilarity(adj)
  modules <- .cluster_by_hclust(dis, min_module_size = min_fgm_size)

  out <- list(
    assignment = data.frame(
      gene = names(modules),
      module = as.integer(modules),
      stringsAsFactors = FALSE
    ),
    module_vector = modules
  )
  if (isTRUE(keep_matrices)) {
    out$adjacency <- adj
    out$dissimilarity <- dis
  }
  out
}

#' Compute module scores
#'
#' @param x Numeric matrix-like object.
#' @param assignments Named integer vector or data frame with `gene` and `module`.
#' @return Sample-by-module score matrix.
#' @export
module_score <- function(x, assignments) {
  x <- as_expression_matrix(x)

  if (is.data.frame(assignments)) {
    if (!all(c("gene", "module") %in% colnames(assignments))) {
      stop("`assignments` data.frame must contain `gene` and `module`.")
    }
    modules <- stats::setNames(as.integer(assignments$module), assignments$gene)
  } else {
    modules <- as.integer(assignments)
    names(modules) <- names(assignments)
  }

  modules <- modules[modules > 0L]
  modules <- modules[names(modules) %in% colnames(x)]
  module_ids <- sort(unique(modules))

  if (length(module_ids) == 0L) {
    return(matrix(numeric(0), nrow = nrow(x), ncol = 0L))
  }

  score <- matrix(NA_real_, nrow = nrow(x), ncol = length(module_ids))
  rownames(score) <- rownames(x)
  colnames(score) <- paste0("M", module_ids)

  for (i in seq_along(module_ids)) {
    id <- module_ids[i]
    genes <- names(modules)[modules == id]
    score[, i] <- rowMeans(x[, genes, drop = FALSE])
  }

  score
}

#' Run nested two-stage workflow
#'
#' @param x Numeric matrix-like object (`samples x genes`).
#' @param mode Similarity mode.
#' @param min_cgm_size Minimum CGM size.
#' @param min_fgm_size Minimum FGM size.
#' @param top_n_genes Number of genes used in stage 1.
#' @param corr_method Correlation method.
#' @return Object of class `NestedWGCNAResult`.
#' @export
run_nested_wgcna <- function(
  x,
  mode = c("paper", "python_compat"),
  min_cgm_size = 100L,
  min_fgm_size = 20L,
  top_n_genes = 2000L,
  corr_method = c("pearson", "spearman")
) {
  mode <- match.arg(mode)
  corr_method <- match.arg(corr_method)
  x_clean <- clear_data(x, min_samples = 2L)

  cgm <- find_cgm(
    x = x_clean,
    mode = mode,
    min_cgm_size = min_cgm_size,
    top_n_genes = top_n_genes,
    corr_method = corr_method,
    keep_matrices = TRUE
  )

  core <- .select_core_genes(cgm$adjacency, cgm$module_vector)
  x_norm <- .normalize_by_core(cgm$matrix, cgm$module_vector, core$core_genes_by_module)

  fgm <- find_fgm(
    x = x_norm,
    mode = mode,
    min_fgm_size = min_fgm_size,
    corr_method = corr_method,
    keep_matrices = FALSE
  )

  out <- list(
    input = list(n_samples = nrow(x_clean), n_genes = ncol(x_clean)),
    params = list(
      mode = mode,
      min_cgm_size = as.integer(min_cgm_size),
      min_fgm_size = as.integer(min_fgm_size),
      top_n_genes = as.integer(top_n_genes),
      corr_method = corr_method
    ),
    cgm = list(assignment = cgm$assignment),
    cgm_core = core$assignment,
    normalized_matrix = x_norm,
    fgm = list(assignment = fgm$assignment),
    module_scores = module_score(x_norm, fgm$module_vector),
    call = match.call()
  )
  class(out) <- "NestedWGCNAResult"
  out
}
