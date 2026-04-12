# Internal helper to prepare module-score correlation matrices.
prepare_module_score_cor <- function(score_df, top_n = 12, corr_method = c("spearman", "pearson")) {
  corr_method <- match.arg(corr_method)

  if (is.null(score_df) || !nrow(score_df) || ncol(score_df) < 2) {
    return(NULL)
  }

  s <- as.matrix(score_df)
  storage.mode(s) <- "double"

  vars <- apply(s, 2, stats::var, na.rm = TRUE)
  vars[!is.finite(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(min(as.integer(top_n), ncol(s)))]
  s <- s[, keep, drop = FALSE]

  r <- suppressWarnings(stats::cor(s, method = corr_method, use = "pairwise.complete.obs"))
  if (is.null(r) || nrow(r) < 2) {
    return(NULL)
  }
  r[!is.finite(r)] <- 0
  diag(r) <- 0

  list(
    matrix = s,
    cor = r,
    modules = colnames(s)
  )
}

# Determine an edge threshold with fallback so the graph is not empty.
resolve_edge_threshold <- function(r, cor_thr = 0.3, quantile_fallback = 0.8) {
  upper <- abs(r[upper.tri(r)])
  if (!length(upper)) {
    return(1)
  }

  thr <- as.numeric(cor_thr)
  if (!is.finite(thr) || thr <= 0) {
    thr <- as.numeric(stats::quantile(upper, probs = quantile_fallback, na.rm = TRUE))
  }

  edge_idx <- which(abs(r) >= thr, arr.ind = TRUE)
  edge_idx <- edge_idx[edge_idx[, 1] < edge_idx[, 2], , drop = FALSE]

  if (!nrow(edge_idx)) {
    thr <- max(0.2, as.numeric(stats::quantile(upper, probs = quantile_fallback, na.rm = TRUE)))
  }

  thr
}

#' Plot a module-score network
#'
#' Draw a clean module network from module-score correlations. The function is
#' designed for NestedWGCNA CGM/FGM score tables and follows WGCNA-style
#' eigengene-network visualization principles.
#'
#' @param score_df Data frame or matrix (`samples x modules`) of module scores.
#' @param module_sizes Optional named numeric vector of module sizes. Names must
#'   match `colnames(score_df)` (for example `CGM_0`, `FGM_2`).
#' @param layer Layer label used in the title.
#' @param top_n Maximum number of modules (highest variance) to include.
#' @param cor_thr Correlation threshold for drawing edges.
#' @param corr_method Correlation method.
#' @param layout_seed Seed used by force-directed layouts.
#' @param show_legend Whether to draw edge-sign legend.
#'
#' @return Invisibly returns the plotted correlation matrix.
#' @export
plot_module_network <- function(
    score_df,
    module_sizes = NULL,
    layer = "CGM",
    top_n = 12,
    cor_thr = 0.3,
    corr_method = c("spearman", "pearson"),
    layout_seed = 42,
    show_legend = TRUE) {
  graphics::par(bg = "white", fg = "#111827", col.main = "#111827", col.lab = "#111827", col.axis = "#111827")

  corr_method <- match.arg(corr_method)
  prep <- prepare_module_score_cor(score_df, top_n = top_n, corr_method = corr_method)
  if (is.null(prep)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste(layer, "network unavailable"))
    return(invisible(NULL))
  }

  r <- prep$cor
  mods <- prep$modules
  thr <- resolve_edge_threshold(r, cor_thr = cor_thr)

  adj <- abs(r) >= thr
  diag(adj) <- FALSE

  node_colors <- grDevices::hcl.colors(length(mods), palette = "Dark 3")
  names(node_colors) <- mods

  default_sizes <- rep(16, length(mods))
  names(default_sizes) <- mods
  if (!is.null(module_sizes)) {
    ms <- module_sizes[mods]
    ms[!is.finite(ms)] <- NA_real_
    if (all(is.na(ms))) {
      node_sizes <- default_sizes
    } else {
      ms0 <- ms
      rng <- range(ms0, na.rm = TRUE)
      if (diff(rng) == 0) {
        node_sizes <- rep(20, length(ms0))
      } else {
        node_sizes <- 12 + 18 * (ms0 - rng[1]) / diff(rng)
      }
      node_sizes[!is.finite(node_sizes)] <- 14
      names(node_sizes) <- mods
    }
  } else {
    node_sizes <- default_sizes
  }

  if (requireNamespace("igraph", quietly = TRUE)) {
    w <- matrix(0, nrow = nrow(r), ncol = ncol(r))
    w[adj] <- abs(r[adj])
    w <- (w + t(w)) / 2

    g <- igraph::graph_from_adjacency_matrix(
      w,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    if (igraph::ecount(g) == 0) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste(layer, "network has no edges at threshold"))
      return(invisible(r))
    }

    el <- igraph::as_edgelist(g, names = FALSE)
    edge_sign <- sign(r[cbind(el[, 1], el[, 2])])
    edge_col <- ifelse(edge_sign >= 0, "#1f77b4AA", "#d62728AA")
    edge_w <- 1 + 5 * igraph::E(g)$weight

    set.seed(as.integer(layout_seed))
    xy <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight, niter = 800)

    graphics::plot(
      g,
      layout = xy,
      vertex.label = mods,
      vertex.label.cex = 0.78,
      vertex.label.color = "#111827",
      vertex.label.family = "sans",
      vertex.size = as.numeric(node_sizes[mods]),
      vertex.color = as.character(node_colors[mods]),
      vertex.frame.color = "#1f2937",
      edge.color = edge_col,
      edge.width = edge_w,
      main = paste0(layer, " module-score network")
    )
    if (isTRUE(show_legend)) {
      graphics::legend(
        "topright",
        legend = c("positive", "negative"),
        col = c("#1f77b4", "#d62728"),
        text.col = "#111827",
        lwd = 2,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(r))
  }

  # Fallback when igraph is unavailable.
  edge_idx <- which(adj, arr.ind = TRUE)
  edge_idx <- edge_idx[edge_idx[, 1] < edge_idx[, 2], , drop = FALSE]
  d <- stats::as.dist(1 - abs(r))
  xy <- try(stats::cmdscale(d, k = 2), silent = TRUE)
  if (inherits(xy, "try-error") || any(!is.finite(xy))) {
    th <- seq(0, 2 * pi, length.out = length(mods) + 1)[-1]
    xy <- cbind(cos(th), sin(th))
  }
  rownames(xy) <- mods

  graphics::plot(
    xy[, 1],
    xy[, 2],
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = paste0(layer, " module-score network")
  )
  if (nrow(edge_idx)) {
    for (i in seq_len(nrow(edge_idx))) {
      a <- edge_idx[i, 1]
      b <- edge_idx[i, 2]
      rc <- r[a, b]
      graphics::segments(
        xy[a, 1], xy[a, 2], xy[b, 1], xy[b, 2],
        col = if (rc >= 0) "#1f77b4AA" else "#d62728AA",
        lwd = 1 + 4 * abs(rc)
      )
    }
  }
  graphics::points(
    xy[, 1],
    xy[, 2],
    pch = 21,
    cex = pmax(0.8, as.numeric(node_sizes[mods]) / 9),
    bg = as.character(node_colors[mods]),
    col = "#1f2937"
  )
  graphics::text(xy[, 1], xy[, 2], labels = mods, pos = 3, cex = 0.78, col = "#111827")
  if (isTRUE(show_legend)) {
    graphics::legend(
      "topright",
      legend = c("positive", "negative"),
      col = c("#1f77b4", "#d62728"),
      text.col = "#111827",
      lwd = 2,
      bty = "n",
      cex = 0.85
    )
  }

  invisible(r)
}

#' Plot a module-score correlation heatmap (WGCNA-style)
#'
#' Plot a dendrogram-ordered module-score correlation heatmap, analogous to
#' module eigengene network heatmaps commonly used in WGCNA tutorials.
#'
#' @param score_df Data frame or matrix (`samples x modules`) of module scores.
#' @param layer Layer label used in the title.
#' @param top_n Maximum number of modules (highest variance) to include.
#' @param corr_method Correlation method.
#'
#' @return Invisibly returns the plotted correlation matrix.
#' @export
plot_module_score_heatmap <- function(
    score_df,
    layer = "CGM",
    top_n = 15,
    corr_method = c("spearman", "pearson")) {
  graphics::par(bg = "white", fg = "#111827", col.main = "#111827", col.lab = "#111827", col.axis = "#111827")

  corr_method <- match.arg(corr_method)
  prep <- prepare_module_score_cor(score_df, top_n = top_n, corr_method = corr_method)
  if (is.null(prep)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste(layer, "heatmap unavailable"))
    return(invisible(NULL))
  }

  r <- prep$cor
  d <- stats::as.dist(1 - abs(r))
  hc <- stats::hclust(d, method = "average")
  ord <- hc$order
  r_ord <- r[ord, ord, drop = FALSE]

  cols <- grDevices::colorRampPalette(
    c("#2c7bb6", "#abd9e9", "#f7f7f7", "#fdae61", "#d7191c")
  )(101)

  stats::heatmap(
    r_ord,
    Rowv = stats::as.dendrogram(hc),
    Colv = stats::as.dendrogram(hc),
    scale = "none",
    symm = TRUE,
    col = cols,
    margins = c(8, 8),
    cexRow = 0.8,
    cexCol = 0.8,
    main = paste0(layer, " module-score correlation heatmap")
  )

  invisible(r_ord)
}
