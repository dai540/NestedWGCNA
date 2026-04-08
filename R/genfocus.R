row_means <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowMeans2(as.matrix(x), na.rm = TRUE))
  }
  rowMeans(x, na.rm = TRUE)
}

#' GenFocus normalization (R rewrite of Python reference logic)
#'
#' @param x Expression matrix (`samples x genes`).
#' @param focus `"eigengene"` or `"gene"`.
#' @param focus_gene Focus gene name when `focus = "gene"`.
#' @param corr_method Correlation method.
#' @param input_scale `"tpm"` or `"fpkm"`.
#' @param corr_thr Correlation threshold for INGS selection.
#' @param CVR_thr CVR threshold for clipping.
#'
#' @return `NULL` if INGS cannot be established; otherwise a list with
#'   normalized matrix and diagnostics.
#' @export
genfocus_normalize <- function(
    x,
    focus = c("eigengene", "gene"),
    focus_gene = NULL,
    corr_method = c("spearman", "pearson"),
    input_scale = c("tpm", "fpkm"),
    corr_thr = 0.9,
    CVR_thr = 0.6) {
  focus <- match.arg(focus)
  corr_method <- match.arg(corr_method)
  input_scale <- match.arg(input_scale)
  x <- as_expression_matrix(x, require_names = FALSE)

  df_tpm <- x
  if (input_scale == "fpkm") {
    totals <- rowSums(df_tpm, na.rm = TRUE)
    totals[totals == 0] <- NA_real_
    df_tpm <- (df_tpm / totals) * 1e6
    df_tpm[!is.finite(df_tpm)] <- 0
  }

  if (focus == "eigengene") {
    z <- scale(df_tpm, center = TRUE, scale = TRUE)
    z[!is.finite(z)] <- 0
    pc <- stats::prcomp(z, center = FALSE, scale. = FALSE)
    eig <- pc$x[, 1]
    cor_vec <- apply(
      df_tpm, 2,
      function(g) stats::cor(g, eig, method = corr_method, use = "pairwise.complete.obs")
    )
  } else {
    if (is.null(focus_gene) || !focus_gene %in% colnames(df_tpm)) {
      stop("`focus_gene` is missing from matrix columns.", call. = FALSE)
    }
    cor_vec <- apply(
      df_tpm, 2,
      function(g) stats::cor(g, df_tpm[, focus_gene], method = corr_method, use = "pairwise.complete.obs")
    )
  }
  cor_vec[!is.finite(cor_vec)] <- 0

  ings <- names(cor_vec[cor_vec > corr_thr])
  if (length(ings) <= 1) {
    return(NULL)
  }

  df_norm <- df_tpm

  # First-pass normalization.
  fings <- row_means(df_tpm[, ings, drop = FALSE])
  non_ings <- setdiff(colnames(df_tpm), ings)
  if (length(non_ings)) {
    df_norm[, non_ings] <- df_norm[, non_ings, drop = FALSE] / fings
  }
  for (g in ings) {
    others <- setdiff(ings, g)
    if (length(others) == 0) {
      next
    }
    fing <- row_means(df_tpm[, others, drop = FALSE])
    df_norm[, g] <- df_norm[, g] / fing
  }

  # CV / CVR.
  cv_before <- rep(NA_real_, length(ings))
  names(cv_before) <- ings
  cv_after <- rep(NA_real_, length(ings))
  names(cv_after) <- ings
  cvr <- rep(NA_real_, length(ings))
  names(cvr) <- ings

  for (g in ings) {
    m1 <- mean(df_tpm[, g], na.rm = TRUE)
    m2 <- mean(df_norm[, g], na.rm = TRUE)
    if (is.finite(m1) && is.finite(m2) && m1 != 0 && m2 != 0) {
      cv_before[g] <- stats::sd(df_tpm[, g], na.rm = TRUE) / m1
      cv_after[g] <- stats::sd(df_norm[, g], na.rm = TRUE) / m2
      cvr[g] <- cv_after[g] / cv_before[g]
    }
  }

  cvr_good <- names(cvr[is.finite(cvr) & cvr < CVR_thr])
  ings_final <- intersect(ings, cvr_good)
  if (length(ings_final) <= 1) {
    return(NULL)
  }

  # Final normalization.
  out <- df_tpm
  fings2 <- row_means(df_tpm[, ings_final, drop = FALSE])
  non_ings2 <- setdiff(colnames(df_tpm), ings_final)
  if (length(non_ings2)) {
    out[, non_ings2] <- out[, non_ings2, drop = FALSE] / fings2
  }
  for (g in ings_final) {
    others <- setdiff(ings_final, g)
    fing <- row_means(df_tpm[, others, drop = FALSE])
    out[, g] <- out[, g] / fing
  }

  out[!is.finite(out)] <- 0
  list(
    ings_genes = ings_final,
    normalized_matrix = out,
    correlation = cor_vec,
    CV = cv_before,
    CVR = cvr,
    focus_type = focus,
    focus_gene = if (focus == "gene") focus_gene else "eigengene"
  )
}
