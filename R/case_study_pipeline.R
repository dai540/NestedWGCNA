read_signature_sets <- function(path, type = c("xcell_csv", "bostongene_tsv")) {
  type <- match.arg(type)
  if (!file.exists(path)) {
    stop("Signature file not found: ", path, call. = FALSE)
  }
  out <- list()
  if (type == "xcell_csv") {
    lines <- readLines(path, warn = FALSE)
    if (length(lines) <= 1) {
      return(out)
    }
    for (i in 2:length(lines)) {
      parts <- trimws(strsplit(lines[[i]], ",", fixed = TRUE)[[1]])
      if (length(parts) < 3 || parts[[1]] == "") {
        next
      }
      genes <- parts[3:length(parts)]
      genes <- genes[nzchar(genes)]
      out[[parts[[1]]]] <- unique(genes)
    }
  } else {
    dat <- utils::read.delim(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(dat))) {
      sig <- trimws(dat[i, 1])
      if (!nzchar(sig)) {
        next
      }
      genes <- trimws(as.character(unlist(dat[i, -(1:2), drop = FALSE])))
      genes <- genes[nzchar(genes)]
      out[[sig]] <- unique(genes)
    }
  }
  out
}

bh_adjust <- function(p) {
  stats::p.adjust(p, method = "BH")
}

hypergeom_enrichment <- function(module_genes, signatures, source, universe_genes) {
  module_genes <- unique(module_genes)
  m <- length(universe_genes)
  k <- length(module_genes)
  rows <- list()
  for (sig in names(signatures)) {
    g <- intersect(signatures[[sig]], universe_genes)
    n <- length(g)
    if (n < 5) {
      next
    }
    x <- length(intersect(module_genes, g))
    if (x == 0) {
      next
    }
    p <- stats::phyper(q = x - 1, m = n, n = m - n, k = k, lower.tail = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      source = source,
      signature = sig,
      M = m,
      K = k,
      n = n,
      k = x,
      p = p,
      overlap_genes = paste(sort(intersect(module_genes, g)), collapse = "|"),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  out$fdr <- bh_adjust(out$p)
  out[order(out$fdr), , drop = FALSE]
}

associate_module_scores <- function(scores, metadata, layer = c("CGM", "FGM")) {
  layer <- match.arg(layer)
  if (is.null(metadata) || !nrow(metadata) || !ncol(scores)) {
    return(data.frame())
  }
  common <- intersect(rownames(scores), rownames(metadata))
  if (length(common) < 10) {
    return(data.frame())
  }
  scores <- scores[common, , drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]

  rows <- list()
  for (ph in colnames(metadata)) {
    trait <- metadata[[ph]]
    if (is.numeric(trait)) {
      for (mod in colnames(scores)) {
        idx <- which(!is.na(trait) & !is.na(scores[[mod]]))
        if (length(idx) < 10) {
          next
        }
        ct <- suppressWarnings(stats::cor.test(scores[[mod]][idx], trait[idx], method = "spearman"))
        rows[[length(rows) + 1]] <- data.frame(
          layer = layer,
          module = mod,
          phenotype = ph,
          test = "spearman",
          levels_or_n = length(idx),
          effect = unname(ct$estimate),
          p = ct$p.value,
          stringsAsFactors = FALSE
        )
      }
    } else {
      g <- as.character(trait)
      tab <- sort(table(g), decreasing = TRUE)
      keep <- names(tab[tab >= 3])
      if (length(keep) < 2 || length(keep) > 12) {
        next
      }
      idx_keep <- which(g %in% keep)
      for (mod in colnames(scores)) {
        values <- scores[[mod]][idx_keep]
        groups <- g[idx_keep]
        split_vals <- split(values, groups)
        if (length(split_vals) < 2) {
          next
        }
        fit <- try(stats::oneway.test(values ~ as.factor(groups), var.equal = FALSE), silent = TRUE)
        if (inherits(fit, "try-error")) {
          next
        }
        rows[[length(rows) + 1]] <- data.frame(
          layer = layer,
          module = mod,
          phenotype = ph,
          test = "anova",
          levels_or_n = length(split_vals),
          effect = unname(fit$statistic),
          p = fit$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(rows)) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  out$fdr <- bh_adjust(out$p)
  out[order(out$fdr), , drop = FALSE]
}

pick_target_cgm <- function(enrichment_tbl, cgm_assign) {
  immune_pattern <- "B-cell|B_cells|T-cell|T_cells|NK|Macroph|Dendritic|Immune|MHC|Cytotoxic|MDSC|Treg|Th1|Th2"
  if (nrow(enrichment_tbl)) {
    cand <- enrichment_tbl[grepl(immune_pattern, enrichment_tbl$signature, ignore.case = TRUE), , drop = FALSE]
    if (nrow(cand)) {
      best <- cand[order(cand$fdr), , drop = FALSE][1, ]
      return(list(module = as.integer(best$module), reason = paste0("immune_signature:", best$signature)))
    }
  }
  non_noise <- cgm_assign[cgm_assign >= 0]
  if (!length(non_noise)) {
    return(list(module = -1L, reason = "none"))
  }
  k <- as.integer(names(sort(table(non_noise), decreasing = TRUE)[1]))
  list(module = k, reason = "largest_cgm")
}

normalize_for_case_study <- function(x) {
  x <- as_expression_matrix(x, require_names = FALSE)
  min_before <- min(x, na.rm = TRUE)
  shift <- 0
  if (min_before < 0) {
    shift <- -min_before + 1e-6
    x <- x + shift
  }
  q99 <- as.numeric(stats::quantile(x, probs = 0.99, na.rm = TRUE))
  log_applied <- FALSE
  if (q99 > 50) {
    x <- log1p(x)
    log_applied <- TRUE
  }
  list(
    matrix = x,
    meta = list(
      input_min_before_shift = min_before,
      shift_applied = shift,
      log1p_applied = as.numeric(log_applied),
      q99_after_shift_before_log = q99
    )
  )
}

#' Run a complete case-study pipeline in R
#'
#' @param expr Expression matrix (`samples x genes`).
#' @param dataset_id Dataset id string.
#' @param source Source description or URL.
#' @param output_dir Output directory.
#' @param metadata Optional sample metadata data.frame (rownames = sample ids).
#' @param mode NestedWGCNA mode.
#' @param cgm_min_cluster_size CGM minimum cluster size.
#' @param fgm_min_cluster_size FGM minimum cluster size.
#' @param top_n_genes Number of variable genes to keep.
#' @param seed Random seed.
#'
#' @return Summary list (also written to `summary.json`).
#' @export
run_case_study <- function(
    expr,
    dataset_id,
    source,
    output_dir,
    metadata = NULL,
    mode = c("paper", "python_compat"),
    cgm_min_cluster_size = 30,
    fgm_min_cluster_size = 10,
    top_n_genes = 1200,
    seed = 42) {
  mode <- match.arg(mode)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Signatures
  shared_dir <- file.path(dirname(output_dir), "_shared_signatures")
  dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
  xcell_path <- file.path(shared_dir, "xCell_signatures.csv")
  boston_path <- file.path(shared_dir, "bostongene_signatures.txt")
  if (!file.exists(xcell_path)) {
    utils::download.file(
      "https://raw.githubusercontent.com/ilyada/NestedWGCNA/main/data/xCell_signatures.csv",
      xcell_path,
      mode = "wb"
    )
  }
  if (!file.exists(boston_path)) {
    utils::download.file(
      "https://raw.githubusercontent.com/ilyada/NestedWGCNA/main/data/bostongene_signatures.txt",
      boston_path,
      mode = "wb"
    )
  }
  xcell <- read_signature_sets(xcell_path, "xcell_csv")
  boston <- read_signature_sets(boston_path, "bostongene_tsv")

  # Normalize + model genes
  n0 <- ncol(expr)
  prep <- normalize_for_case_study(expr)
  x <- prep$matrix
  x <- top_variable_genes(x, n = top_n_genes)

  # Stage 1 to select target module.
  cgm_raw <- find_cgm(
    x = x,
    mode = mode,
    corr_method = "pearson",
    min_cluster_size = cgm_min_cluster_size,
    seed = seed
  )
  cgm_assign <- stats::setNames(as.integer(cgm_raw$clusters), colnames(x))
  universe <- colnames(x)

  cgm_enr <- data.frame()
  for (k in sort(unique(cgm_assign[cgm_assign >= 0]))) {
    genes <- names(cgm_assign)[cgm_assign == k]
    tmp <- rbind(
      hypergeom_enrichment(genes, xcell, "xCell", universe),
      hypergeom_enrichment(genes, boston, "BostonGene", universe)
    )
    if (nrow(tmp)) {
      tmp$module <- k
      cgm_enr <- rbind(cgm_enr, tmp)
    }
  }
  target <- pick_target_cgm(cgm_enr, cgm_assign)

  # Full nested run with target module fixed.
  res <- run_nested_wgcna(
    x = x,
    mode = mode,
    min_cgm_size = cgm_min_cluster_size,
    min_fgm_size = fgm_min_cluster_size,
    target_cgm = target$module,
    seed = seed
  )

  fgm_assign <- res$fgm$assignment
  fgm_enr <- data.frame()
  for (k in sort(unique(fgm_assign[fgm_assign >= 0]))) {
    genes <- names(fgm_assign)[fgm_assign == k]
    tmp <- rbind(
      hypergeom_enrichment(genes, xcell, "xCell", universe),
      hypergeom_enrichment(genes, boston, "BostonGene", universe)
    )
    if (nrow(tmp)) {
      tmp$fgm <- k
      fgm_enr <- rbind(fgm_enr, tmp)
    }
  }

  # Associations
  assoc_cgm <- associate_module_scores(res$module_scores$cgm, metadata, "CGM")
  assoc_fgm <- associate_module_scores(res$module_scores$fgm, metadata, "FGM")
  assoc <- rbind(assoc_cgm, assoc_fgm)

  # Save outputs
  utils::write.table(
    data.frame(gene = names(cgm_assign), cgm = cgm_assign, cgm_core = res$cgm$core_assignment),
    file = file.path(output_dir, "cgm_assignments.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(gene = names(fgm_assign), fgm = fgm_assign, fgm_core = res$fgm$core_assignment),
    file = file.path(output_dir, "fgm_assignments.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(cgm_enr, file = file.path(output_dir, "cgm_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(fgm_enr, file = file.path(output_dir, "fgm_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(res$module_scores$cgm, file = file.path(output_dir, "cgm_scores.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  utils::write.table(res$module_scores$fgm, file = file.path(output_dir, "fgm_scores.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  utils::write.table(assoc, file = file.path(output_dir, "phenotype_associations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Size plots
  cgm_size <- sort(table(cgm_assign[cgm_assign >= 0]))
  if (length(cgm_size)) {
    grDevices::png(file.path(output_dir, "cgm_module_sizes.png"), width = 1000, height = 620, res = 140)
    graphics::barplot(cgm_size, col = "#2a9d8f", main = paste0("CGM sizes (", dataset_id, ")"), xlab = "CGM id", ylab = "genes")
    grDevices::dev.off()
  }
  fgm_size <- sort(table(fgm_assign[fgm_assign >= 0]))
  if (length(fgm_size)) {
    grDevices::png(file.path(output_dir, "fgm_module_sizes.png"), width = 1000, height = 620, res = 140)
    graphics::barplot(fgm_size, col = "#e76f51", main = paste0("FGM sizes (", dataset_id, ")"), xlab = "FGM id", ylab = "genes")
    grDevices::dev.off()
  }

  summary <- list(
    dataset_id = dataset_id,
    source = source,
    random_state = seed,
    input_meta = prep$meta,
    n_samples = nrow(x),
    n_genes_input = n0,
    n_genes_model = ncol(x),
    cgm_min_cluster_size = cgm_min_cluster_size,
    cgm_n_modules = length(unique(cgm_assign[cgm_assign >= 0])),
    cgm_assigned_genes = sum(cgm_assign >= 0),
    cgm_noise_genes = sum(cgm_assign < 0),
    target_cgm = target$module,
    target_cgm_reason = target$reason,
    target_cgm_size = sum(cgm_assign == target$module),
    genfocus_corr_thr = res$params$genfocus_corr_thr_used,
    genfocus_ings_n = if (is.null(res$genfocus)) 0L else length(res$genfocus$ings_genes),
    fgm_min_cluster_size_requested = fgm_min_cluster_size,
    fgm_min_cluster_size_used = fgm_min_cluster_size,
    fgm_n_modules = length(unique(fgm_assign[fgm_assign >= 0])),
    fgm_assigned_genes = sum(fgm_assign >= 0),
    fgm_noise_genes = sum(fgm_assign < 0)
  )
  if (nrow(cgm_enr)) {
    top <- cgm_enr[order(cgm_enr$fdr), c("module", "source", "signature", "k", "K", "fdr"), drop = FALSE]
    top <- do.call(rbind, lapply(split(top, top$module), function(z) z[seq_len(min(3, nrow(z))), , drop = FALSE]))
    summary$top_cgm_enrichment <- lapply(seq_len(nrow(top)), function(i) as.list(top[i, , drop = FALSE]))
  }
  if (nrow(fgm_enr)) {
    top <- fgm_enr[order(fgm_enr$fdr), c("fgm", "source", "signature", "k", "K", "fdr"), drop = FALSE]
    top <- do.call(rbind, lapply(split(top, top$fgm), function(z) z[seq_len(min(3, nrow(z))), , drop = FALSE]))
    summary$top_fgm_enrichment <- lapply(seq_len(nrow(top)), function(i) as.list(top[i, , drop = FALSE]))
  }
  if (nrow(assoc)) {
    top_assoc <- assoc[order(assoc$fdr), , drop = FALSE]
    top_assoc <- top_assoc[seq_len(min(15, nrow(top_assoc))), c("layer", "module", "phenotype", "test", "effect", "fdr"), drop = FALSE]
    summary$top_phenotype_associations <- lapply(seq_len(nrow(top_assoc)), function(i) as.list(top_assoc[i, , drop = FALSE]))
  }

  jsonlite::write_json(
    summary,
    path = file.path(output_dir, "summary.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null",
    digits = 10
  )
  summary
}
