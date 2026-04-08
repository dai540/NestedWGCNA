#' @export
print.NestedWGCNAResult <- function(x, ...) {
  cat("<NestedWGCNAResult>\n")
  cat("samples:", x$input_info$n_samples, "\n")
  cat("genes:", x$input_info$n_genes, "\n")
  cat("mode:", x$params$mode, "\n")
  cat("CGM modules:", length(unique(x$cgm$assignment[x$cgm$assignment >= 0])), "\n")
  cat("FGM modules:", length(unique(x$fgm$assignment[x$fgm$assignment >= 0])), "\n")
  invisible(x)
}

#' Summarize a NestedWGCNA result
#'
#' @param object A `NestedWGCNAResult`.
#' @param ... Unused.
#'
#' @return A named list summary.
#' @export
summary.NestedWGCNAResult <- function(object, ...) {
  list(
    n_samples = object$input_info$n_samples,
    n_genes = object$input_info$n_genes,
    mode = object$params$mode,
    cgm_n_modules = length(unique(object$cgm$assignment[object$cgm$assignment >= 0])),
    fgm_n_modules = length(unique(object$fgm$assignment[object$fgm$assignment >= 0])),
    target_cgm = object$params$target_cgm
  )
}
