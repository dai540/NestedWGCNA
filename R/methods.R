#' @export
print.NestedWGCNAResult <- function(x, ...) {
  cat("<NestedWGCNAResult>\n", sep = "")
  cat("Samples:", x$input$n_samples, "\n")
  cat("Genes:", x$input$n_genes, "\n")
  cat("CGM modules:", length(unique(x$cgm$assignment$module[x$cgm$assignment$module > 0])), "\n")
  cat("FGM modules:", length(unique(x$fgm$assignment$module[x$fgm$assignment$module > 0])), "\n")
  invisible(x)
}

#' Summarize nested analysis result
#'
#' @param object `NestedWGCNAResult` object.
#' @param ... Unused.
#' @return A summary list.
#' @export
summary.NestedWGCNAResult <- function(object, ...) {
  cgm_modules <- sort(unique(object$cgm$assignment$module[object$cgm$assignment$module > 0]))
  fgm_modules <- sort(unique(object$fgm$assignment$module[object$fgm$assignment$module > 0]))
  core_n <- nrow(object$cgm_core)

  out <- list(
    n_samples = object$input$n_samples,
    n_genes = object$input$n_genes,
    n_cgm_modules = length(cgm_modules),
    n_fgm_modules = length(fgm_modules),
    n_core_genes = core_n
  )
  class(out) <- "summary.NestedWGCNAResult"
  out
}

#' @export
print.summary.NestedWGCNAResult <- function(x, ...) {
  cat("NestedWGCNA summary\n")
  cat("Samples:", x$n_samples, "\n")
  cat("Genes:", x$n_genes, "\n")
  cat("CGM modules:", x$n_cgm_modules, "\n")
  cat("FGM modules:", x$n_fgm_modules, "\n")
  cat("Core genes:", x$n_core_genes, "\n")
  invisible(x)
}
