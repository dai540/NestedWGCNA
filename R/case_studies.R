#' List bundled case-study ids
#'
#' @return Character vector of case-study dataset ids.
#' @export
available_case_studies <- function() {
  c("tcga_blca", "all_leukemia", "bladderbatch")
}

#' Read case-study summary JSON
#'
#' @param dataset Dataset id. One of `available_case_studies()`.
#'
#' @return A named list parsed from `summary.json`.
#' @export
case_study_summary <- function(dataset = "tcga_blca") {
  if (!dataset %in% available_case_studies()) {
    stop("Unknown dataset: ", dataset, call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required.", call. = FALSE)
  }
  p <- system.file("extdata", "case_studies", dataset, "summary.json", package = "NestedWGCNA")
  if (p == "") {
    stop("summary.json not found for dataset: ", dataset, call. = FALSE)
  }
  jsonlite::fromJSON(p)
}

#' Read case-study output table
#'
#' @param dataset Dataset id. One of `available_case_studies()`.
#' @param table_name Output table basename without extension.
#'
#' @return A data frame.
#' @export
case_study_table <- function(dataset = "tcga_blca", table_name = "phenotype_associations") {
  if (!dataset %in% available_case_studies()) {
    stop("Unknown dataset: ", dataset, call. = FALSE)
  }
  valid <- c(
    "cgm_assignments",
    "fgm_assignments",
    "cgm_enrichment",
    "fgm_enrichment",
    "cgm_scores",
    "fgm_scores",
    "phenotype_associations"
  )
  if (!table_name %in% valid) {
    stop("Unknown table_name: ", table_name, call. = FALSE)
  }
  p <- system.file(
    "extdata",
    "case_studies",
    dataset,
    paste0(table_name, ".tsv"),
    package = "NestedWGCNA"
  )
  if (p == "") {
    stop("Table not found: ", table_name, " for dataset: ", dataset, call. = FALSE)
  }
  utils::read.delim(p, check.names = FALSE)
}
