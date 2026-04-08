# NestedWGCNA

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/NestedWGCNA/)
[![R-CMD-check](https://github.com/dai540/NestedWGCNA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/NestedWGCNA/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/dai540/NestedWGCNA)](https://github.com/dai540/NestedWGCNA/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`NestedWGCNA` is an R package for two-stage gene co-expression network analysis.
It provides an end-to-end workflow to discover coarse-grained modules (CGMs),
normalize within selected module cores (GenFocus), and discover fine-grained
modules (FGMs):

<https://dai540.github.io/NestedWGCNA/>

- `find_cgm()` for coarse-grained module discovery
- `genfocus_normalize()` for core-based normalization
- `find_fgm()` for fine-grained module discovery in a selected CGM

The package is intentionally focused. It does not implement raw-count
normalization or heavy downstream clinical modeling. Instead, it standardizes
the two-stage network workflow, result objects, case-study outputs, and
reproducible tutorials.

It is designed for analysts working with expression matrices in
`samples/cells x genes` format, especially bulk RNA-seq and pseudobulk-style
inputs.

## Installation

Install from GitHub:

```r
install.packages("pak")
pak::pak("dai540/NestedWGCNA")
```

Or:

```r
install.packages("remotes")
remotes::install_github("dai540/NestedWGCNA")
```

Or install from a source tarball:

```r
install.packages("path/to/NestedWGCNA_<version>.tar.gz", repos = NULL, type = "source")
```

Then load the package:

```r
library(NestedWGCNA)
```

Optional dependencies for clustering and tutorials:

```r
install.packages(c("uwot", "dbscan", "matrixStats"))
BiocManager::install(c("ALL", "bladderbatch", "Biobase"))
```

## Citation

If you use `NestedWGCNA`, cite the package as:

> Dai (2026). *NestedWGCNA: Two-Stage Gene Co-Expression Network Analysis*. R package.
> <https://dai540.github.io/NestedWGCNA/>

You can also retrieve the citation from R:

```r
citation("NestedWGCNA")
```

## Documentation

<https://dai540.github.io/NestedWGCNA/>

## What NestedWGCNA does

`NestedWGCNA` does four things.

- Computes adjacency/dissimilarity matrices with explicit reproducibility modes
- Discovers CGMs and core genes using UMAP + HDBSCAN + core decomposition
- Runs nested FGM discovery with GenFocus-aware fallback behavior
- Exports standardized case-study tables, enrichment results, and figures

In practice, the package is doing this:

- `compute_adjacency()` and `compute_dissimilarity()` define the network space
- `find_cgm()` finds CGMs and core assignments
- `find_fgm()` and `run_nested_wgcna()` run stage-2 decomposition
- `module_score()` computes per-sample module scores
- `run_case_study()` generates reproducible real-data outputs

The package supports two main method modes:

- `mode = "paper"`: adjacency `r^2`, dissimilarity `sqrt(1 - r^2)`
- `mode = "python_compat"`: adjacency `|r|`, dissimilarity `sqrt(1 - |r|)`

## Main functions

- `run_nested_wgcna()`
- `find_cgm()`
- `find_fgm()`
- `genfocus_normalize()`
- `compute_adjacency()`
- `compute_dissimilarity()`
- `run_case_study()`
- `module_score()`
- `available_case_studies()`
- `case_study_summary()`
- `case_study_table()`

## Stable outputs

Each case study generates a stable output layout:

- `summary.json`
- `cgm_assignments.tsv`
- `fgm_assignments.tsv`
- `cgm_scores.tsv`
- `fgm_scores.tsv`
- `cgm_enrichment.tsv`
- `fgm_enrichment.tsv`
- `phenotype_associations.tsv`
- `cgm_module_sizes.png`
- `fgm_module_sizes.png`

## Example

```r
library(NestedWGCNA)

set.seed(42)
expr <- matrix(
  abs(rnorm(80 * 500)),
  nrow = 80,
  ncol = 500,
  dimnames = list(paste0("S", seq_len(80)), paste0("G", seq_len(500)))
)

res <- run_nested_wgcna(
  x = expr,
  mode = "paper",
  min_cgm_size = 30,
  min_fgm_size = 10,
  top_n_genes = 400,
  seed = 42
)

summary(res)
```

## Tutorials

The tutorial site is organized around:

- package article (method and API)
- TCGA-BLCA case study
- ALL leukemia case study
- bladderbatch case study

Run the real-data tutorial pipelines:

```r
system("Rscript inst/scripts/run_tcga_blca_case_study.R")
system("Rscript inst/scripts/run_all_leukemia_case_study.R")
system("Rscript inst/scripts/run_bladderbatch_case_study.R")
```

## What NestedWGCNA cannot do yet

- It does not perform raw-count normalization or alignment/quantification
- It does not guarantee exact UMAP/HDBSCAN label identity across environments
- It does not include heavy downstream survival/response modeling workflows
- It does not include controlled-access cohort bundles

## Package layout

- `R/preprocess.R`: input validation and filtering helpers
- `R/adjacency.R`: adjacency computation
- `R/dissimilarity.R`: dissimilarity computation
- `R/clustering.R`: UMAP/HDBSCAN wrappers and CGM discovery
- `R/genfocus.R`: GenFocus normalization
- `R/pipeline.R`: nested workflow and scoring
- `R/case_study_pipeline.R`: real-data case-study runner
- `vignettes/`: package article and real-data tutorials
- `inst/extdata/case_studies/`: bundled case-study outputs
