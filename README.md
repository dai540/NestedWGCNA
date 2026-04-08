# NestedWGCNA

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/NestedWGCNA/)
[![R-CMD-check](https://github.com/dai540/NestedWGCNA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/NestedWGCNA/actions/workflows/R-CMD-check.yaml)
[![pkgdown deploy](https://github.com/dai540/NestedWGCNA/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/dai540/NestedWGCNA/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`NestedWGCNA` is an R package for two-stage gene co-expression analysis:

- Stage 1: discover coarse-grained modules (CGMs)
- Stage 2: normalize by CGM core genes (GenFocus) and discover fine-grained modules (FGMs)

<https://dai540.github.io/NestedWGCNA/>

This package is a full R implementation. The prior Python implementation is not used at runtime.

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

Then load:

```r
library(NestedWGCNA)
```

Core clustering backends:

```r
install.packages(c("uwot", "dbscan", "matrixStats"))
```

## What NestedWGCNA does

`NestedWGCNA` standardizes:

- adjacency and dissimilarity calculation with explicit mode control
- CGM discovery with UMAP + HDBSCAN
- core decomposition and core erosion
- GenFocus normalization
- nested CGM -> FGM workflow
- module score calculation
- case-study runner with enrichment and phenotype association outputs

The package provides two main method modes:

- `mode = "paper"`: adjacency \(r^2\), dissimilarity \(\sqrt{1-r^2}\)
- `mode = "python_compat"`: adjacency \(|r|\), dissimilarity \(\sqrt{1-|r|}\)

## Main functions

- `run_nested_wgcna()`
- `find_cgm()`
- `find_fgm()`
- `genfocus_normalize()`
- `compute_adjacency()`
- `compute_dissimilarity()`
- `run_case_study()`
- `available_case_studies()`
- `case_study_summary()`
- `case_study_table()`

## Quick example

```r
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

names(res$module_scores)
```

## Real-data tutorials (multiple datasets)

Each tutorial uses downloaded real data and includes background, objective, methods, results, discussion, and interpretation.

1. TCGA-BLCA bulk RNA-seq:
```r
system("Rscript inst/scripts/run_tcga_blca_case_study.R")
```
2. Bioconductor ALL leukemia:
```r
system("Rscript inst/scripts/run_all_leukemia_case_study.R")
```
3. Bioconductor bladderbatch:
```r
system("Rscript inst/scripts/run_bladderbatch_case_study.R")
```

Outputs are written under:

- `inst/extdata/case_studies/tcga_blca/`
- `inst/extdata/case_studies/all_leukemia/`
- `inst/extdata/case_studies/bladderbatch/`

Each case-study directory contains:

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

## Documentation

Main pkgdown articles:

- package article
- tutorial: TCGA-BLCA
- tutorial: ALL leukemia
- tutorial: bladderbatch

Build docs locally:

```r
pkgdown::build_site()
```

## Author

- Dai

## Citation

If you use `NestedWGCNA`, cite:

> Dai (2026). *NestedWGCNA: Two-Stage Gene Co-Expression Network Analysis*. R package.
> <https://dai540.github.io/NestedWGCNA/>
