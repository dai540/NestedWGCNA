# NestedWGCNA

[![R-CMD-check](https://github.com/dai540/NestedWGCNA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/NestedWGCNA/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`NestedWGCNA` is a computation-first R package for two-stage gene co-expression analysis.
This release is intentionally minimal: it focuses on network construction, module discovery,
and nested decomposition only.

## Design policy (v1.00)

- Focus only on computation.
- No visualization helpers in the package API.
- No bundled case-study pipelines or tutorial data.

## Installation

```r
install.packages("remotes")
remotes::install_github("dai540/NestedWGCNA")
```

## Core functions

- `compute_adjacency()`
- `compute_dissimilarity()`
- `run_gene_umap()`
- `run_gene_hdbscan()`
- `core_decompose_weighted()`
- `erode_to_core()`
- `find_cgm()`
- `genfocus_normalize()`
- `find_fgm()`
- `run_nested_wgcna()`
- `module_score()`

## Minimal example

```r
library(NestedWGCNA)

set.seed(42)
x <- matrix(
  abs(rnorm(80 * 400)),
  nrow = 80,
  ncol = 400,
  dimnames = list(paste0("S", 1:80), paste0("G", 1:400))
)

res <- run_nested_wgcna(
  x = x,
  mode = "paper",
  min_cgm_size = 30,
  min_fgm_size = 10,
  top_n_genes = 300,
  seed = 42
)

summary(res)
```

## Reproducibility modes

- `mode = "paper"`: adjacency `r^2`, dissimilarity `sqrt(1 - r^2)`
- `mode = "python_compat"`: adjacency `|r|`, dissimilarity `sqrt(1 - |r|)`
