# NestedWGCNA

`NestedWGCNA` is a minimal R package for **two-stage gene co-expression analysis**.
It is intentionally rebuilt from scratch with a strict small-footprint policy:
only computation-focused code is kept, heavy bundled data is removed, and the
directory layout is minimized.

## Package Scope

This package focuses only on the core computational path:

1. Build a gene-gene adjacency matrix from a sample-by-gene matrix.
2. Convert adjacency to a dissimilarity matrix.
3. Detect **coarse-grained modules (CGM)**.
4. Select module core genes by within-module connectivity.
5. Normalize expression by module core baselines.
6. Detect **fine-grained modules (FGM)** on normalized expression.
7. Calculate module scores.

Excluded by design:

- Large demo or real datasets.
- Plot-heavy helper APIs.
- Domain-specific downstream analyses (enrichment, survival, etc.).

## Design Principles

- Minimal dependencies (base R + `stats`).
- Deterministic and lightweight defaults.
- Clear input contract: matrix shape is `samples x genes`.
- Explicit reproducibility modes:
  - `mode = "paper"` uses `r^2` adjacency.
  - `mode = "python_compat"` uses `|r|` adjacency.
- Small repository size:
  - No bundled bulky files.
  - No downloaded external datasets.
  - No committed build artifacts (`*.tar.gz`, `*.Rcheck`, `tmp_*`, generated `docs/`).

## Installation

```r
install.packages("remotes")
remotes::install_github("dai540/NestedWGCNA")
```

## Input Requirements

- Numeric matrix or data frame.
- Rows are samples (or cells), columns are genes.
- At least 2 samples and 2 genes.
- Constant or non-finite genes are removed by `clear_data()`.

## Main API

- `as_expression_matrix()`
- `clear_data()`
- `top_variable_genes()`
- `compute_adjacency()`
- `compute_dissimilarity()`
- `find_cgm()`
- `find_fgm()`
- `module_score()`
- `run_nested_wgcna()`

## Quick Start

```r
library(NestedWGCNA)

set.seed(1)
x <- matrix(
  rnorm(80 * 600),
  nrow = 80,
  ncol = 600,
  dimnames = list(paste0("S", seq_len(80)), paste0("G", seq_len(600)))
)

res <- run_nested_wgcna(
  x = x,
  mode = "paper",
  min_cgm_size = 60,
  min_fgm_size = 20,
  top_n_genes = 400
)

summary(res)
head(res$fgm$assignment)
```

## Output Object

`run_nested_wgcna()` returns an object of class `NestedWGCNAResult` containing:

- `input`: cleaned matrix dimensions.
- `params`: run parameters.
- `cgm`: CGM assignments and matrices.
- `cgm_core`: selected core genes per CGM.
- `normalized_matrix`: core-normalized matrix.
- `fgm`: FGM assignments and matrices.
- `module_scores`: sample-level module score matrix.

## pkgdown Documentation Structure

This package includes pkgdown configuration with four article groups:

- Getting Started
- Guides
- Tutorials
- Reference

Each group has at least one article in `vignettes/`, and the function reference
is generated from roxygen documentation.

## Reproducibility Notes

- Clustering is performed with base `hclust` to keep dependency and runtime cost low.
- Module labels are deterministic for fixed input and parameters.
- This package does not download data automatically.

## Minimal Directory Policy

The repository is intentionally small and contains only files needed for:

- package build/check,
- test execution,
- pkgdown configuration,
- concise vignettes based on simulated data.

Temporary files and build artifacts are excluded by `.gitignore`.
