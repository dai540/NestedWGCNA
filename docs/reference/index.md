# Package index

## Package

- [`NestedWGCNA-package`](https://dai540.github.io/NestedWGCNA/reference/NestedWGCNA-package.md)
  [`NestedWGCNA`](https://dai540.github.io/NestedWGCNA/reference/NestedWGCNA-package.md)
  : NestedWGCNA

## End-to-end pipeline

- [`run_nested_wgcna()`](https://dai540.github.io/NestedWGCNA/reference/run_nested_wgcna.md)
  : Run NestedWGCNA end-to-end
- [`find_cgm()`](https://dai540.github.io/NestedWGCNA/reference/find_cgm.md)
  : Discover coarse-grained modules (CGMs)
- [`find_fgm()`](https://dai540.github.io/NestedWGCNA/reference/find_fgm.md)
  : Discover FGMs inside a CGM
- [`genfocus_normalize()`](https://dai540.github.io/NestedWGCNA/reference/genfocus_normalize.md)
  : GenFocus normalization (R rewrite of Python reference logic)
- [`module_score()`](https://dai540.github.io/NestedWGCNA/reference/module_score.md)
  : Calculate per-sample module scores

## Network primitives

- [`compute_adjacency()`](https://dai540.github.io/NestedWGCNA/reference/compute_adjacency.md)
  : Compute adjacency matrix
- [`compute_dissimilarity()`](https://dai540.github.io/NestedWGCNA/reference/compute_dissimilarity.md)
  : Compute dissimilarity matrix
- [`core_decompose_weighted()`](https://dai540.github.io/NestedWGCNA/reference/core_decompose_weighted.md)
  : Weighted core decomposition
- [`erode_to_core()`](https://dai540.github.io/NestedWGCNA/reference/erode_to_core.md)
  : Erode clusters to core genes

## Clustering wrappers

- [`run_gene_umap()`](https://dai540.github.io/NestedWGCNA/reference/run_gene_umap.md)
  : Run UMAP on genes using precomputed dissimilarity
- [`run_gene_hdbscan()`](https://dai540.github.io/NestedWGCNA/reference/run_gene_hdbscan.md)
  : Run HDBSCAN on gene embedding

## Input preprocessing

- [`as_expression_matrix()`](https://dai540.github.io/NestedWGCNA/reference/as_expression_matrix.md)
  : Coerce input to numeric expression matrix

- [`clear_data()`](https://dai540.github.io/NestedWGCNA/reference/clear_data.md)
  :

  Basic cleaning (Python `clear_data` equivalent)

- [`bootstrap_filter()`](https://dai540.github.io/NestedWGCNA/reference/bootstrap_filter.md)
  :

  Bootstrap-style filter (Python `bootstrap_filter` equivalent)

- [`classification_err()`](https://dai540.github.io/NestedWGCNA/reference/classification_err.md)
  : Classification error of a vector

- [`top_variable_genes()`](https://dai540.github.io/NestedWGCNA/reference/top_variable_genes.md)
  : Select top variable genes

## Result object

- [`summary(`*`<NestedWGCNAResult>`*`)`](https://dai540.github.io/NestedWGCNA/reference/summary.NestedWGCNAResult.md)
  : Summarize a NestedWGCNA result
