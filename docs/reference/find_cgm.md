# Discover coarse-grained modules (CGMs)

Discover coarse-grained modules (CGMs)

## Usage

``` r
find_cgm(
  x,
  mode = c("paper", "python_compat"),
  corr_method = c("pearson", "spearman"),
  min_cluster_size = 30,
  n_components = NULL,
  seed = 42
)
```

## Arguments

- x:

  Expression matrix (`samples x genes`).

- mode:

  Similarity mode.

- corr_method:

  Correlation method.

- min_cluster_size:

  Minimum cluster size.

- n_components:

  UMAP dimensions; if `NULL`, computed as in Python code.

- seed:

  Random seed.

## Value

List with CGM assignments, core assignments, embedding, and matrices.
