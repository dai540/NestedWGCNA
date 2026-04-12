# Run UMAP on genes using precomputed dissimilarity

Run UMAP on genes using precomputed dissimilarity

## Usage

``` r
run_gene_umap(dissimilarity, n_components = 10, n_neighbors = 30, seed = 42)
```

## Arguments

- dissimilarity:

  Gene-gene dissimilarity matrix.

- n_components:

  UMAP embedding dimensions.

- n_neighbors:

  UMAP neighbors.

- seed:

  Random seed.

## Value

Numeric embedding matrix.
