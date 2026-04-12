# Run HDBSCAN on gene embedding

Run HDBSCAN on gene embedding

## Usage

``` r
run_gene_hdbscan(embedding, min_cluster_size = 30, min_samples = NULL)
```

## Arguments

- embedding:

  UMAP embedding.

- min_cluster_size:

  Minimum cluster size.

- min_samples:

  Ignored in current implementation; kept for API parity.

## Value

Integer labels (`-1` for noise, `0..k-1` for clusters).
