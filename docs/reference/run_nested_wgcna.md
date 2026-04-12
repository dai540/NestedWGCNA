# Run NestedWGCNA end-to-end

Run NestedWGCNA end-to-end

## Usage

``` r
run_nested_wgcna(
  x,
  mode = c("paper", "python_compat"),
  min_cgm_size = 100,
  min_fgm_size = 10,
  corr_method = c("pearson", "spearman"),
  genfocus_corr_method = c("spearman", "pearson"),
  target_cgm = NULL,
  top_n_genes = NULL,
  seed = 42
)
```

## Arguments

- x:

  Expression matrix (`samples x genes`).

- mode:

  Similarity mode.

- min_cgm_size:

  Minimum CGM cluster size.

- min_fgm_size:

  Minimum FGM cluster size.

- corr_method:

  Correlation method for CGM stage.

- genfocus_corr_method:

  Correlation method for GenFocus.

- target_cgm:

  Optional target CGM id for FGM stage. If `NULL`, uses largest
  non-noise CGM.

- top_n_genes:

  Optional variable-gene filter size.

- seed:

  Random seed.

## Value

Object of class `NestedWGCNAResult`.
