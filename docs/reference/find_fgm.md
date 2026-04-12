# Discover FGMs inside a CGM

Discover FGMs inside a CGM

## Usage

``` r
find_fgm(
  x,
  mode = c("paper", "python_compat"),
  min_cluster_size = 10,
  genfocus_corr_method = c("spearman", "pearson"),
  genfocus_corr_thr = 0.9,
  genfocus_CVR_thr = 0.6,
  seed = 42
)
```

## Arguments

- x:

  Expression matrix (`samples x genes`) for a single target CGM.

- mode:

  Similarity mode.

- min_cluster_size:

  Minimum FGM cluster size.

- genfocus_corr_method:

  Correlation method for GenFocus.

- genfocus_corr_thr:

  Correlation threshold for INGS.

- genfocus_CVR_thr:

  CVR clipping threshold.

- seed:

  Random seed.

## Value

List with FGM assignments and normalization output.
