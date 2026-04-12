# Compute adjacency matrix

Compute adjacency matrix

## Usage

``` r
compute_adjacency(
  x,
  method = c("pearson", "spearman"),
  mode = c("paper", "python_compat", "custom"),
  signed = FALSE,
  power = NULL
)
```

## Arguments

- x:

  Expression matrix (`samples x genes`).

- method:

  Correlation method (`"pearson"` or `"spearman"`).

- mode:

  Similarity mode: `"paper"`, `"python_compat"`, or `"custom"`.

- signed:

  Whether to use signed adjacency.

- power:

  Optional power parameter. If `NULL`, defaults by mode.

## Value

Gene-gene adjacency matrix.
