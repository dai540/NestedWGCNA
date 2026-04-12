# Compute dissimilarity matrix

Compute dissimilarity matrix

## Usage

``` r
compute_dissimilarity(
  adjacency,
  method = c("paper_metric", "python_compat", "tom", "mot", "custom"),
  custom_fun = NULL
)
```

## Arguments

- adjacency:

  Adjacency matrix.

- method:

  One of `"paper_metric"`, `"python_compat"`, `"tom"`, `"mot"`, or
  `"custom"`.

- custom_fun:

  Optional function for custom dissimilarity.

## Value

Gene-gene dissimilarity matrix.
