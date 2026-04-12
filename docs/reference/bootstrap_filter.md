# Bootstrap-style filter (Python `bootstrap_filter` equivalent)

Keeps genes with classification error larger than `threshold`.

## Usage

``` r
bootstrap_filter(x, threshold = 1/exp(1))
```

## Arguments

- x:

  Expression matrix (`samples x genes`).

- threshold:

  Threshold (default `1 / exp(1)`).

## Value

Filtered matrix.
