# Calculate per-sample module scores

Calculate per-sample module scores

## Usage

``` r
module_score(x, assignment, prefix = "M_")
```

## Arguments

- x:

  Expression matrix (`samples x genes`).

- assignment:

  Integer module assignment vector for genes (`-1` = noise).

- prefix:

  Prefix for score column names.

## Value

Data frame (`samples x modules`) of module scores.
