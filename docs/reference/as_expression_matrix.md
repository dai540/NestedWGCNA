# Coerce input to numeric expression matrix

Coerce input to numeric expression matrix

## Usage

``` r
as_expression_matrix(x, require_names = TRUE)
```

## Arguments

- x:

  Matrix-like object (`matrix` or `data.frame`) with rows as samples and
  columns as genes.

- require_names:

  If `TRUE`, requires row and column names.

## Value

Numeric matrix.
