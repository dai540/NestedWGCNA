# Basic cleaning (Python `clear_data` equivalent)

Basic cleaning (Python `clear_data` equivalent)

## Usage

``` r
clear_data(x, fillna = NULL)
```

## Arguments

- x:

  Numeric expression matrix (`samples x genes`).

- fillna:

  If `NULL`, drop genes containing `NA`; otherwise fill `NA` with this
  value.

## Value

Cleaned matrix.
