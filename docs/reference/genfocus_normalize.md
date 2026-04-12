# GenFocus normalization (R rewrite of Python reference logic)

GenFocus normalization (R rewrite of Python reference logic)

## Usage

``` r
genfocus_normalize(
  x,
  focus = c("eigengene", "gene"),
  focus_gene = NULL,
  corr_method = c("spearman", "pearson"),
  input_scale = c("tpm", "fpkm"),
  corr_thr = 0.9,
  CVR_thr = 0.6
)
```

## Arguments

- x:

  Expression matrix (`samples x genes`).

- focus:

  `"eigengene"` or `"gene"`.

- focus_gene:

  Focus gene name when `focus = "gene"`.

- corr_method:

  Correlation method.

- input_scale:

  `"tpm"` or `"fpkm"`.

- corr_thr:

  Correlation threshold for INGS selection.

- CVR_thr:

  CVR threshold for clipping.

## Value

`NULL` if INGS cannot be established; otherwise a list with normalized
matrix and diagnostics.
