# NestedWGCNA

`NestedWGCNA` is a two-stage gene co-expression workflow for bulk transcriptome data.

- Stage 1: discover coarse-grained modules (CGMs).
- Stage 2: normalize within CGM cores (GenFocus) and discover fine-grained modules (FGMs).

This repository is organized as an R package project with reproducible, real-data tutorial pipelines.

## What is included

- R package scaffold (`DESCRIPTION`, `R/`, `tests/`, `vignettes/`).
- Python reference code (`inst/python_ref/utils.py`, `inst/python_ref/GenFocus.py`).
- Shared case-study engine (`inst/scripts/run_matrix_case_study.py`).
- Multiple executed case studies in `inst/extdata/case_studies/`.

## Real-data tutorials (multiple datasets)

All tutorials are based on executed analyses and include background, objective, results, discussion, and interpretation.

1. TCGA-BLCA RNA-seq (downloaded from upstream NestedWGCNA data URL)
```bash
python inst/scripts/run_tcga_blca_case_study.py
```
2. ALL leukemia microarray (downloaded as Bioconductor package data)
```bash
Rscript inst/scripts/run_all_leukemia_case_study.R
```
3. bladderbatch bladder cancer microarray (downloaded as Bioconductor package data)
```bash
Rscript inst/scripts/run_bladderbatch_case_study.R
```

Outputs are written to:

- `inst/extdata/case_studies/tcga_blca/`
- `inst/extdata/case_studies/all_leukemia/`
- `inst/extdata/case_studies/bladderbatch/`

Each directory contains `summary.json`, module assignments, module scores, phenotype associations, enrichment tables, and figures.

## Documentation / pkgdown

Main articles in `vignettes/`:

- `package-article.Rmd`
- `tutorial-tcga-blca.Rmd`
- `tutorial-all-leukemia.Rmd`
- `tutorial-bladderbatch.Rmd`

Build site:

```r
pkgdown::build_site()
```

## Notes

- Runtime compatibility patches are applied for modern `scikit-learn` + UMAP.
- GenFocus set operations are patched for deterministic outputs.
- `mode = "paper"` vs `mode = "python_compat"` policy is documented in the package article and method spec.

## Citation

Dyugay et al. *Improved gene co-expression network analysis and its application for biomarker discovery of the response to immunotherapy.*
