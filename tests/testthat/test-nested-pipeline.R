test_that("run_nested_wgcna returns expected structure", {
  skip_if_not_installed("uwot")
  skip_if_not_installed("dbscan")

  set.seed(7)
  x <- matrix(
    abs(rnorm(80 * 220)),
    nrow = 80,
    ncol = 220,
    dimnames = list(paste0("S", seq_len(80)), paste0("G", seq_len(220)))
  )

  res <- run_nested_wgcna(
    x = x,
    mode = "paper",
    min_cgm_size = 20,
    min_fgm_size = 8,
    top_n_genes = 200,
    seed = 42
  )

  expect_s3_class(res, "NestedWGCNAResult")
  expect_true(all(c("cgm", "fgm") %in% names(res$module_scores)))
  expect_equal(length(res$cgm$assignment), 200)
  expect_equal(length(res$fgm$assignment), 200)
})
