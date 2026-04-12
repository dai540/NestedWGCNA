test_that("adjacency and dissimilarity are valid", {
  set.seed(1)
  x <- matrix(rnorm(20 * 100), nrow = 20, ncol = 100)
  colnames(x) <- paste0("g", seq_len(ncol(x)))

  a <- compute_adjacency(x, mode = "paper")
  d <- compute_dissimilarity(a)

  expect_equal(dim(a), c(100, 100))
  expect_equal(dim(d), c(100, 100))
  expect_equal(unname(diag(a)), rep(1, 100))
  expect_equal(unname(diag(d)), rep(0, 100))
})

test_that("nested workflow returns result class", {
  set.seed(2)
  x <- matrix(rnorm(40 * 300), nrow = 40, ncol = 300)
  colnames(x) <- paste0("g", seq_len(ncol(x)))

  res <- run_nested_wgcna(
    x,
    mode = "python_compat",
    min_cgm_size = 30,
    min_fgm_size = 10,
    top_n_genes = 200
  )

  expect_s3_class(res, "NestedWGCNAResult")
  expect_true(is.matrix(res$module_scores))
})
