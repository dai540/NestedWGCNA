test_that("adjacency differs by mode", {
  set.seed(1)
  x <- matrix(
    rnorm(40 * 30),
    nrow = 40,
    ncol = 30,
    dimnames = list(paste0("S", seq_len(40)), paste0("G", seq_len(30)))
  )

  a_paper <- compute_adjacency(x, mode = "paper")
  a_py <- compute_adjacency(x, mode = "python_compat")

  expect_equal(dim(a_paper), c(30, 30))
  expect_equal(dim(a_py), c(30, 30))
  expect_true(all(diag(a_paper) == 1))
  expect_true(all(diag(a_py) == 1))
  expect_true(any(abs(a_paper - a_py) > 1e-8))
})

test_that("dissimilarity matrix basic properties", {
  set.seed(2)
  x <- matrix(
    rnorm(50 * 25),
    nrow = 50,
    ncol = 25,
    dimnames = list(paste0("S", seq_len(50)), paste0("G", seq_len(25)))
  )
  a <- compute_adjacency(x, mode = "paper")
  d <- compute_dissimilarity(a, method = "paper_metric")

  expect_equal(dim(d), c(25, 25))
  expect_true(max(abs(d - t(d))) < 1e-10)
  expect_true(all(diag(d) == 0))
  expect_true(all(d >= 0))
})
