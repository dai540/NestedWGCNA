test_that("plot_module_network returns a correlation matrix", {
  set.seed(1)
  scores <- data.frame(
    CGM_0 = rnorm(40),
    CGM_1 = rnorm(40),
    CGM_2 = rnorm(40),
    CGM_3 = rnorm(40)
  )

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file, width = 600, height = 500)
  on.exit({
    grDevices::dev.off()
    unlink(png_file)
  }, add = TRUE)

  r <- plot_module_network(scores, layer = "CGM", top_n = 4, cor_thr = 0.2)
  expect_true(is.matrix(r))
  expect_equal(dim(r), c(4, 4))
})

test_that("plot_module_score_heatmap returns a reordered matrix", {
  set.seed(2)
  scores <- data.frame(
    FGM_0 = rnorm(35),
    FGM_1 = rnorm(35),
    FGM_2 = rnorm(35),
    FGM_3 = rnorm(35),
    FGM_4 = rnorm(35)
  )

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file, width = 700, height = 600)
  on.exit({
    grDevices::dev.off()
    unlink(png_file)
  }, add = TRUE)

  r <- plot_module_score_heatmap(scores, layer = "FGM", top_n = 5)
  expect_true(is.matrix(r))
  expect_equal(nrow(r), ncol(r))
  expect_equal(nrow(r), 5)
})
