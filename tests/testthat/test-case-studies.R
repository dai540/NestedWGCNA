test_that("case study registry is stable", {
  expect_setequal(
    available_case_studies(),
    c("tcga_blca", "all_leukemia", "bladderbatch")
  )
})

test_that("case study summary is readable", {
  s <- case_study_summary("tcga_blca")
  expect_true(is.list(s))
  expect_equal(s$dataset_id, "tcga_blca")
  expect_true(is.numeric(s$n_samples))
  expect_true(is.numeric(s$cgm_n_modules))
})

test_that("case study table loader works", {
  tab <- case_study_table("bladderbatch", "phenotype_associations")
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("layer", "module", "phenotype", "fdr") %in% colnames(tab)))
})
