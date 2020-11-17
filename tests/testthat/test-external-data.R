context("test external counts table")

test_that("test matrix exists", {
  pbmc_mat <- get_test_counts_matrix()
  expect_is(pbmc_mat, "matrix")
})

test_that("test matrix is of reasonable size", {
  pbmc_mat <- get_test_counts_matrix()
  expect_gt(nrow(pbmc_mat), 100)
  expect_gt(ncol(pbmc_mat), 50)
  expect_lt(nrow(pbmc_mat), 10000)
  expect_lt(ncol(pbmc_mat), 10000)
})

test_that("test matrix contains some expected genes", {
  pbmc_mat <- get_test_counts_matrix()
  gene_names <- rownames(pbmc_mat)
  expect_match(gene_names, "CD79A", fixed = TRUE, all = FALSE)
  expect_match(gene_names, "CD79B", fixed = TRUE, all = FALSE)
  expect_match(gene_names, "CD3D", fixed = TRUE, all = FALSE)
  expect_match(gene_names, "LYZ", fixed = TRUE, all = FALSE)
})
