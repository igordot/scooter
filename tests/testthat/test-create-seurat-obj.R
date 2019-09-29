context("test create seurat obj")

test_that("seurat obj can be created", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  expect_s4_class(pbmc_obj, "Seurat")
})

test_that("Seurat object can be filtered", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  pbmc_obj <- filter_data(pbmc_obj, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10)
  expect_lt(ncol(pbmc_obj), 60)
})


# test_that("multi-sample seurat obj can be created", {
#   pbmc_mat <- get_test_counts_matrix()
#   pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat, min_genes = 10, out_dir = "x")
#   expect_s4_class(pbmc_obj, "seurat")
# })
