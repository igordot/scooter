context("test utils functions")

test_that("merge_metadata can merge a Seurat object and a dataframe", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat,
                                out_dir = "x")
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata1 = pbmc_obj,
                                    metadata2 = metadata,
                                    log_file = NULL,
                                    write = FALSE)
  expect_equal(ncol(merged_metadata), 9)
})

test_that("merge_metadata can merge two dataframes", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat,
                                out_dir = "x")
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata1 = metadata,
                                    metadata2 = metadata,
                                    log_file = NULL,
                                    write = FALSE)
  expect_equal(ncol(merged_metadata), 9)
})

test_that("merge_metadata can merge a dataframe and a tibble", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat,
                                out_dir = "x")
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata1 = metadata,
                                    metadata2 = as_tibble(metadata),
                                    log_file = NULL,
                                    write = FALSE)
  expect_equal(ncol(merged_metadata), 9)
})

