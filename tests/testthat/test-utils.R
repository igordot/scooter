context("test utils functions")

test_that("merge_metadata can merge a Seurat object and a dataframe", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat,
                                min_genes = 10, out_dir = "x")
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata_one = pbmc_obj,
                                    metadata_two = metadata,
                                    log_file = NULL,
                                    write = FALSE,
                                    proj_name = "",
                                    label = "",
                                    out_dir = ".")
  expect_equal(ncol(merged_metadata), 9)
})
