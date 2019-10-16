context("test utils functions")

test_that("merge_metadata can merge a Seurat object and a dataframe", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata1 = pbmc_obj,
                                    metadata2 = metadata,
                                    log_file = NULL)
  expect_equal(ncol(merged_metadata), 9)
})

test_that("merge_metadata can merge two dataframes", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata1 = metadata,
                                    metadata2 = metadata,
                                    log_file = NULL)
  expect_equal(ncol(merged_metadata), 9)
})

test_that("merge_metadata can merge a dataframe and a tibble", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  metadata <- pbmc_obj@meta.data
  merged_metadata <- merge_metadata(metadata1 = metadata,
                                    metadata2 = as_tibble(metadata),
                                    log_file = NULL)
  expect_equal(ncol(merged_metadata), 9)
})

test_that("as_data_frame_seurat converts metadata to dataframe", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                         package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  s_obj <- add_seurat_assay(seurat_obj = s_obj,
                            assay = "ADT",
                            counts_matrix = counts$`Antibody Capture`,
                            log_file = NULL)
  s_obj_filt <- filter_data(data = s_obj,
                            log_file = NULL,
                            min_genes = NULL,
                            max_genes = NULL,
                            max_mt = 10)
  s_obj_norm <- normalize_data(data = s_obj_filt,
                               method = "sct",
                               assay = "RNA")
  s_obj_dr <- run_dr(s_obj_norm, dr_method = "pca",
                     prefix = "test", var_features = TRUE,
                     num_pcs = 20, assay = "SCT")

  s_obj_dr <- run_dr(s_obj_dr, dr_method = "umap",
                     prefix = "test", reduction = "pcatest",
                     num_dim_use = 20, assay = "SCT", num_neighbors = 6)

  metadata <- as_data_frame_seurat(seurat_obj = s_obj_dr,
                       metadata = TRUE)

  expect_gt(ncol(metadata), ncol(s_obj_dr@meta.data))
})


test_that("as_data_frame_seurat converts reduction to dataframe", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                         package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  s_obj <- add_seurat_assay(seurat_obj = s_obj,
                            assay = "ADT",
                            counts_matrix = counts$`Antibody Capture`,
                            log_file = NULL)
  s_obj_filt <- filter_data(data = s_obj,
                            log_file = NULL,
                            min_genes = NULL,
                            max_genes = NULL,
                            max_mt = 10)
  s_obj_norm <- normalize_data(data = s_obj_filt,
                               method = "sct",
                               assay = "RNA")
  s_obj_dr <- run_dr(s_obj_norm, dr_method = "pca",
                     prefix = "test", var_features = TRUE,
                     num_pcs = 20, assay = "SCT")

  s_obj_dr <- run_dr(s_obj_dr, dr_method = "umap",
                     prefix = "test", reduction = "pcatest",
                     num_dim_use = 20, assay = "SCT", num_neighbors = 6)

  metadata <- as_data_frame_seurat(seurat_obj = s_obj_dr,
                                   metadata = FALSE,
                                   reduction = "pcatest")

  expect_equal(ncol(metadata), 21)
})

test_that("as_data_frame_seurat converts assay to dataframe", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                         package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  s_obj <- add_seurat_assay(seurat_obj = s_obj,
                            assay = "ADT",
                            counts_matrix = counts$`Antibody Capture`,
                            log_file = NULL)
  s_obj_filt <- filter_data(data = s_obj,
                            log_file = NULL,
                            min_genes = NULL,
                            max_genes = NULL,
                            max_mt = 10)
  s_obj_norm <- normalize_data(data = s_obj_filt,
                               method = "sct",
                               assay = "RNA")
  s_obj_dr <- run_dr(s_obj_norm, dr_method = "pca",
                     prefix = "test", var_features = TRUE,
                     num_pcs = 20, assay = "SCT")

  s_obj_dr <- run_dr(s_obj_dr, dr_method = "umap",
                     prefix = "test", reduction = "pcatest",
                     num_dim_use = 20, assay = "SCT", num_neighbors = 6)

  metadata <- as_data_frame_seurat(seurat_obj = s_obj_dr,
                                   metadata = FALSE,
                                   assay = "RNA",
                                   slot = "counts",
                                   features = c("MAP4", "FTL"))

  expect_equal(ncol(metadata), 3)
})

test_that("as_data_frame_seurat converts metadata, reduction, and assay to dataframe", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                         package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  s_obj <- add_seurat_assay(seurat_obj = s_obj,
                            assay = "ADT",
                            counts_matrix = counts$`Antibody Capture`,
                            log_file = NULL)
  s_obj_filt <- filter_data(data = s_obj,
                            log_file = NULL,
                            min_genes = NULL,
                            max_genes = NULL,
                            max_mt = 10)
  s_obj_norm <- normalize_data(data = s_obj_filt,
                               method = "sct",
                               assay = "RNA")
  s_obj_dr <- run_dr(s_obj_norm, dr_method = "pca",
                     prefix = "test", var_features = TRUE,
                     num_pcs = 20, assay = "SCT")

  s_obj_dr <- run_dr(s_obj_dr, dr_method = "umap",
                     prefix = "test", reduction = "pcatest",
                     num_dim_use = 20, assay = "SCT", num_neighbors = 6)

  metadata <- as_data_frame_seurat(seurat_obj = s_obj_dr,
                                   metadata = TRUE,
                                   assay = "RNA",
                                   slot = "counts",
                                   features = c("MAP4", "FTL"),
                                   reduction = "pcatest")

  expect_equal(ncol(metadata), 31)
})
