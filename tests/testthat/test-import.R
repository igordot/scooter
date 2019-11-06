test_that("import_mtx can read in 10x data", {
  counts <- import_mtx(data_path = system.file("extdata",
                                                  "outs/filtered_feature_bc_matrix",
                                                  package = "scooter"),
                          gene_column = 2)
  truth <- readRDS(system.file("extdata",
                               "import_matrix_counts.rds",
                               package = "scooter"))
  expect_identical(counts, truth)
})

test_that("load_sample_counts_matrix can read in 10x data", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata",
                                                         "",
                                                         package = "scooter"))
  expect_identical(names(counts), c("Antibody Capture", "Gene Expression"))
})

test_that("load_sample_counts_matrix can read in tsv Antibody Capture file", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "HTO.tsv",
                                                                 package = "scooter"))
  expect_identical(names(counts), c("Antibody Capture"))
})

test_that("load_sample_counts_matrix can read in csv Antibody Capture file", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata",
                                                                 "HTO.csv",
                                                                 package = "scooter"))
  expect_identical(names(counts), c("Antibody Capture"))
})

test_that("Seurat object can be created from RNA data", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                                  package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  # only 462 genes are above zero in the test data
  expect_s4_class(s_obj, "Seurat")
})

test_that("Seurat object can be created from ADT data", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                                  package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Antibody Capture`,
                             assay = "ADT", log_file = NULL)
  expect_s4_class(s_obj, "Seurat")
})

test_that("Seurat obj can be created", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  expect_s4_class(pbmc_obj, "Seurat")
})

test_that("Assay can be added to Seurat object",{
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                                  package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  s_obj <- add_seurat_assay(seurat_obj = s_obj,
                            assay = "ADT",
                            counts_matrix = counts$`Antibody Capture`,
                            log_file = NULL)
  expect_s4_class(s_obj, "Seurat")
})

test_that("Seurat object can be filtered", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  pbmc_obj <- filter_data(pbmc_obj, log_file = NULL,
                          min_genes = NULL, max_genes = NULL, max_mt = 10)
  expect_lt(ncol(pbmc_obj), 60)
})

test_that("Seurat object can be filtered", {
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
  expect_lt(ncol(s_obj_filt), ncol(s_obj))
})

test_that("Seurat object can be log normalized", {
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
                               method = "log_norm",
                               assay = "RNA")
  expect_gt(ncol(s_obj_norm@assays$RNA@scale.data), 180)
})

test_that("Seurat object can be SC transform", {
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
  expect_gt(ncol(s_obj_norm@assays$SCT@scale.data), 180)
})

