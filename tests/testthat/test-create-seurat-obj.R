context("test create seurat obj")

test_that("seurat obj can be created", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  expect_s4_class(pbmc_obj, "Seurat")
})

test_that("Seurat object can be filtered", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat)
  pbmc_obj <- filter_data(pbmc_obj, log_file = NULL,
                          min_genes = NULL, max_genes = NULL, max_mt = 10)
  expect_lt(ncol(pbmc_obj), 60)
})

test_that("Read in 10x", {
  counts <- import_matrix(data_path = system.file("extdata",
                                                  "outs/filtered_feature_bc_matrix",
                                                  package = "scooter"),
                          gene.column = 2)
  truth <- readRDS(system.file("extdata",
                               "import_matrix_counts.rds",
                               package = "scooter"))
  expect_identical(counts, truth )
})

test_that("Read in csv", {
  counts <- load_sample_counts_matrix(sample_names = "test",
                                      AC_text_file = system.file("extdata",
                                                  "HTO.csv",
                                                  package = "scooter"),
                                      delim = ",")

})

test_that("Read in tsv", {
  counts <- load_sample_counts_matrix(sample_names = "test",
                                      AC_text_file = system.file("extdata",
                                                                 "HTO.tsv",
                                                                 package = "scooter"),
                                      delim = "\t")

})

test_that("Read in 10x", {
  counts <- load_sample_counts_matrix(sample_names = "test",
                                      data_path_10x = system.file("extdata",
                                                                  "",
                                                                  package = "scooter"))

})

# test_that("multi-sample seurat obj can be created", {
#   pbmc_mat <- get_test_counts_matrix()
#   pbmc_obj <- create_seurat_obj(counts_matrix = pbmc_mat, min_genes = 10, out_dir = "x")
#   expect_s4_class(pbmc_obj, "seurat")
# })
