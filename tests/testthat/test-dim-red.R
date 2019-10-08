test_that("pca can be calcuated", {
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
  pc <- run_pca(data = t(s_obj_filt@assays$RNA@counts),
                num_pcs = 20,
                prefix = "PC_")
  expect_equal(nrow(pc$feature.loadings), nrow(s_obj_filt))
})

test_that("tsne can be calcuated", {
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
  pc <- run_pca(data = t(s_obj_filt@assays$RNA@counts),
                num_pcs = 20,
                prefix = "PC_")
  tsne <- run_tsne(data = pc$cell.embeddings,
                   seed.use = 1,
                   dim.embed = 2,
                   prefix = "tSNE_")
  expect_equal(nrow(tsne), ncol(s_obj_filt))
})

test_that("umap can be calcuated", {
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
  pc <- run_pca(data = t(s_obj_filt@assays$RNA@counts),
                num_pcs = 20,
                prefix = "PC_")
  umap <- run_umap(data = pc$cell.embeddings,
                   num_neighbors = 10,
                   prefix = "UMAP_")
  expect_equal(nrow(umap), ncol(s_obj_filt))
})

test_that("run dr", {
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
  expect_equal(names(s_obj_dr@reductions), "pcatest")
})

test_that("run dr umap", {
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
                     num_dim_use = 20, assay = "SCT", num_neighbor = 10)

  expect_equal(names(s_obj_dr@reductions), c("pcatest", "umaptest"))
})


