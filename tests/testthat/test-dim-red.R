test_that("pca can be calcuated", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      data_path_10x = system.file("extdata", "",
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
                                      data_path_10x = system.file("extdata", "",
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
                                      data_path_10x = system.file("extdata", "",
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


