test_that("pca can be calcuated", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  s_obj <- add_seurat_assay(
    x = s_obj,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  )
  s_obj_filt <- filter_cells(
    x = s_obj,
    log_file = NULL,
    min_counts = 1,
    min_genes = 1,
    max_genes = NULL,
    max_mt = 10
  )
  pc <- run_pca(
    x = t(GetAssayData(s_obj_filt, assay = "RNA", layer = "counts")),
    num_pcs = 20
  )
  expect_equal(nrow(pc$feature.loadings), nrow(s_obj_filt))
})

test_that("tsne can be calcuated", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  s_obj <- add_seurat_assay(
    x = s_obj,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  )
  s_obj_filt <- filter_cells(
    x = s_obj,
    log_file = NULL,
    min_counts = 1,
    min_genes = 1,
    max_genes = NULL,
    max_mt = 10
  )
  pc <- run_pca(
    x = t(GetAssayData(s_obj_filt, assay = "RNA", layer = "counts")),
    num_pcs = 20
  )
  tsne <- run_tsne(
    x = pc$cell.embeddings,
    seed.use = 1,
    dim.embed = 2
  )
  expect_equal(nrow(tsne), ncol(s_obj_filt))
})

test_that("umap can be calcuated", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  s_obj <- add_seurat_assay(
    x = s_obj,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  )
  s_obj_filt <- filter_cells(
    x = s_obj,
    log_file = NULL,
    min_counts = 1,
    min_genes = 1,
    max_genes = NULL,
    max_mt = 10
  )
  pc <- run_pca(
    x = t(GetAssayData(s_obj_filt, assay = "RNA", layer = "counts")),
    num_pcs = 20
  )
  umap <- run_umap(
    x = pc$cell.embeddings,
    num_neighbors = 10
  )
  expect_equal(nrow(umap), ncol(s_obj_filt))
})

test_that("run dr", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  s_obj <- add_seurat_assay(
    x = s_obj,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  )
  s_obj_filt <- filter_cells(
    x = s_obj,
    log_file = NULL,
    min_counts = 1,
    min_genes = 1,
    max_genes = NULL,
    max_mt = 10
  )
  s_obj_norm <- normalize_counts(
    x = s_obj_filt,
    method = "sct",
    assay = "RNA"
  )
  s_obj_dr <- run_dr(
    s_obj_norm,
    dr_method = "pca",
    suffix = "test",
    var_features = TRUE,
    num_pcs = 20,
    assay = "SCT"
  )
  expect_equal(names(s_obj_dr@reductions), "pcatest")
})

test_that("run dr umap", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  s_obj <- add_seurat_assay(
    x = s_obj,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  )
  s_obj_filt <- filter_cells(
    x = s_obj,
    log_file = NULL,
    min_counts = 1,
    min_genes = 1,
    max_genes = NULL,
    max_mt = 10
  )
  s_obj_norm <- normalize_counts(
    x = s_obj_filt,
    method = "sct",
    assay = "RNA"
  )
  s_obj_dr <- run_dr(
    s_obj_norm,
    dr_method = "pca",
    suffix = "test",
    var_features = TRUE,
    num_pcs = 20,
    assay = "SCT"
  )

  s_obj_dr <- run_dr(
    s_obj_dr,
    dr_method = "umap",
    suffix = "test",
    reduction = "pcatest",
    num_dim = 20,
    assay = "SCT",
    num_neighbors = 6
  )

  expect_equal(names(s_obj_dr@reductions), c("pcatest", "umaptest"))
})
