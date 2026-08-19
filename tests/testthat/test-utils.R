context("test utils functions")

test_that("merge_metadata can merge a Seurat object and a dataframe", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)
  extra_metadata <- data.frame(
    cell = rownames(pbmc_obj@meta.data),
    extra_col = seq_len(nrow(pbmc_obj@meta.data))
  )
  merged_metadata <- merge_metadata(
    x = pbmc_obj,
    y = extra_metadata,
    log_file = NULL
  )
  # counts track the width of initialize_seurat_object()'s metadata (9 columns) plus the one new
  # column from extra_metadata ("cell" itself becomes rownames, not a column, on the Seurat path)
  expect_equal(ncol(merged_metadata), 10)
})

test_that("merge_metadata can merge two dataframes", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)
  metadata <- pbmc_obj@meta.data
  extra_metadata <- data.frame(
    cell = rownames(metadata),
    extra_col = seq_len(nrow(metadata))
  )
  merged_metadata <- merge_metadata(
    x = metadata,
    y = extra_metadata,
    log_file = NULL
  )
  # cell + the 9 original columns + the one new column from extra_metadata
  expect_equal(ncol(merged_metadata), 11)
})

test_that("merge_metadata can merge a dataframe and a tibble", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)
  metadata <- pbmc_obj@meta.data
  extra_metadata <- tibble(
    cell = rownames(metadata),
    extra_col = seq_len(nrow(metadata))
  )
  merged_metadata <- merge_metadata(
    x = metadata,
    y = extra_metadata,
    log_file = NULL
  )
  expect_equal(ncol(merged_metadata), 11)
})

test_that("merge_metadata errors when x and y share column names besides cell", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)
  metadata <- pbmc_obj@meta.data
  expect_error(
    merge_metadata(x = metadata, y = metadata, log_file = NULL),
    "share column names"
  )
})

test_that("merge_metadata warns when x and y do not fully overlap on cell", {
  x <- data.frame(cell = c("a", "b", "c"), val_x = 1:3)
  y <- data.frame(cell = c("b", "c", "d"), val_y = 1:3)
  expect_warning(
    merge_metadata(x = x, y = y, log_file = NULL),
    "cells do not fully overlap"
  )
})

test_that("as_data_frame_seurat converts metadata to dataframe", {
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
    min_genes = 1,
    min_counts = 1,
    max_genes = NULL,
    max_mt = 10
  )
  s_obj_norm <- normalize_counts(x = s_obj_filt, method = "sct", assay = "RNA")
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

  metadata <- as_data_frame_seurat(x = s_obj_dr, metadata = TRUE)

  expect_gt(ncol(metadata), ncol(s_obj_dr@meta.data))
})


test_that("as_data_frame_seurat converts reduction to dataframe", {
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
    min_genes = 1,
    min_counts = 1,
    max_genes = NULL,
    max_mt = 10
  )
  s_obj_norm <- normalize_counts(x = s_obj_filt, method = "sct", assay = "RNA")
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

  metadata <- as_data_frame_seurat(
    x = s_obj_dr,
    metadata = FALSE,
    reduction = "pcatest"
  )

  expect_equal(ncol(metadata), 21)
})

test_that("as_data_frame_seurat converts assay to dataframe", {
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
    min_genes = 1,
    min_counts = 1,
    max_genes = NULL,
    max_mt = 10
  )
  s_obj_norm <- normalize_counts(x = s_obj_filt, method = "sct", assay = "RNA")
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

  metadata <- as_data_frame_seurat(
    x = s_obj_dr,
    metadata = FALSE,
    assay = "RNA",
    slot = "counts",
    features = c("MAP4", "FTL")
  )

  expect_equal(ncol(metadata), 3)
})

test_that("as_data_frame_seurat converts metadata, reduction, and assay to dataframe", {
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
    min_genes = 1,
    min_counts = 1,
    max_genes = NULL,
    max_mt = 10
  )
  s_obj_norm <- normalize_counts(x = s_obj_filt, method = "sct", assay = "RNA")
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

  metadata <- as_data_frame_seurat(
    x = s_obj_dr,
    metadata = TRUE,
    assay = "RNA",
    slot = "counts",
    features = c("MAP4", "FTL"),
    reduction = "pcatest"
  )

  expect_equal(ncol(metadata), 36)
})

test_that("as_data_frame_seurat converts metadata, reduction, and assay to dataframe", {
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
    min_genes = 1,
    min_counts = 1,
    max_genes = NULL,
    max_mt = 10
  )

  # normalize_counts() already selects variable features and scales for "log", matching "sct",
  # so nothing further is needed before run_dr() — same as the sct tests above
  s_obj_norm <- normalize_counts(x = s_obj_filt, method = "log", assay = "RNA")

  s_obj_dr <- run_dr(
    s_obj_norm,
    dr_method = "pca",
    suffix = "test",
    var_features = TRUE,
    num_pcs = 20,
    assay = "RNA"
  )

  s_obj_dr <- run_dr(
    s_obj_dr,
    dr_method = "umap",
    suffix = "test",
    reduction = "pcatest",
    num_dim = 20,
    assay = "RNA",
    num_neighbors = 6
  )

  metadata <- as_data_frame_seurat(
    x = s_obj_dr,
    metadata = TRUE,
    assay = "RNA",
    slot = "scale.data",
    features = c("MAP4", "FTL"),
    reduction = c("pcatest", "umaptest")
  )

  expect_equal(ncol(metadata), 36)
})

test_that("metadata and reduction embeddings are collected into one table", {
  test_data <- SeuratObject::pbmc_small

  meta_tbl <- save_metadata(
    test_data,
    reduction = c("tsne", "pca"),
    file = NULL
  )

  # one row per cell, joined on a `cell` column
  expect_s3_class(meta_tbl, "tbl_df")
  expect_equal(nrow(meta_tbl), ncol(test_data))
  expect_true("cell" %in% colnames(meta_tbl))
  expect_setequal(meta_tbl$cell, colnames(test_data))

  # metadata columns match the Seurat object as-is, with embeddings joined in
  expect_setequal(
    setdiff(colnames(meta_tbl), colnames(test_data@meta.data)),
    c("cell", grep("^tSNE_|^PC_", colnames(meta_tbl), value = TRUE))
  )
  expect_true(any(grepl("^tSNE_", colnames(meta_tbl))))
  expect_true(any(grepl("^PC_", colnames(meta_tbl))))

  # embeddings are rounded, metadata counts are not
  expect_equal(meta_tbl$tSNE_1, round(meta_tbl$tSNE_1, 3))
  expect_equal(sort(meta_tbl$nCount_RNA), sort(unname(test_data$nCount_RNA)))
})

test_that("save_metadata skips reductions that are absent", {
  test_data <- SeuratObject::pbmc_small
  # umap is not in pbmc_small, so it must be skipped rather than error
  meta_tbl <- save_metadata(
    test_data,
    reduction = c("tsne", "umap"),
    file = NULL
  )

  expect_equal(nrow(meta_tbl), ncol(test_data))
  expect_true(any(grepl("^tSNE_", colnames(meta_tbl))))
  # the umap key is lower case ("umap_1"), so match case insensitively
  expect_false(any(grepl("umap", colnames(meta_tbl), ignore.case = TRUE)))
})

test_that("save_metadata writes a csv when file is not NULL", {
  test_data <- SeuratObject::pbmc_small
  csv_file <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(csv_file))

  meta_tbl <- save_metadata(test_data, reduction = "tsne", file = csv_file)

  expect_true(file.exists(csv_file))
  csv_tbl <- readr::read_csv(csv_file, show_col_types = FALSE)
  expect_equal(nrow(csv_tbl), nrow(meta_tbl))
})

test_that("save_counts collects the counts matrix into one table", {
  test_data <- SeuratObject::pbmc_small
  csv_file <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(csv_file))

  counts_tbl <- save_counts(test_data, file = csv_file)

  # one row per gene, one column per cell plus the gene names
  expect_s3_class(counts_tbl, "tbl_df")
  expect_equal(nrow(counts_tbl), nrow(test_data))
  expect_equal(ncol(counts_tbl), ncol(test_data) + 1)
  expect_setequal(counts_tbl$gene, rownames(test_data))

  expect_true(file.exists(csv_file))
  csv_tbl <- readr::read_csv(csv_file, show_col_types = FALSE)
  expect_equal(nrow(csv_tbl), nrow(counts_tbl))
})

test_that("save_counts requires a file to be named", {
  # reported before any of the work, rather than after densifying the matrix
  expect_error(save_counts(SeuratObject::pbmc_small), "\"file\" is missing")
})

test_that("save_counts leaves whole numbers unrounded", {
  test_data <- SeuratObject::pbmc_small
  raw_counts <- SeuratObject::GetAssayData(test_data, layer = "counts")

  # digits = -2 would round the counts to the nearest hundred if rounding were applied
  counts_tbl <- save_counts(
    test_data,
    layer = "counts",
    digits = -2,
    file = NULL
  )

  expect_equal(sum(counts_tbl[, -1]), sum(raw_counts))
})

test_that("save_counts rounds values that are not whole numbers", {
  test_data <- SeuratObject::pbmc_small

  data_tbl <- save_counts(test_data, layer = "data", digits = 1, file = NULL)
  data_values <- unlist(data_tbl[, -1], use.names = FALSE)

  # the normalized layer really is fractional, and comes back rounded
  expect_true(any(data_values != trunc(data_values)))
  expect_equal(data_values, round(data_values, 1))
})

test_that("save_counts refuses a matrix over the 2^31 element limit", {
  # 2.5e9 elements but only 50000 nonzeros, so the fixture is cheap to build
  n <- 50000
  big_obj <- SeuratObject::CreateSeuratObject(Matrix::sparseMatrix(
    i = seq_len(n),
    j = seq_len(n),
    x = 1,
    dims = c(n, n),
    dimnames = list(paste0("g", seq_len(n)), paste0("c", seq_len(n)))
  ))

  # the element limit is reported ahead of the memory one, since no machine gets around it
  expect_error(
    save_counts(big_obj, layer = "counts", file = NULL),
    "2\\^31 element limit"
  )
})

test_that("save_counts rejoins the layers left split by a merge", {
  pbmc_mat <- get_test_counts_matrix()
  half <- ncol(pbmc_mat) %/% 2
  merged_obj <- merge(
    initialize_seurat_object(counts_matrix = pbmc_mat[, 1:half]),
    initialize_seurat_object(
      counts_matrix = pbmc_mat[, (half + 1):ncol(pbmc_mat)]
    )
  )

  # merge() leaves the layers split, which GetAssayData() refuses to read
  expect_gt(
    length(SeuratObject::Layers(merged_obj[["RNA"]], search = "counts")),
    1
  )

  counts_tbl <- save_counts(merged_obj, layer = "counts", file = NULL)

  expect_equal(ncol(counts_tbl), ncol(merged_obj) + 1)
  expect_setequal(counts_tbl$gene, rownames(merged_obj))
})

test_that("resolve_seurat_object passes an in-memory Seurat object through unchanged", {
  pbmc_obj <- initialize_seurat_object(counts_matrix = get_test_counts_matrix())
  expect_identical(resolve_seurat_object(pbmc_obj), pbmc_obj)
})

test_that("resolve_seurat_object reads a .qs2 file", {
  pbmc_obj <- initialize_seurat_object(counts_matrix = get_test_counts_matrix())
  qs2_file <- tempfile(fileext = ".qs2")
  on.exit(unlink(qs2_file))
  qs2::qs_save(pbmc_obj, file = qs2_file)

  expect_equal(ncol(resolve_seurat_object(qs2_file)), ncol(pbmc_obj))
})

test_that("resolve_seurat_object reads an .rds file", {
  pbmc_obj <- initialize_seurat_object(counts_matrix = get_test_counts_matrix())
  rds_file <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_file))
  saveRDS(pbmc_obj, file = rds_file)

  expect_equal(ncol(resolve_seurat_object(rds_file)), ncol(pbmc_obj))
})

test_that("resolve_seurat_object errors on a missing file", {
  expect_error(resolve_seurat_object("does-not-exist.qs2"), "path not found")
})

test_that("resolve_seurat_object errors on an unsupported extension", {
  txt_file <- tempfile(fileext = ".txt")
  on.exit(unlink(txt_file))
  writeLines("not a seurat object", txt_file)

  expect_error(resolve_seurat_object(txt_file), "unsupported file extension")
})

test_that("resolve_seurat_object errors when the file does not contain a Seurat object", {
  qs2_file <- tempfile(fileext = ".qs2")
  on.exit(unlink(qs2_file))
  qs2::qs_save(list(not = "a seurat object"), file = qs2_file)

  expect_error(
    resolve_seurat_object(qs2_file),
    "did not contain a Seurat object"
  )
})

test_that("resolve_seurat_object reads the single qs2/rds file in a directory", {
  pbmc_obj <- initialize_seurat_object(counts_matrix = get_test_counts_matrix())
  analysis_dir <- tempfile()
  dir.create(analysis_dir)
  on.exit(unlink(analysis_dir, recursive = TRUE))
  qs2::qs_save(pbmc_obj, file.path(analysis_dir, "seurat_obj.qs2"))

  expect_equal(ncol(resolve_seurat_object(analysis_dir)), ncol(pbmc_obj))
})

test_that("resolve_seurat_object errors on a directory with no qs2/rds file", {
  empty_dir <- tempfile()
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE))

  expect_error(resolve_seurat_object(empty_dir), "no .rds or .qs2 file")
})

test_that("resolve_seurat_object errors on a directory with multiple qs2/rds files", {
  pbmc_obj <- initialize_seurat_object(counts_matrix = get_test_counts_matrix())
  analysis_dir <- tempfile()
  dir.create(analysis_dir)
  on.exit(unlink(analysis_dir, recursive = TRUE))
  qs2::qs_save(pbmc_obj, file.path(analysis_dir, "seurat_obj.qs2"))
  saveRDS(pbmc_obj, file.path(analysis_dir, "backup.rds"))

  expect_error(resolve_seurat_object(analysis_dir), "multiple .rds/.qs2 files")
})
