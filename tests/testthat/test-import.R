test_that("import_mtx can read in 10x data", {
  counts <- import_mtx(
    data_path = system.file(
      "extdata",
      "outs/filtered_feature_bc_matrix",
      package = "scooter"
    ),
    gene_column = 2
  )
  truth <- readRDS(system.file(
    "extdata",
    "import_matrix_counts.rds",
    package = "scooter"
  ))
  expect_identical(counts, truth)
})

test_that("read_counts_file can read in 10x data", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  expect_identical(names(counts), c("Antibody Capture", "Gene Expression"))
})

test_that("read_counts_file can read in tsv Antibody Capture file", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "HTO.tsv", package = "scooter")
  )
  expect_identical(
    rownames(counts$`Antibody Capture`),
    c("Hashtag1", "Hashtag2")
  )
})

test_that("scooter will error is there is no valid path", {
  expect_error(read_counts_file(
    sample_name = "test",
    path = system.file("extdata-test", "", package = "scooter")
  ))
})

test_that("scooter will error is there is no mtx in directory", {
  expect_error(read_counts_file(
    sample_name = "test",
    path = system.file("R", "", package = "scooter")
  ))
})

test_that("read_counts_file errors and lists files when multiple h5 files are found", {
  h5_dir <- file.path(tempdir(), "multi-h5")
  dir.create(h5_dir, showWarnings = FALSE)
  file.create(file.path(h5_dir, "sample1.h5"))
  file.create(file.path(h5_dir, "sample2.h5"))
  expect_error(
    read_counts_file(sample_name = "test", path = h5_dir),
    "multiple h5 files found"
  )
})

test_that("read_counts_file errors when no mtx or h5 file is found", {
  empty_dir <- file.path(tempdir(), "no-counts")
  dir.create(empty_dir, showWarnings = FALSE)
  expect_error(
    read_counts_file(sample_name = "test", path = empty_dir),
    "does not contain matrix.mtx or an h5 file"
  )
})

test_that("read_counts_file can read in csv Antibody Capture file", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "HTO.csv", package = "scooter")
  )
  expect_identical(names(counts), c("Antibody Capture"))
})

test_that("scooter will remove -1 from barcodes from txt file", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "pbmc.txt", package = "scooter")
  )
  expect_that(
    colnames(counts$`Gene Expression`)[1],
    equals("test:GGAATCTGCTTAGG")
  )
})

test_that("Seurat object can be created from RNA data", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  # only 462 genes are above zero in the test data
  expect_s4_class(s_obj, "Seurat")
})

test_that("Seurat object can be created from ADT data", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Antibody Capture`,
    assay = "ADT",
    log_file = NULL
  )
  expect_s4_class(s_obj, "Seurat")
})

test_that("Seurat obj can be created", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)
  expect_s4_class(pbmc_obj, "Seurat")
})

test_that("Assay can be added to Seurat object", {
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
  expect_s4_class(s_obj, "Seurat")
})

test_that("Assay cannot be added to a random object", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  expect_error(add_seurat_assay(
    x = counts,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  ))
})

test_that("Seurat object can be filtered", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)
  pbmc_obj <- filter_cells(
    pbmc_obj,
    log_file = NULL,
    min_genes = NULL,
    max_genes = NULL,
    max_mt = 10
  )
  expect_lt(ncol(pbmc_obj), 60)
})

test_that("Seurat object can be filtered", {
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
  expect_lt(ncol(s_obj_filt), ncol(s_obj))
})

test_that("Seurat object can be log normalized", {
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
  s_obj_norm <- normalize_counts(
    x = s_obj_filt,
    method = "log",
    assay = "RNA"
  )

  # normalize_counts() now matches SCTransform()'s one-call pattern: normalize, select variable
  # features, and scale, all in one call
  expect_true("data" %in% SeuratObject::Layers(s_obj_norm[["RNA"]]))
  expect_true("scale.data" %in% SeuratObject::Layers(s_obj_norm[["RNA"]]))
  expect_gt(length(SeuratObject::VariableFeatures(s_obj_norm)), 0)
  expect_gt(
    ncol(GetAssayData(s_obj_norm, assay = "RNA", layer = "scale.data")),
    150
  )
})

test_that("Seurat object can be SC transform", {
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
  s_obj_norm <- normalize_counts(
    x = s_obj_filt,
    method = "sct",
    assay = "RNA"
  )
  expect_gt(
    ncol(GetAssayData(s_obj_norm, assay = "SCT", layer = "scale.data")),
    150
  )
})

test_that("mad_outlier_bounds returns median +/- nmads * MAD", {
  v <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
  b <- mad_outlier_bounds(v, nmads = 3, log = FALSE)
  expect_named(b, c("min", "max"))
  expect_equal(b[["min"]], median(v) - 3 * mad(v))
  expect_equal(b[["max"]], median(v) + 3 * mad(v))
})

test_that("mad_outlier_bounds on the log scale is positive and back-transformed", {
  set.seed(1)
  v <- round(exp(rnorm(1000, mean = log(500), sd = 0.4)))
  b <- mad_outlier_bounds(v, nmads = 3, log = TRUE)
  lx <- log(v)
  expect_equal(b[["min"]], exp(median(lx) - 3 * mad(lx)))
  expect_equal(b[["max"]], exp(median(lx) + 3 * mad(lx)))
  # lower bound stays positive and below the median, upper stays above it
  expect_gt(b[["min"]], 0)
  expect_lt(b[["min"]], median(v))
  expect_gt(b[["max"]], median(v))
})

test_that("filter_cells floors the data-driven lower cutoffs at 300 genes / 500 counts", {
  # a wide gene-count distribution whose MAD lower bounds fall below the floors
  set.seed(1)
  ngenes <- 1000
  ncells <- 400
  mat <- matrix(
    0L,
    ngenes,
    ncells,
    dimnames = list(paste0("g", seq_len(ngenes)), paste0("c", seq_len(ncells)))
  )
  n_expr <- pmin(
    pmax(round(exp(rnorm(ncells, mean = log(400), sd = 0.6))), 20),
    ngenes
  )
  for (j in seq_len(ncells)) {
    mat[sample(ngenes, n_expr[j]), j] <- rpois(n_expr[j], 3) + 1L
  }
  obj <- suppressWarnings(initialize_seurat_object(counts_matrix = mat))
  expect_gte(nrow(obj), 1000) # whole-transcriptome, so the 300-gene floor applies

  # confirm the MAD lower bounds are actually below the floors on this distribution
  expect_lt(
    mad_outlier_bounds(obj$detected_genes, nmads = 3, log = TRUE)[["min"]],
    300
  )
  expect_lt(
    mad_outlier_bounds(obj$total_counts, nmads = 3, log = TRUE)[["min"]],
    500
  )

  # the fixed cutoff is applied after the MAD stage, so it binds where the MAD bound was looser
  kept <- suppressWarnings(filter_cells(obj, max_mt = 100, log_file = NULL))
  expect_gte(min(kept$detected_genes), 500)

  # min_genes = NULL leaves only the MAD stage, which lets those same cells through
  kept_mad_only <- suppressWarnings(
    filter_cells(obj, min_genes = NULL, max_mt = 100, log_file = NULL)
  )
  expect_lt(min(kept_mad_only$detected_genes), 500)
  expect_gt(ncol(kept_mad_only), ncol(kept))
})

test_that("filter_cells num_mads controls how aggressive the outlier stage is", {
  set.seed(3)
  n <- 600
  metadata <- data.frame(
    total_counts = round(exp(rnorm(n, log(3000), 0.5))),
    detected_genes = round(exp(rnorm(n, log(1200), 0.5))),
    pct_mito = 0,
    row.names = paste0("c", seq_len(n))
  )

  # only the outlier stage, so the cell counts isolate its effect
  strict <- filter_cells(
    metadata,
    num_mads = 1,
    min_genes = NULL,
    log_file = NULL
  )
  loose <- filter_cells(
    metadata,
    num_mads = 5,
    min_genes = NULL,
    log_file = NULL
  )
  none <- filter_cells(
    metadata,
    num_mads = NULL,
    min_genes = NULL,
    log_file = NULL
  )

  expect_lt(length(strict), length(loose))
  expect_lte(length(loose), length(none))
  expect_length(none, n)
})

test_that("filter_cells stops when too few cells survive filtering", {
  pbmc_obj <- suppressWarnings(initialize_seurat_object(
    counts_matrix = get_test_counts_matrix()
  ))

  # the default min_cells (50) is not cleared by this fixture at the default min_genes (500)
  expect_error(
    filter_cells(pbmc_obj, log_file = NULL, max_mt = 100),
    "below the minimum of 50 required to proceed"
  )

  # raising min_genes only shrinks the surviving set further, so the same cutoff still trips
  expect_error(
    filter_cells(pbmc_obj, log_file = NULL, max_mt = 100, min_cells = 1000),
    "below the minimum of 1000 required to proceed"
  )

  # min_cells = NULL is an explicit opt-out
  kept <- filter_cells(
    pbmc_obj,
    log_file = NULL,
    max_mt = 100,
    min_cells = NULL
  )
  expect_gt(length(kept), 0)
})

test_that("filter_cells applies min_genes regardless of the panel size", {
  # a targeted panel: few genes total, with a wide per-cell distribution
  set.seed(2)
  ngenes <- 500
  ncells <- 400
  mat <- matrix(
    0L,
    ngenes,
    ncells,
    dimnames = list(paste0("g", seq_len(ngenes)), paste0("c", seq_len(ncells)))
  )
  n_expr <- pmin(
    pmax(round(exp(rnorm(ncells, mean = log(150), sd = 0.6))), 10),
    ngenes
  )
  for (j in seq_len(ncells)) {
    mat[sample(ngenes, n_expr[j]), j] <- rpois(n_expr[j], 3) + 1L
  }
  obj <- suppressWarnings(initialize_seurat_object(counts_matrix = mat))
  expect_lt(nrow(obj), 1000) # a targeted panel: no cell can reach 500 genes

  gene_lower <- round(mad_outlier_bounds(
    obj$detected_genes,
    nmads = 3,
    log = TRUE
  )[["min"]])
  expect_lt(gene_lower, 300) # the MAD lower bound is well below the default cutoff

  # there is no longer a panel exemption: the default min_genes still applies, and on data at this
  # scale it removes almost everything (min_cells = NULL: that near-total removal is the point of
  # this assertion, not a production run the 50-cell floor should guard)
  kept_default <- filter_cells(
    obj@meta.data,
    max_mt = 100,
    min_cells = NULL,
    log_file = NULL
  )
  kept_no_min <- filter_cells(
    obj@meta.data,
    min_genes = NULL,
    max_mt = 100,
    log_file = NULL
  )
  expect_lt(length(kept_default), 0.1 * length(kept_no_min))
  expect_gte(min(obj@meta.data[kept_default, "detected_genes"]), 500)

  # a panel needs min_genes set for its scale (or NULL), leaving the MAD stage in charge
  kept <- suppressWarnings(filter_cells(
    obj,
    min_genes = NULL,
    max_mt = 100,
    log_file = NULL
  ))
  expect_gte(min(kept$detected_genes), gene_lower)
  expect_lt(min(kept$detected_genes), 300)
})

test_that("filter_cells defaults to MAD gene and count cutoffs when none are given", {
  pbmc_mat <- get_test_counts_matrix()
  pbmc_obj <- suppressWarnings(initialize_seurat_object(
    counts_matrix = pbmc_mat
  ))

  # the resolved default cutoffs should match a 3-MAD log-scale rule on genes and counts
  expected_genes <- mad_outlier_bounds(
    pbmc_obj$detected_genes,
    nmads = 3,
    log = TRUE
  )
  expected_counts <- mad_outlier_bounds(
    pbmc_obj$total_counts,
    nmads = 3,
    log = TRUE
  )
  # min_cells = NULL: the fixture is too small to clear the 50-cell production floor, and this
  # test is only checking the MAD math, not simulating a viable run
  kept <- filter_cells(
    pbmc_obj,
    log_file = NULL,
    max_mt = 100,
    min_cells = NULL
  )

  expect_gte(min(kept$detected_genes), round(expected_genes[["min"]]))
  expect_lte(max(kept$detected_genes), round(expected_genes[["max"]]))
  expect_gte(min(kept$total_counts), round(expected_counts[["min"]]))
  expect_lte(max(kept$total_counts), round(expected_counts[["max"]]))
})

test_that("initialize_seurat_object adds the QC metrics the QC plots expect", {
  pbmc_mat <- get_test_counts_matrix()
  s_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)

  expect_true(all(
    c("detected_genes", "total_counts", "pct_mito") %in%
      colnames(s_obj@meta.data)
  ))

  # aliases of the Seurat originals, which filter_cells() still reads
  expect_equal(s_obj$total_counts, s_obj$total_counts)
  expect_equal(s_obj$detected_genes, s_obj$detected_genes)

  # so the object plots straight away, with no renaming at the call site
  expect_s3_class(
    plot_metrics_distribution(
      s_obj,
      metrics = c("detected_genes", "total_counts", "pct_mito")
    ),
    "ggplot"
  )
  expect_s3_class(plot_metrics_correlations(s_obj), "ggplot")
})

test_that("initialize_seurat_object adds sample_name as an orig.ident alias", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts[["Gene Expression"]],
    assay = "RNA"
  )

  expect_true("sample_name" %in% colnames(s_obj@meta.data))
  expect_equal(s_obj$sample_name, s_obj$orig.ident)
  expect_true(all(s_obj$sample_name == "test"))
})

test_that("filter_cells keeps the non-RNA assays", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts[["Gene Expression"]],
    assay = "RNA"
  )
  s_obj <- add_seurat_assay(
    s_obj,
    assay = "ADT",
    counts_matrix = counts[["Antibody Capture"]]
  )
  adt_features <- rownames(s_obj@assays$ADT)

  s_filt <- filter_cells(s_obj, min_genes = 1, min_counts = 1, max_mt = 10)

  # subset(features =) drops an assay entirely when none of its features are listed
  expect_true("ADT" %in% names(s_filt@assays))
  expect_setequal(rownames(s_filt@assays$ADT), adt_features)
  expect_equal(ncol(s_filt@assays$ADT), ncol(s_filt))
})

test_that("create_seurat_object goes from a sample path to a clusterable object", {
  s_obj <- suppressWarnings(create_seurat_object(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter"),
    min_genes = 1,
    min_counts = 1,
    max_mt = 10,
    num_pcs = 10
  ))

  # antibody capture is picked up automatically and survives the filtering
  expect_setequal(names(s_obj@assays), c("RNA", "ADT"))

  # every layer the downstream steps need, produced in one call
  expect_true(all(
    c("counts", "data", "scale.data") %in%
      SeuratObject::Layers(s_obj[["RNA"]])
  ))
  expect_gt(length(SeuratObject::VariableFeatures(s_obj)), 0)
  expect_setequal(names(s_obj@reductions), c("pca", "tsne", "umap"))
  expect_equal(ncol(SeuratObject::Embeddings(s_obj, "pca")), 10)

  # QC metrics are present, so the object is ready for the QC plots
  expect_true(all(
    c("detected_genes", "total_counts", "pct_mito") %in%
      colnames(s_obj@meta.data)
  ))

  # and it is genuinely ready for clustering
  expect_s4_class(
    calculate_clusters(
      s_obj,
      assay = "RNA",
      reduction = "pca",
      num_dim = 10,
      res = 0.5
    ),
    "Seurat"
  )
})

test_that("create_seurat_object runs PCA on the SCT assay for normalization_method = 'sct'", {
  s_obj <- suppressWarnings(create_seurat_object(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter"),
    min_genes = 1,
    min_counts = 1,
    max_mt = 10,
    normalization_method = "sct",
    num_pcs = 10
  ))

  # SCTransform() creates and activates "SCT"; normalize_counts() does not run its extra
  # variable-feature/scaling steps for method = "sct" since SCTransform() already did that
  expect_true("SCT" %in% names(s_obj@assays))
  expect_equal(SeuratObject::DefaultAssay(s_obj), "SCT")
  expect_setequal(names(s_obj@reductions), c("pca", "tsne", "umap"))
})

test_that("create_seurat_object errors when there is no gene expression data", {
  # a flat file with few features is inferred to be antibody capture, so this path produces a
  # counts list with no "Gene Expression" element at all
  expect_error(
    create_seurat_object(
      sample_name = "test",
      path = system.file("extdata", "HTO.tsv", package = "scooter")
    ),
    "no gene expression data"
  )
})
