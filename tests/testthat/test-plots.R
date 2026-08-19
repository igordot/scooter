context("test plotting-related functions")

test_that("color scheme", {
  colors <- get_color_scheme()
  # colors <- getOption("scooter.colors_samples")
  expect_type(colors, "character")
  expect_gt(length(colors), 50)
  colors <- get_color_scheme(type = "samples")
  # colors <- getOption("scooter.colors_clusters")
  expect_type(colors, "character")
  expect_gt(length(colors), 10)
})

test_that("QC plots cover every metric and share their column check", {
  test_data <- SeuratObject::pbmc_small
  test_data$total_counts <- test_data$nCount_RNA
  test_data$detected_genes <- test_data$nFeature_RNA
  # a real spread, not a constant: FeatureScatter titles each panel with a cor(), which warns on
  # a zero-variance input
  test_data$pct_mito <- 100 *
    test_data$total_counts /
    max(test_data$total_counts)
  test_data$pct_ribo <- 100 *
    test_data$detected_genes /
    max(test_data$detected_genes)
  test_data$pct_hb <- 100 *
    test_data$total_counts /
    max(test_data$total_counts)

  # one scatter panel per pair of metrics, one violin panel per metric; plot_metrics_distribution()
  # has no default metric set - the caller always states it explicitly
  five_metrics <- c(
    "detected_genes",
    "total_counts",
    "pct_mito",
    "pct_ribo",
    "pct_hb"
  )
  expect_length(plot_metrics_correlations(test_data)$layers, 3)
  expect_length(
    plot_metrics_distribution(test_data, metrics = five_metrics)$layers,
    5
  )

  # both report a missing metadata column rather than failing inside Seurat
  expect_error(
    plot_metrics_correlations(
      test_data,
      pairs = list(c("total_counts", "absent"))
    ),
    "metadata columns not found: absent"
  )
  expect_error(
    plot_metrics_distribution(
      test_data,
      metrics = c("total_counts", "absent")
    ),
    "metadata columns not found: absent"
  )
})

test_that("dimensionality reduction plot point size", {
  point_size <- get_dr_point_size(10)
  expect_type(point_size, "double")
  expect_gt(point_size, 0)
  expect_lt(point_size, 10)
  point_size <- get_dr_point_size(1000)
  expect_type(point_size, "double")
  expect_gt(point_size, 0)
  expect_lt(point_size, 10)
  point_size <- get_dr_point_size(100000)
  expect_type(point_size, "double")
  expect_gt(point_size, 0)
  expect_lt(point_size, 10)

  # larger datasets get a smaller point size
  expect_gt(get_dr_point_size(10), get_dr_point_size(100000))

  # accepts a Seurat object (uses ncol) or a data frame (uses nrow), matching the count form
  seurat_obj <- SeuratObject::pbmc_small
  expect_equal(
    get_dr_point_size(seurat_obj),
    get_dr_point_size(ncol(seurat_obj))
  )
  df <- data.frame(cell = seq_len(2000))
  expect_equal(get_dr_point_size(df), get_dr_point_size(nrow(df)))
  # a tibble (class length > 1) must not error
  tbl <- tibble::tibble(cell = seq_len(20000))
  expect_equal(get_dr_point_size(tbl), get_dr_point_size(20000))
})

test_that("filter_cells on a Seurat object always writes its QC plots", {
  pbmc_mat <- get_test_counts_matrix()
  s_obj <- initialize_seurat_object(counts_matrix = pbmc_mat)

  # automatic output is the point (see setup.R for the working-directory isolation this relies on)
  filter_cells(
    s_obj,
    min_genes = 1,
    min_counts = 1,
    max_mt = 10,
    log_file = NULL
  )
  # each plot type writes its own file, flat into the working directory
  expect_true(all(
    c(
      "metrics-distribution-unfiltered.png",
      "metrics-distribution-filtered.png",
      "metrics-correlations-filtered.png"
    ) %in%
      list.files()
  ))
})

test_that("the QC plot functions save themselves when given a file", {
  test_data <- SeuratObject::pbmc_small
  test_data$total_counts <- test_data$nCount_RNA
  test_data$detected_genes <- test_data$nFeature_RNA
  set.seed(1)
  test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
  test_data$pct_ribo <- stats::runif(ncol(test_data), 0, 10)
  three_metrics <- c("detected_genes", "total_counts", "pct_mito")
  out <- file.path(tempdir(), "qc-selfsave")
  on.exit(unlink(out, recursive = TRUE))

  # the plot is returned either way; a path also writes it, creating the directory
  expect_s3_class(
    plot_metrics_distribution(test_data, metrics = three_metrics),
    "ggplot"
  )
  expect_false(dir.exists(out))

  p <- plot_metrics_distribution(
    test_data,
    metrics = three_metrics,
    file = file.path(out, "dist.png")
  )
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(file.path(out, "dist.png")))

  plot_metrics_correlations(test_data, file = file.path(out, "corr.png"))
  expect_true(file.exists(file.path(out, "corr.png")))
})
