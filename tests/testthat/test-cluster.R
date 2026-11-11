test_that("cluster stats and expression functions write their tables", {
  withr::local_dir(withr::local_tempdir())
  test_data <- SeuratObject::pbmc_small
  test_data <- RunUMAP(test_data, dims = 1:5, verbose = FALSE)
  test_data$total_counts <- test_data$nCount_RNA
  test_data$detected_genes <- test_data$nFeature_RNA
  test_data$pct_mito <- 0
  test_data$orig.ident <- test_data$groups

  calculate_cluster_stats(test_data, label = "test")
  expect_true(file.exists("metadata-test.csv.gz"))
  expect_true(file.exists("summary-test.csv"))
  expect_true(file.exists("summary-test-per-sample.csv"))

  # tsne and umap are each optional
  umap_only <- test_data
  umap_only[["tsne"]] <- NULL
  calculate_cluster_stats(umap_only, label = "umap")
  metadata_umap <- readr::read_csv(
    "metadata-umap.csv.gz",
    show_col_types = FALSE
  )
  expect_false(any(grepl("tsne", colnames(metadata_umap), ignore.case = TRUE)))
  expect_true(any(grepl("umap", colnames(metadata_umap), ignore.case = TRUE)))

  umap_only[["umap"]] <- NULL
  calculate_cluster_stats(umap_only, label = "none")
  expect_false(file.exists("metadata-none.csv.gz"))
  expect_true(file.exists("summary-none.csv"))

  calculate_cluster_expression(test_data, label = "test")
  expect_true(file.exists("expression-mean-test.csv"))
  expect_true(file.exists("expression-mean-test-per-sample.csv"))
  expect_true(file.exists("expression-raw-test-per-sample.csv"))
  expression_raw <- readr::read_csv(
    "expression-raw-test.csv",
    show_col_types = FALSE
  )
  expect_equal(
    sum(expression_raw[, -1]),
    sum(SeuratObject::LayerData(test_data, layer = "counts"))
  )
})

test_that("calculate_cluster_markers writes markers and heatmaps", {
  withr::local_dir(withr::local_tempdir())
  test_data <- SeuratObject::pbmc_small
  # plot_cluster_markers_top() plots on the umap reduction
  test_data <- RunUMAP(test_data, dims = 1:5, verbose = FALSE)

  suppressWarnings(
    calculate_cluster_markers(test_data, label = "test", test = "wilcox")
  )
  all_markers <- readr::read_csv(
    "markers-global/markers-test-wilcox-all.csv",
    show_col_types = FALSE
  )
  expect_named(
    all_markers,
    c("cluster", "gene", "log2FC", "p_val", "p_val_adj")
  )
  expect_true(file.exists(
    "markers-global/markers-test-wilcox-heatmap-top5.csv"
  ))

  suppressWarnings(
    calculate_cluster_markers(
      test_data,
      label = "test",
      test = "wilcox",
      pairwise = TRUE
    )
  )
  expect_true(
    file.exists("markers-pairwise/markers-test-pairwise-wilcox-all.csv")
  )
})

test_that("differential_expression_per_cluster writes DE tables", {
  withr::local_dir(withr::local_tempdir())
  test_data <- SeuratObject::pbmc_small

  suppressWarnings(
    differential_expression_per_cluster(
      test_data,
      cluster_column = "RNA_snn_res.1",
      group_column = "groups"
    )
  )
  de_files <- list.files("differential-expression", pattern = "-all\\.csv$")
  expect_gt(length(de_files), 0)
})
