test_that("Seurat object can be split by batch", {
  # same fixture strategy Seurat itself uses for integration: pbmc_small converted to a v5 assay
  # and split by the `groups` column (g1: 44 cells, g2: 36 cells)
  test_data <- SeuratObject::pbmc_small
  # real objects always carry pct_mito (initialize_seurat_object() adds it); pbmc_small has none,
  # so give it one here for realism, even though normalize_counts() does not regress it out
  set.seed(1)
  test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
  suppressWarnings(
    test_data[["RNA"]] <- SeuratObject::CreateAssay5Object(
      counts = SeuratObject::LayerData(
        test_data,
        assay = "RNA",
        layer = "counts"
      )
    )
  )

  split_obj <- suppressMessages(
    split_layers_by_batch(
      test_data,
      batch_var = "groups",
      num_variable_genes = 100,
      num_dim = 10
    )
  )

  # the layers are split, not the object
  expect_s4_class(split_obj, "Seurat")
  expect_setequal(
    SeuratObject::Layers(split_obj[["RNA"]], search = "counts"),
    c("counts.g1", "counts.g2")
  )
  # splitting must not lose or duplicate cells
  expect_equal(ncol(split_obj), ncol(test_data))
  # the reduction that IntegrateLayers() corrects
  expect_true("pca" %in% SeuratObject::Reductions(split_obj))
})

test_that("variable features are reported per batch", {
  test_data <- SeuratObject::pbmc_small
  # real objects always carry pct_mito (initialize_seurat_object() adds it); pbmc_small has none,
  # so give it one here for realism, even though normalize_counts() does not regress it out
  set.seed(1)
  test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
  suppressWarnings(
    test_data[["RNA"]] <- SeuratObject::CreateAssay5Object(
      counts = SeuratObject::LayerData(
        test_data,
        assay = "RNA",
        layer = "counts"
      )
    )
  )
  split_obj <- suppressMessages(
    split_layers_by_batch(
      test_data,
      batch_var = "groups",
      num_variable_genes = 100,
      num_dim = 10
    )
  )

  var_genes <- variable_features_by_batch(split_obj)

  expect_length(var_genes, 2)
  expect_setequal(names(var_genes), c("g1", "g2"))
  # every batch contributes features, and they are real genes
  expect_true(all(lengths(var_genes) > 0))
  expect_true(all(unlist(var_genes) %in% rownames(split_obj)))
})

test_that("variable features by batch requires a split object", {
  expect_error(
    variable_features_by_batch(SeuratObject::pbmc_small),
    "split_layers_by_batch"
  )
})

test_that("Seurat object layers can be integrated", {
  test_data <- SeuratObject::pbmc_small
  # real objects always carry pct_mito (initialize_seurat_object() adds it); pbmc_small has none,
  # so give it one here for realism, even though normalize_counts() does not regress it out
  set.seed(1)
  test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
  suppressWarnings(
    test_data[["RNA"]] <- SeuratObject::CreateAssay5Object(
      counts = SeuratObject::LayerData(
        test_data,
        assay = "RNA",
        layer = "counts"
      )
    )
  )
  split_obj <- suppressMessages(
    split_layers_by_batch(
      test_data,
      batch_var = "groups",
      num_variable_genes = 100,
      num_dim = 10
    )
  )

  int_obj <- suppressWarnings(suppressMessages(
    integrate_layers(
      split_obj,
      num_dim = 10,
      int_reduction = "rpca",
      k_weight = 20
    )
  ))

  # v5 integration produces a reduction, not an "integrated" assay
  expect_true("rpca" %in% SeuratObject::Reductions(int_obj))
  expect_false("integrated" %in% SeuratObject::Assays(int_obj))
  expect_equal(
    nrow(SeuratObject::Embeddings(int_obj, "rpca")),
    ncol(test_data)
  )

  # layers are rejoined, so GetAssayData() works again
  expect_setequal(
    SeuratObject::Layers(int_obj[["RNA"]], search = "counts"),
    "counts"
  )
  expect_equal(
    ncol(Seurat::GetAssayData(int_obj, assay = "RNA", layer = "counts")),
    ncol(test_data)
  )
})

test_that("k_weight is reduced to fit the smallest batch", {
  # pbmc_small$orig.ident is a single value, so the batch sizes have to come from the layers:
  # the smallest `groups` batch is 36, well under the default k_weight of 100
  test_data <- SeuratObject::pbmc_small
  # real objects always carry pct_mito (initialize_seurat_object() adds it); pbmc_small has none,
  # so give it one here for realism, even though normalize_counts() does not regress it out
  set.seed(1)
  test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
  suppressWarnings(
    test_data[["RNA"]] <- SeuratObject::CreateAssay5Object(
      counts = SeuratObject::LayerData(
        test_data,
        assay = "RNA",
        layer = "counts"
      )
    )
  )
  split_obj <- suppressMessages(
    split_layers_by_batch(
      test_data,
      batch_var = "groups",
      num_variable_genes = 100,
      num_dim = 10
    )
  )

  # would error with "k.weight is set larger than the number of cells" if not clamped
  int_obj <- suppressWarnings(suppressMessages(
    integrate_layers(split_obj, num_dim = 10, int_reduction = "rpca")
  ))
  expect_true("rpca" %in% SeuratObject::Reductions(int_obj))
})

test_that("integration rejects invalid arguments", {
  expect_error(integrate_layers(list(), num_dim = 2), "too few dims")
  expect_error(integrate_layers(list(), num_dim = 99), "too many dims")
  expect_error(
    integrate_layers(list(), num_dim = 10, int_reduction = "nonsense"),
    "is not valid"
  )
})

test_that("merging requires at least two objects", {
  expect_error(
    merge_seurat_objects(list(SeuratObject::pbmc_small)),
    "at least 2 samples"
  )
})

test_that("merging accepts a mix of objects and analysis directory paths", {
  make_sample <- function(sample_name) {
    test_data <- SeuratObject::pbmc_small
    # real objects always carry these (initialize_seurat_object() adds them); pbmc_small has none
    set.seed(1)
    test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
    test_data$total_counts <- test_data$nCount_RNA
    test_data$detected_genes <- test_data$nFeature_RNA
    suppressWarnings(
      test_data[["RNA"]] <- SeuratObject::CreateAssay5Object(
        counts = SeuratObject::LayerData(
          test_data,
          assay = "RNA",
          layer = "counts"
        )
      )
    )
    test_data$orig.ident <- factor(sample_name)
    SeuratObject::RenameCells(test_data, add.cell.id = sample_name)
  }

  obj1 <- make_sample("s1")
  obj2 <- make_sample("s2")

  # obj2 only exists on disk, as an analysis dir would - obj1 stays in memory
  analysis_dir <- file.path(tempdir(), "s2-analysis-dir")
  dir.create(analysis_dir)
  qs2::qs_save(obj2, file.path(analysis_dir, "seurat_obj.qs2"))

  merged <- merge_seurat_objects(list(obj1, analysis_dir), min_cells = 1)

  expect_s4_class(merged, "Seurat")
  expect_equal(ncol(merged), ncol(obj1) + ncol(obj2))
  expect_setequal(as.character(unique(merged$orig.ident)), c("s1", "s2"))

  expect_error(
    merge_seurat_objects(list(obj1, file.path(tempdir(), "no-such-dir"))),
    "path not found"
  )
})

test_that("integration stops when a batch is too small to weight", {
  test_data <- SeuratObject::pbmc_small
  # real objects always carry pct_mito (initialize_seurat_object() adds it); pbmc_small has none,
  # so give it one here for realism, even though normalize_counts() does not regress it out
  set.seed(1)
  test_data$pct_mito <- stats::runif(ncol(test_data), 0, 10)
  suppressWarnings(
    test_data[["RNA"]] <- SeuratObject::CreateAssay5Object(
      counts = SeuratObject::LayerData(
        test_data,
        assay = "RNA",
        layer = "counts"
      )
    )
  )
  # 80 cells split 70/10, so the smaller batch falls under the 25 cell floor
  test_data$tiny <- rep(c("big", "small"), times = c(70, 10))
  split_obj <- suppressMessages(
    split_layers_by_batch(
      test_data,
      batch_var = "tiny",
      num_variable_genes = 100,
      num_dim = 10
    )
  )

  expect_error(
    integrate_layers(split_obj, num_dim = 10, int_reduction = "rpca"),
    "too few to integrate"
  )
})
