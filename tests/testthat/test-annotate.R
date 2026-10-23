context("test annotation functions")

test_that("annotate_pan_human_azimuth() converts False to NA and writes outputs", {
  skip_if_not_installed("AzimuthAPI")

  query <- SeuratObject::pbmc_small
  query[["umap"]] <- CreateDimReducObject(
    embeddings = matrix(
      rnorm(2 * ncol(query)),
      ncol = 2,
      dimnames = list(colnames(query), c("UMAP_1", "UMAP_2"))
    ),
    key = "UMAP_",
    assay = "RNA"
  )

  # a stand-in for the server: 5 cells have an inconsistent hierarchy
  fake_cloud_azimuth <- function(object, assay) {
    cells <- colnames(object)
    inconsistent <- seq_along(cells) <= 5
    refined <- ifelse(inconsistent, "False", rep(c("A", "B"), length.out = 80))
    labels <- rep(c("A1", "B1"), length.out = 80)
    object$full_hierarchical_labels <- paste0("root|", labels)
    object$final_level_labels <- labels
    object$final_level_confidence <- runif(80)
    object$full_consistent_hierarchy <- !inconsistent
    object$azimuth_broad <- refined
    object$azimuth_medium <- refined
    object$azimuth_fine <- refined
    object$azimuth_label <- labels
    object[["azimuth_embed"]] <- CreateDimReducObject(
      embeddings = matrix(
        rnorm(2 * ncol(object)),
        ncol = 2,
        dimnames = list(colnames(object), c("AZ_1", "AZ_2"))
      ),
      key = "AZ_",
      assay = "RNA"
    )
    object
  }
  testthat::local_mocked_bindings(
    CloudAzimuth = fake_cloud_azimuth,
    .package = "AzimuthAPI"
  )

  annotated <- annotate_pan_human_azimuth(query, verbose = FALSE)

  expect_equal(sum(is.na(annotated$pan_human_azimuth_broad)), 5)
  expect_false(any(annotated$pan_human_azimuth_fine %in% "False", na.rm = TRUE))
  expect_equal(
    annotated$pan_human_azimuth_label,
    annotated$pan_human_azimuth_final_level_labels
  )
  expect_false("azimuth_embed" %in% Reductions(annotated))
  expect_equal(
    setdiff(colnames(annotated@meta.data), colnames(query@meta.data)),
    c(
      "pan_human_azimuth_full_hierarchical_labels",
      "pan_human_azimuth_final_level_labels",
      "pan_human_azimuth_final_level_confidence",
      "pan_human_azimuth_full_consistent_hierarchy",
      "pan_human_azimuth_broad",
      "pan_human_azimuth_medium",
      "pan_human_azimuth_fine",
      "pan_human_azimuth_label"
    )
  )

  predictions <- readr::read_csv(
    "annotation/annotation-pan-human-azimuth.csv.gz",
    show_col_types = FALSE
  )
  expect_equal(
    colnames(predictions),
    c(
      "cell",
      "full_hierarchical_labels",
      "final_level_labels",
      "final_level_confidence",
      "full_consistent_hierarchy",
      "azimuth_broad",
      "azimuth_medium",
      "azimuth_fine",
      "azimuth_label"
    )
  )
  expect_equal(nrow(predictions), ncol(query))
  expect_equal(sum(predictions$azimuth_broad %in% "False"), 5)
  expect_equal(
    predictions$azimuth_label,
    unname(annotated$pan_human_azimuth_label)
  )
  expect_true(all(
    file.exists(paste0(
      "dr-umap-",
      c(
        "pan_human_azimuth_label",
        "pan_human_azimuth_broad",
        "pan_human_azimuth_medium",
        "pan_human_azimuth_fine",
        "pan_human_azimuth_final_level_labels",
        "pan_human_azimuth_final_level_confidence",
        "pan_human_azimuth_full_consistent_hierarchy",
        "pan_human_azimuth_full_hierarchical_labels"
      ),
      ".png"
    ))
  ))
})

test_that("annotate_pan_human_azimuth() errors before upload without a umap", {
  skip_if_not_installed("AzimuthAPI")

  expect_error(
    annotate_pan_human_azimuth(SeuratObject::pbmc_small),
    "umap"
  )
})

test_that("annotate_pan_human_azimuth() errors before upload if already run", {
  skip_if_not_installed("AzimuthAPI")

  query <- SeuratObject::pbmc_small
  query[["umap"]] <- CreateDimReducObject(
    embeddings = matrix(
      rnorm(2 * ncol(query)),
      ncol = 2,
      dimnames = list(colnames(query), c("UMAP_1", "UMAP_2"))
    ),
    key = "UMAP_",
    assay = "RNA"
  )
  query$pan_human_azimuth_label <- "A"

  expect_error(annotate_pan_human_azimuth(query), "pan_human_azimuth_label")
})
