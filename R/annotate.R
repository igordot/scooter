#' Transfer labels from a reference Seurat object
#'
#' Reference-based annotation wrapping Seurat's reference-mapping workflow
#' ([Seurat::FindTransferAnchors()] + [Seurat::TransferData()]) to annotate
#' `query` with labels (cell type, or any other categorical metadata) from an
#' already integrated and annotated reference object. Optionally also projects
#' `query` onto the reference's UMAP via [Seurat::MapQuery()].
#'
#' @param query Query Seurat object.
#' @param ref Reference Seurat object - already integrated and annotated.
#' @param query_label_col Metadata column name added to `query` for the
#'   transferred labels.
#' @param ref_label_col Metadata column in `ref` holding the labels to
#'   transfer (cell type, or any other category).
#' @param query_assay Assay in `query` to use for anchor finding. `"RNA"`
#'   (log-normalized) or `"SCT"`, depending on how `query` was built. No
#'   default, like the params below: this and `ref_assay`/`normalization`/
#'   `num_dim` all describe how a specific object was built.
#' @param ref_assay Assay in `ref` to use for anchor finding.
#' @param normalization Normalization method shared by `ref` and `query`
#'   (`"log"`/`"LogNormalize"` or `"sct"`/`"SCT"`).
#'   [Seurat::FindTransferAnchors()] takes one value for both objects, so
#'   `ref` and `query` must actually share a method.
#' @param ref_reduction Reduction in `ref` to project `query` onto.
#' @param num_dim Number of dimensions to use for anchor finding and label
#'   transfer.
#' @param anchor_reduction Dimensional reduction workflow
#'   [Seurat::FindTransferAnchors()] uses to find anchors.
#' @param ref_umap If `TRUE`, also runs [Seurat::MapQuery()] and adds a
#'   `"ref.umap"` reduction to `query`, projecting it onto `ref`'s UMAP. `ref`
#'   needs a UMAP built with `return.model = TRUE`, which this function does
#'   not check for - `MapQuery()` errors on its own if it is missing.
#' @param verbose Print progress messages.
#'
#' @return `query` with `query_label_col` (a factor) added to its metadata,
#'   and a `"ref.umap"` reduction added if `ref_umap = TRUE`. Also writes
#'   the predictions table (`cell`/`predicted.id`/`prediction.score.max`) to
#'   `annotation/annotation-transfer-label-<query_label_col>.csv.gz`, and a
#'   `dr-<reduction>-<query_label_col>.png` UMAP colored by the transferred
#'   labels for each of `"umap"`/`"ref.umap"` already present on `query`.
#'
#' @importFrom readr write_csv
#' @export
transfer_labels <- function(
  query,
  ref,
  query_label_col,
  ref_label_col,
  query_assay,
  ref_assay,
  normalization,
  ref_reduction = "pca",
  num_dim,
  anchor_reduction = "pcaproject",
  ref_umap = FALSE,
  verbose = TRUE
) {
  if (!ref_label_col %in% colnames(ref@meta.data)) {
    stop(glue(
      "ref_label_col '{ref_label_col}' not found; available: ",
      "{toString(colnames(ref@meta.data))}"
    ))
  }

  normalization <- switch(
    normalization,
    log = ,
    LogNormalize = "LogNormalize",
    sct = ,
    SCT = "SCT",
    stop(glue(
      "normalization must be one of: log, LogNormalize, sct, SCT"
    ))
  )

  if (verbose) {
    message(glue("transferring labels from reference column: {ref_label_col}"))
    ref_labels <- sort(unique(as.character(ref[[ref_label_col]][[1]])))
    message(glue("reference labels: {toString(ref_labels)}"))
  }

  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    normalization.method = normalization,
    reference.assay = ref_assay,
    query.assay = query_assay,
    dims = 1:num_dim,
    reduction = anchor_reduction,
    reference.reduction = ref_reduction
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref@meta.data[[ref_label_col]],
    dims = 1:num_dim
  )

  if (verbose) {
    message(glue(
      "predicted labels:\n",
      "{paste(capture.output(table(predictions$predicted.id)), ",
      "collapse = '\n')}"
    ))
    message(glue(
      "cells annotated: {nrow(predictions)}, ",
      "distinct labels: {n_distinct(predictions$predicted.id)}"
    ))
  }

  predictions_tbl <- as_tibble(predictions, rownames = "cell")
  predictions_tbl <- select(
    predictions_tbl,
    .data$cell,
    .data$predicted.id,
    .data$prediction.score.max
  )

  annotation_dir <- "annotation"
  if (!dir.exists(annotation_dir)) {
    dir.create(annotation_dir)
  }
  write_csv(
    predictions_tbl,
    glue("{annotation_dir}/annotation-transfer-label-{query_label_col}.csv.gz")
  )

  new_metadata <- setNames(
    data.frame(
      as.factor(predictions_tbl$predicted.id),
      row.names = predictions_tbl$cell
    ),
    query_label_col
  )
  query <- AddMetaData(query, metadata = new_metadata)

  if (ref_umap) {
    query <- MapQuery(
      anchorset = anchors,
      reference = ref,
      query = query,
      reference.reduction = ref_reduction,
      reduction.model = "umap"
    )
  }

  umap_reductions <- intersect(c("umap", "ref.umap"), Reductions(query))
  for (umap_reduction in umap_reductions) {
    plot_dr_group(
      query,
      reduction = umap_reduction,
      group_by = query_label_col,
      file_prefix = glue("dr-{umap_reduction}-{query_label_col}"),
      file_format = "png"
    )
  }

  query
}
