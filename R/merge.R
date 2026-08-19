#' Merge multiple Seurat objects into one.
#'
#' Merges a list of Seurat objects, rejoins the split layers that Seurat 5
#' leaves behind, and removes genes that are detected in too few cells.
#' Features of the non-RNA assays (such as ADT and HTO) are always kept.
#'
#' Elements of `seurat_object_list` may be Seurat objects, or paths to an
#' analysis directory containing a `seurat_obj.qs2` (as written by the
#' CLI). Which one it is does not matter to the caller: this function
#' resolves each element to an object first.
#'
#' Once every input is an object, this function strips each one's own
#' reductions, scaled data, and clustering columns before merging. A
#' per-sample PCA, tSNE/UMAP, or clustering — written by
#' [create_seurat_object()]/[cluster_seurat_object()] — is meaningless once
#' merged. [normalize_counts()]/[run_pca()] recompute all of it downstream
#' anyway.
#'
#' Writes `metrics-distribution.png`, a per-sample QC violin plot, to the
#' working directory.
#'
#' @param seurat_object_list List of Seurat objects and/or analysis directory
#'   paths (mixing both is fine).
#' @param min_cells Keep genes detected in at least this many cells. If `NULL`,
#'   10 is used, or 0.1% of the cells when there are more than 50,000 of them.
#' @param log_file Filename for the log file.
#'
#' @return A merged Seurat object.
#'
#' @importFrom Matrix rowSums colSums
#' @export
merge_seurat_objects <- function(
  seurat_object_list,
  min_cells = NULL,
  log_file = NULL
) {
  if (length(seurat_object_list) < 2) {
    stop("must have at least 2 samples to merge")
  }

  # Each element may already be a Seurat object, or a path to an analysis
  # directory holding a seurat_obj.qs2 (as written by the CLI).
  # resolve_seurat_object() finds a directory's one qs2/rds file itself, so
  # elements pass straight through.
  #
  # Logging stays here, not inside the utility: only a wrapper like this
  # one knows what is worth logging.
  seurat_object_list <- lapply(seurat_object_list, function(elem) {
    if (is.character(elem)) {
      write_message(
        glue("loading seurat object: {find_seurat_file(elem)}"),
        log_file
      )
    }
    obj <- resolve_seurat_object(elem)
    write_message(glue("cells: {ncol(obj)}"), log_file)
    write_message(glue("assays: {toString(Assays(obj))}"), log_file)
    obj
  })

  # Each object may carry its own per-sample PCA/tSNE/UMAP/clustering. None
  # of that is meaningful once merged, so this strips it before merging,
  # rather than let it leak into the combined object.
  #
  # Variable features are not cleared: Assay5 keeps them in the assay
  # meta.data. normalize_counts() recomputes them downstream.
  seurat_object_list <- lapply(seurat_object_list, function(so) {
    so <- DietSeurat(
      so,
      layers = c("counts", "data"),
      dimreducs = NULL,
      graphs = NULL
    )
    so@meta.data <- so@meta.data |>
      select(-starts_with("snn_res."), -starts_with("res."))
    so
  })

  # Seurat 5 keeps the per-sample layers split after a merge (counts.1,
  # counts.2, ...), which makes GetAssayData() error out, so recombine them
  # into single layers. Nothing downstream consumes split layers here: the
  # batch split for integration is redone by split_layers_by_batch().
  merged_obj <- merge(
    seurat_object_list[[1]],
    seurat_object_list[2:length(seurat_object_list)]
  )
  merged_obj <- SeuratObject::JoinLayers(merged_obj)
  rm(seurat_object_list)

  write_message(glue("merged input cells: {ncol(merged_obj)}"), log_file)
  write_message(glue("merged input genes: {nrow(merged_obj)}"), log_file)

  # filter poorly expressed genes
  detected <- Matrix::rowSums(
    GetAssayData(merged_obj, assay = "RNA", layer = "counts") > 0
  )
  if (!length(min_cells)) {
    min_cells <- 10
    if (ncol(merged_obj) > 50000) min_cells <- ncol(merged_obj) * 0.001
  }
  filtered_genes <- detected[detected >= min_cells] |> names() |> sort()

  # keep all HTO and ADT features if present
  for (extra_assay in c("HTO", "ADT")) {
    if (extra_assay %in% names(merged_obj@assays)) {
      filtered_genes <- c(
        filtered_genes,
        rownames(merged_obj@assays[[extra_assay]])
      )
    }
  }
  merged_obj <- subset(merged_obj, features = filtered_genes)

  # encode sample name as factor (also sets alphabetical sample order)
  merged_obj@meta.data$orig.ident <- factor(merged_obj@meta.data$orig.ident)
  merged_obj <- set_identity(x = merged_obj, identity_column = "orig.ident")

  write_message(glue("merged cells: {ncol(merged_obj)}"), log_file)
  write_message(glue("merged genes: {nrow(merged_obj)}"), log_file)

  # per-sample minimum diagnostics, useful for spotting a single thin sample
  # dragging the gene/cell filter above down
  counts <- GetAssayData(merged_obj, assay = "RNA", layer = "counts")
  write_message(
    glue("min cells per gene: {min(Matrix::rowSums(counts > 0))}"),
    log_file
  )
  write_message(
    glue("min genes per cell: {min(Matrix::colSums(counts > 0))}"),
    log_file
  )
  write_message(
    glue("min counts per cell: {round(min(Matrix::colSums(counts)), 3)}"),
    log_file
  )

  # a merge always combines 2+ samples, so pct_ribo/pct_hb are left out - not
  # worth the width
  plot_metrics_distribution(
    merged_obj,
    metrics = c("detected_genes", "total_counts", "pct_mito"),
    file = "metrics-distribution.png",
    width = metrics_plot_width(merged_obj)
  )

  return(merged_obj)
}
