# Wrappers for pipeline

#' Create a usable Seurat object from a sample.
#'
#' Everything needed to take a sample from a counts file to an object ready
#' for clustering, in one call.
#'
#' The function runs these steps:
#'
#' 1. Read the counts ([read_counts_file()]).
#' 2. Build the object ([initialize_seurat_object()]).
#' 3. Attach any antibody capture data as an "ADT" assay
#'    ([add_seurat_assay()]).
#' 4. Apply the quality cutoffs ([filter_cells()]).
#' 5. Normalize, select variable features, and scale
#'    ([normalize_counts()], which does all three for either method).
#' 6. Run PCA ([run_pca()]).
#' 7. Compute a first-pass tSNE/UMAP ([run_tsne()], [run_umap()]).
#'
#' This function attaches the ADT assay *before* filtering, on purpose.
#' [add_seurat_assay()] restricts the object to the barcodes both matrices
#' share. Attaching ADT after filtering would apply the cutoffs to a
#' different set of cells.
#'
#' PCA, and the tSNE/UMAP computed from it, runs on whichever assay
#' normalization left active. `SCTransform()` creates and activates "SCT".
#' The log path leaves "RNA" active. So `dr_assay` is decided directly from
#' `normalization_method`, not read back off the object. Reading it back off
#' `DefaultAssay()` would reintroduce the implicit dependency an explicit
#' `assay` argument exists to avoid.
#'
#' The tSNE/UMAP computed here are a first-pass look, capped at `num_dim`
#' dimensions. [cluster_seurat_object()] recomputes both once a clustering
#' resolution's own `num_dim` is chosen. So this pair is not the final one
#' used for the actual clustering plots.
#'
#' @param path Path to a 10x directory or a flat counts file.
#' @param sample_name Sample/library name, used as the cell ID prefix.
#' @param num_mads,min_genes,max_genes,min_counts,max_counts,max_mt Quality
#'   cutoffs, passed to [filter_cells()] — see there for how the outlier stage
#'   and the fixed cutoffs combine.
#' @param normalization_method Normalization method, passed to
#'   [normalize_counts()] (as its own `method` argument): "log" or "sct".
#' @param num_variable_genes Number of variable features.
#' @param num_pcs Principal components to compute. Reduced automatically when
#'   the object is too small to support that many.
#' @param num_dim Number of dimensions for the first-pass tSNE/UMAP, capped at
#'   `num_pcs`.
#' @param log_file Filename for the log file.
#'
#' @return A Seurat object with normalized data, variable features, scaled
#'   data, and "pca", "tsne", and "umap" reductions — ready for
#'   [cluster_seurat_object()]. Also writes the variable-feature, PCA, tSNE,
#'   and UMAP diagnostic plots (via [normalize_counts()], [run_pca()],
#'   [run_tsne()], and [run_umap()]) and the metadata+embeddings table (via
#'   [save_metadata()], as `metadata.csv.gz`) to the working directory.
#'
#' @export
create_seurat_object <- function(
  path,
  sample_name,
  num_mads = 3,
  min_genes = 500,
  max_genes = NULL,
  min_counts = NULL,
  max_counts = NULL,
  max_mt = 10,
  normalization_method = "log",
  num_variable_genes = 3000,
  num_pcs = 50,
  num_dim = 30,
  log_file = NULL
) {
  counts_list <- read_counts_file(
    path = path,
    sample_name = sample_name,
    log_file = log_file
  )

  # a flat file with few features is inferred to be antibody capture, so gene
  # expression is not guaranteed to be there; say so rather than let
  # initialize_seurat_object() choke on a NULL
  if (!"Gene Expression" %in% names(counts_list)) {
    stop(glue(
      "no gene expression data in {path}, found: ",
      "{paste(names(counts_list), collapse = ', ')}"
    ))
  }

  so <- initialize_seurat_object(
    counts_matrix = counts_list[["Gene Expression"]],
    assay = "RNA",
    log_file = log_file
  )

  if ("Antibody Capture" %in% names(counts_list)) {
    so <- add_seurat_assay(
      so,
      assay = "ADT",
      counts_matrix = counts_list[["Antibody Capture"]],
      log_file = log_file
    )
  }

  message("\n\n ===== filter_cells() ===== \n\n")
  so <- filter_cells(
    so,
    num_mads = num_mads,
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    max_mt = max_mt,
    log_file = log_file
  )

  # normalizes, selects variable features, and scales in one call, for either
  # method
  message(glue(
    "\n\n ===== normalize_counts() ({normalization_method}) ===== \n\n"
  ))
  so <- normalize_counts(
    so,
    method = normalization_method,
    num_variable_genes = num_variable_genes,
    assay = "RNA",
    log_file = log_file
  )

  # SCTransform() creates and activates "SCT"; the log path leaves "RNA" active
  # - decided from `normalization_method` directly, not read back off
  # DefaultAssay()
  dr_assay <- if (normalization_method == "sct") "SCT" else "RNA"

  # RunPCA() needs strictly fewer components than either dimension of the
  # scaled matrix
  num_pcs <- min(num_pcs, ncol(so) - 1, nrow(so) - 1)

  so <- run_pca(
    so,
    num_pcs = num_pcs,
    assay = dr_assay,
    var_features = TRUE
  )

  # rough first-pass tSNE/UMAP for a quick look; cluster_seurat_object()
  # computes the real ones once a clustering resolution's num_dim is chosen
  num_dim <- min(as.integer(num_dim), num_pcs)
  so <- run_tsne(
    so,
    reduction = "pca",
    num_dim = num_dim,
    assay = dr_assay,
    file_format = "png"
  )
  so <- run_umap(
    so,
    reduction = "pca",
    num_dim = num_dim,
    assay = dr_assay,
    file_format = "png"
  )

  save_metadata(so)

  so
}
