#' Split a Seurat object into per-batch layers and prepare it for integration.
#'
#' Drops the scaled data, reductions, graphs, and clustering columns, then
#' splits the assay layers by a metadata column — splitting the layers rather
#' than the object is what the Seurat 5 `IntegrateLayers()` workflow expects.
#'
#' Everything after the split runs through the same verbs as the single-sample
#' path ([normalize_counts()], [run_pca()]), so the two workflows cannot drift
#' apart on parameters. [integrate_layers()] then consumes the resulting "pca"
#' reduction.
#'
#' @param x Seurat object.
#' @param batch_var Metadata column to split the layers by.
#' @param assay Assay to split.
#' @param num_variable_genes Number of variable features. Matches the
#'   single-sample default.
#' @param num_dim Number of principal components to compute.
#' @param log_file Filename for the log file.
#'
#' @return A Seurat object with split layers and a "pca" reduction.
#'
#' @export
split_layers_by_batch <- function(
  x,
  batch_var = "orig.ident",
  assay = "RNA",
  num_variable_genes = 3000,
  num_dim = 50,
  log_file = NULL
) {
  if (!batch_var %in% colnames(x@meta.data)) {
    stop(glue("batch variable {batch_var} is not in the metadata"))
  }

  # command line arguments end up as characters
  num_dim <- as.integer(num_dim)
  num_variable_genes <- as.integer(num_variable_genes)

  # clean up object (drops scale.data, reductions, and graphs)
  x <- DietSeurat(
    x,
    layers = c("counts", "data"),
    dimreducs = NULL,
    graphs = NULL
  )
  x@meta.data <- x@meta.data |> select(-starts_with("snn_res."))
  x@meta.data <- x@meta.data |> select(-starts_with("res."))

  batch_sizes <- table(x@meta.data[[batch_var]])
  for (batch_name in names(batch_sizes)) {
    write_message(
      glue("batch {batch_name} cells: {batch_sizes[[batch_name]]}"),
      log_file
    )
  }

  # split the layers by batch, so every step below is computed per batch
  x[[assay]] <- split(x[[assay]], f = x@meta.data[[batch_var]])

  # normalize_counts() now normalizes, selects variable features, and scales in
  # one call (matching the single-sample path exactly, including
  # vars.to.regress = NULL), so a second, separate variance step is no longer
  # needed here
  x <- normalize_counts(
    x,
    method = "log",
    num_variable_genes = num_variable_genes,
    assay = assay,
    log_file = log_file
  )

  # suffix = "" keeps the reduction named "pca", which is what
  # IntegrateLayers() looks for
  x <- run_pca(x, num_pcs = num_dim, assay = assay, var_features = TRUE)

  return(x)
}
