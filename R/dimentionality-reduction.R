#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param data A tibble with metadata.
#' @param num_pcs Minimum number of genes per cell.
#' @param num_dim Maximim number of genes per cell.
#' @param num_neighbors Maximum percentage of mitochondrial reads per cell.
#' @param log_file log file.
#' @param prefix suffix.
#'
#' @return list of dimensionality reduced/
#'
#' @import dplyr
#' @import Seurat
#' @export
run_dimensionality_reduction<- function(data, num_pcs, num_dim, num_neighbors, log_file = NULL, prefix = NULL){

  message_str <- "\n\n ========== dimensionality reduction ========== \n\n"
  write_message(message_str, log_file)

  prefix <- ifelse(is.null(prefix), "", prefix)

  pca_out <- run_dr(data = data, dr_method = "pca",
                    num_dim = num_dim, prefix = paste0("PC", prefix, "_"))

  feature.loadings <- pca_out[[1]]
  cell.embeddings <- pca_out[[2]]
  sdev = pca_out[[3]]
  pca_obj = pca_out[[4]]

  tsne_out <- run_dr(data = cell.embeddings, dr_method = "tsne",
                     num_dim = num_dim, prefix = paste0("tSNE", prefix, "_"))

  umap_out<- run_dr(data = cell.embeddings, num_dim = num_dim,
                    num_neighbors = num_neighbors, prefix = paste0("UMAP", prefix, "_"))

  return(list(feature.loadings = feature.loadings,
              cell.embeddings = cell.embeddings,
              sdev = sdev,
              pca_obj = pca_obj,
              tsne_out = tsne_out,
              umap_out = umap_uwot_out))
}

#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param data Data to use for dimensionality reduction.
#' @param dr_method Dimensionality reduction method.
#' @param assay .
#' @param features .
#' @param graph .
#' @param num_dim .
#'
#' @return list of dimensionality reduced/
#'
#' @import dplyr
#' @import Seurat
#' @export
run_dr <- function(data, dr_method, ...) {
  UseMethod("run_dr")
}

#' @export
run_dr.default <- function(data, dr_method, ...) {
  if (dr_method == "pca") {
    out = run_pca(data = data, ...)
  } else if (dr_method == "tsne") {
    out = run_tsne(data = data, ...)
  } else if (dr_method == "umap") {
    out = run_umap(data = data, ...)
  } else {
    stop("Unknown dimensionality reduction method.")
  }
}

#' @export
run_dr.Seurat <- function(data, dr_method, assay, features, graph, num_dim, ...) {

  if (!is.null(x = features)) {
    data.use <- t(x = GetAssayData(object = data,
                                   slot = 'data',
                                   assay = assay)[features, ])
  } else if (!is.null(x = num_dim)) {
    data.use <- Embeddings(data[[reduction]])[, num_dim]
    assay <- DefaultAssay(object = data[[reduction]])
  } else if (!is.null(x = graph)) {
    data.use <- data[[graph]]
  } else {
    stop("Please specify one of num_dim, features, or graph")
  }
  run_dr(data = data.use, dr_method = dr_method, ...)
}

#' Run PCA
#'
#' @param data A tibble with metadata.
#' @param num_dim Maximim number of genes per cell.
#' @param prefix suffix.
#'
#' @return pca.
#'
#' @import irlba
#' @export

run_pca <- function(data, num_dim, prefix = "PC_"){

  npcs <- min(num_dim, nrow(data))
  pca_out <- irlba::prcomp_irlba(t(data), n = npcs)

  feature.loadings <- pca_out$rotation
  sdev <- pca_out$sdev
  total.variance <- sum(sdev)
  cell.embeddings <- pca_out$x

  rownames(x = feature.loadings) <- rownames(data)
  colnames(x = feature.loadings) <- paste0(prefix, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(data)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)

  return(list(feature.loadings, cell.embeddings, sdev, pca_out))
}

#' Run TSNE
#'
#' @param data A tibble with metadata.
#' @param seed.use .
#' @param dim.embed .
#' @param prefix suffix.
#'
#' @return pca.
#'
#' @import irlba
#'
#' @export
run_tsne <- function(data, seed.use = 1,  dim.embed = 2, prefix = "tSNE_"){

  set.seed(seed = seed.use)
  tsne.data <- Rtsne(data, dims = dim.embed)$Y
  colnames(x = tsne.data) <- paste0(prefix, 1:ncol(x = tsne.data))
  rownames(x = tsne.data) <- rownames(x = data)

  return(tsne.data)
}


#' Run UMAP
#'
#' @param data A tibble with metadata.
#' @param num_neighbors num neighbors.
#' @param prefix suffix.
#'
#' @return umap
#'
#' @importFrom uwot umap
#' @export
run_umap <- function(data, num_neighbors, prefix = "UMAP_") {

  umap_out <- umap(data, n_neighbors = num_neighbors)

  colnames(umap_out) <- paste0(prefix,
                               dim_red_suffix,
                               1:ncol(umap_out))

  rownames(umap_out) <- rownames(cell.embeddings)

  return(umap_out)
}

