#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param data Data.
#' @param assay Assay in Seurat object.
#' @param var_features Boolean. To use variable features.
#' @param num_pcs Number of PCs to calculate.
#' @param num_dim Number of PCs to use for tsne and umap.
#' @param num_neighbors Number of neighbors to use for umap.
#' @param log_file log file.
#' @param prefix prefix for dimensionality reduction labels.
#'
#' @return list of dimensionality reductions.
#'
#' @import dplyr
#' @import Seurat
#' @export
run_dimensionality_reduction <- function(data, assay = "RNA", var_features = TRUE, num_pcs, num_dim, num_neighbors, log_file = NULL, prefix = NULL){

  message_str <- "\n\n ========== dimensionality reduction ========== \n\n"
  write_message(message_str, log_file)

  prefix <- ifelse(is.null(prefix), "", prefix)

  pca_out <- run_dr(data = data, dr_method = "pca",
                    var_features = var_features,
                    assay = assay,
                    num_pcs = num_pcs,
                    prefix = paste0("PC", prefix, "_"))

  feature.loadings <- pca_out$feature.loadings
  cell.embeddings <- pca_out$cell.embeddings
  sdev = pca_out$sdev
  pca_obj = pca_out$pca_out

  tsne_out <- run_dr(data = cell.embeddings[,1:num_dim], dr_method = "tsne",
                      prefix = paste0("tSNE", prefix, "_"))

  umap_out<- run_dr(data = cell.embeddings[, 1:num_dim], dr_method = "umap",
                    num_neighbors = num_neighbors, prefix = paste0("UMAP", prefix, "_"))

  return(list(feature.loadings = feature.loadings,
              cell.embeddings = cell.embeddings,
              sdev = sdev,
              pca_obj = pca_obj,
              tsne_out = tsne_out,
              umap_out = umap_out))
}

add_dim_red_seurat <- function(seurat_obj, dim_red_list, prefix = NULL){

  pca.dim.reduc <- new(Class = "DimReduc",
                       cell.embeddings =  as.matrix(dim_red_list$cell.embeddings),
                       feature.loadings = as.matrix(dim_red_list$feature.loadings),
                       assay.used = "RNA",
                       stdev = dim_red_list$sdev,
                       key = paste0("PC", prefix))

  tsne.dim.reduc <- new(Class = "DimReduc",
                        cell.embeddings =  as.matrix(dim_red_list$tsne_out),
                        assay.used = "RNA",
                        key = paste0("tSNE", prefix))

  umap.dim.reduc <- new(Class = "DimReduc",
                        cell.embeddings =  as.matrix(dim_red_list$umap_out),
                        assay.used = "RNA",
                        key = paste0("UMAP", prefix))


  prefix <- ifelse(is.null(prefix), "", prefix)

  seurat_obj[[paste("pca", prefix, sep = "")]] <- pca.dim.reduc
  seurat_obj[[paste("tsne", prefix, sep = "")]] <- tsne.dim.reduc
  seurat_obj[[paste("umap", prefix, sep = "")]] <- umap.dim.reduc

  return(seurat_obj)
}

#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param data Data to use for dimensionality reduction.
#' @param dr_method Dimensionality reduction method.
#' @param assay .
#' @param features .
#' @param graph .
#' @param num_dim .
#' @param var_features .
#'
#' @return list of dimensionality reduced/
#'
#' @import dplyr
#' @importFrom Seurat Embeddings DefaultAssay GetAssayData VariableFeatures
#' @export
#' @seealso run_pca(), run_tsne(), run_umap()
run_dr <- function(data, dr_method, ...) {
  UseMethod("run_dr")
}

#' @export
run_dr.default <- function(data, dr_method = c("pca", "tsne", "umap"), ...) {

  dr_method = match.arg(dr_method)

  if (dr_method == "pca") {

    out = run_pca(data = data, ...)

  } else if (dr_method == "tsne") {

    out = run_tsne(data = data, ...)

  } else if (dr_method == "umap") {

    out = run_umap(data = data, ...)

  }
}

#' @export
run_dr.Seurat <- function(data, dr_method, assay = NULL, var_features = FALSE, features = NULL, graph = NULL, num_dim_use = NULL, ...) {

  if(var_features){
    features <- VariableFeatures(data[[assay]])
    data.use <- t(GetAssayData(object = data,
                               slot = 'data',
                               assay = assay)[features, ])
  # if the features are specified, get features from the specified assay
  } else if (!is.null(features)) {

    data.use <- t(GetAssayData(object = data,
                                   slot = 'data',
                                   assay = assay)[features, ])
  # if the number of dimensions is specified, use those dimensions from pca
  } else if (!is.null(num_dim_use)) {

    data.use <- Embeddings(data[[reduction]])[, num_dim_use]
    assay <- DefaultAssay(object = data[[reduction]])

  # if a graph is specified, use the graph
  } else if (!is.null(graph)) {
    data.use <- data[[graph]]
  } else {
    stop("Please specify one of num_dim_use, features, or graph")
  }
  # run the specified dimensionality reduction method with the specified data
  run_dr(data = data.use, dr_method = dr_method, ...)
}

#' Run PCA
#'
#' @param data A tibble with metadata.
#' @param num_pcs Maximim number of genes per cell.
#' @param prefix suffix.
#'
#' @return named list of feature loadings, cell embeddings, sdev, output from pca
#'
#' @import irlba
#' @export
run_pca <- function(data, num_pcs, prefix = "PC_"){

  # make sure the number of pcs isn't too large
  npcs <- min(num_pcs, nrow(data))

  # calculate PCA
  pca_out <- irlba::prcomp_irlba(data, n = npcs)

  # Extract information
  feature.loadings <- pca_out$rotation
  sdev <- pca_out$sdev
  total.variance <- sum(sdev)
  cell.embeddings <- pca_out$x

  rownames(feature.loadings) <- colnames(data)
  colnames(feature.loadings) <- paste0(prefix, 1:npcs)
  rownames(cell.embeddings) <- rownames(data)
  colnames(cell.embeddings) <- paste0(prefix, 1:npcs)

  return(list(feature.loadings = feature.loadings,
              cell.embeddings = cell.embeddings,
              sdev = sdev,
              pca_out = pca_out))
}

#' Run TSNE
#'
#' @param data Data to run tsne on.
#' @param seed.use seed to use.
#' @param dim.embed Number of tsne embeddings to return.
#' @param prefix suffix.
#'
#' @return tsne.
#'
#' @import Rtsne
#'
#' @export
run_tsne <- function(data, seed.use = 1,  dim.embed = 2, prefix = "tSNE_"){

  set.seed(seed = seed.use)
  # run tsne
  tsne.data <- Rtsne(data, dims = dim.embed)$Y

  colnames(x = tsne.data) <- paste0(prefix, 1:ncol(x = tsne.data))
  rownames(x = tsne.data) <- rownames(x = data)

  return(tsne.data)
}


#' Run UMAP
#'
#' @param data Data to run UMAP on.
#' @param num_neighbors num neighbors.
#' @param prefix suffix.
#'
#' @return umap
#'
#' @importFrom uwot umap
#' @export
run_umap <- function(data, num_neighbors, prefix = "UMAP_") {

  # run umap
  umap_out <- umap(data, n_neighbors = num_neighbors)

  colnames(umap_out) <- paste0(prefix,
                               1:ncol(umap_out))

  rownames(umap_out) <- rownames(data)

  return(umap_out)
}

