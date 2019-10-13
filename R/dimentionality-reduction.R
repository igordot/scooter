#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param data Data to use for dimensionality reduction.
#' @param dr_method Dimensionality reduction method.
#' @param prefix .
#' @param assay .
#' @param features .
#' @param graph .
#' @param num_dim .
#' @param min_dist .
#' @param var_features .
#'
#' @return list of dimensionality reduced/
#'
#' @import dplyr
#' @importFrom Seurat Embeddings DefaultAssay GetAssayData VariableFeatures
#' @export
#' @seealso run_pca(), run_tsne(), run_umap()
run_dr <- function(data, dr_method, prefix, assay = NULL,
                   var_features = FALSE, features = NULL,
                   graph = NULL, num_dim_use = NULL, reduction = NULL,
                   num_neighbors = NULL, num_pcs= NULL,...) {
  UseMethod("run_dr")
}

#' @export
run_dr.default <- function(data, dr_method = c("pca", "tsne", "umap"), prefix, assay = NULL,
                           var_features = FALSE, features = NULL,
                           graph = NULL, num_dim_use = NULL, reduction = NULL,
                           num_neighbors = NULL, num_pcs= NULL, ...) {

  dr_method = match.arg(dr_method)

  if (dr_method == "pca") {

    out = run_pca(data = data, prefix = prefix, num_pcs = num_pcs)

  } else if (dr_method == "tsne") {

    out = run_tsne(data = data, prefix = prefix)

  } else if (dr_method == "umap") {

    out = run_umap(data = data, prefix = prefix, num_neighbors = num_neighbors)

  }
}

#' @export
run_dr.Seurat <- function(data, dr_method, prefix, assay = NULL,
                          var_features = FALSE, features = NULL,
                          graph = NULL, num_dim_use = NULL, reduction = NULL,
                          num_neighbors = NULL, num_pcs= NULL, ...) {

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

    data.use <- Embeddings(data[[reduction]])[,1:num_dim_use]
    assay <- DefaultAssay(object = data[[reduction]])

  # if a graph is specified, use the graph
  } else if (!is.null(graph)) {
    data.use <- data[[graph]]
  } else {
    stop("Please specify one of num_dim_use, features, or graph")
  }

  # run the specified dimensionality reduction method with the specified data
  dim_reduction <- run_dr(data = data.use, dr_method = dr_method,
                          prefix = prefix, assay = assay,
                          var_features = var_features, features = features,
                          graph = graph, num_dim_use = num_dim_use, reduction = reduction,
                          num_neighbors = num_neighbors, num_pcs= num_pcs, ...)

  if(dr_method == "pca") {
    cell.embeddings <- as.matrix(dim_reduction$cell.embeddings)
    colnames(cell.embeddings) <- paste0("PC", prefix, 1:ncol(cell.embeddings))
    feature.loadings <- as.matrix(dim_reduction$feature.loadings)
    colnames(feature.loadings) <- paste0("PC", prefix, 1:ncol(feature.loadings))
    dim_red_class <- new(Class = "DimReduc",
                         cell.embeddings =  cell.embeddings,
                         feature.loadings = feature.loadings,
                         assay.used = assay,
                         stdev = dim_reduction$sdev,
                         key = paste0("PC", prefix))
  } else if(dr_method == "tsne") {
    cell.embeddings <- as.matrix(dim_reduction)
    colnames(cell.embeddings) <- paste0("tSNE", prefix, 1:ncol(cell.embeddings))
    dim_red_class <- new(Class = "DimReduc",
                          cell.embeddings =  cell.embeddings,
                          assay.used = assay,
                          key = paste0("tSNE", prefix))
  } else if(dr_method == "umap") {
    cell.embeddings <- as.matrix(dim_reduction)
    colnames(cell.embeddings) <- paste0("UMAP", prefix, 1:ncol(cell.embeddings))
    dim_red_class <- new(Class = "DimReduc",
                          cell.embeddings =  cell.embeddings,
                          assay.used = assay,
                          key = paste0("UMAP", prefix))
  }
  data[[glue("{dr_method}{prefix}")]] <- dim_red_class

  return(data)
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
run_umap <- function(data, num_neighbors, min_dist = 0.3, prefix = "UMAP_") {
  # run umap
  umap_out <- umap(data, n_neighbors = num_neighbors, min_dist = min_dist)

  colnames(umap_out) <- paste0(prefix,
                               1:ncol(umap_out))

  rownames(umap_out) <- rownames(data)

  return(umap_out)
}

