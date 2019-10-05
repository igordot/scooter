#' Filter out cells based on minimum and maximum number of genes and maximum percentage
#' mitochondrial reads. If cutoffs are not provided, the min_genes will be the 0.02 quantile, and the max genes will be 0.98 quantile and the mitochondrial percentage will be 10%.
#'
#' @param data A tibble with metadata.
#' @param min_genes Minimum number of genes per cell.
#' @param max_genes Maximim number of genes per cell.
#' @param max_mt Maximum percentage of mitochondrial reads per cell.
#' @param log_file log file.
#'
#' @return Filtered data
#'
#' @import dplyr
#' @importFrom stats quantile
#' @export
filter_data <- function(data, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10) {
  UseMethod("filter_data")
}

#' @return cells to keep
#' @export
filter_data.default <- function(data, log_file = NULL, min_genes, max_genes, max_mt) {

  cells_subset =
    data %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    filter(.data$nFeature_RNA > min_genes,
           .data$nFeature_RNA < max_genes,
           .data$percent.mito < max_mt) %>%
    pull(.data$cell)

  return(cells_subset)
}

#' @export
filter_data.Seurat <- function(data, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10) {

  s_obj = data
  message_str <- glue("\n\n ========== filter data matrix ========== \n\n
                      unfiltered min genes: {min(s_obj$nFeature_RNA)}
                      unfiltered max genes: {max(s_obj$nFeature_RNA)}
                      unfiltered mean num genes: {round(mean(s_obj$nFeature_RNA), 3)}
                      unfiltered median num genes: {median(s_obj$nFeature_RNA)}")
  write_message(message_str, log_file)

  # convert arguments to integers (command line arguments end up as characters)
  min_genes = as.numeric(min_genes)
  max_genes = as.numeric(max_genes)
  max_mt = as.numeric(max_mt)

  # default cutoffs (gene numbers rounded to nearest 10)
  # as.numeric() converts NULLs to 0 length numerics, so can't use is.null()
  if (!length(min_genes)) min_genes = s_obj$nFeature_RNA %>%
    quantile(0.02, names = FALSE) %>%
    round(-1)

  if (!length(max_genes)) max_genes = s_obj$nFeature_RNA %>%
    quantile(0.98, names = FALSE) %>%
    round(-1)

  if (!length(max_mt)) max_mt = 10

  message_str <- glue("min genes cutoff: {min_genes}
                      max genes cutoff: {max_genes}
                      max mitochondrial percentage cutoff: {max_mt}
                      imported cells: {ncol(s_obj)}
                      imported genes: {nrow(s_obj)}")
  write_message(message_str, log_file)

  # filter
  cells_subset <- filter_data(s_obj@meta.data,
                              log_file,
                              min_genes,
                              max_genes,
                              max_mt)
  # subset based on the cells that passed the filtering
  s_obj = subset(s_obj, cells = cells_subset)

  message_str <- glue("filtered cells: {ncol(s_obj)}
                      filtered genes: {nrow(s_obj)}")
  write_message(message_str, log_file)

  return(s_obj)
}

#' Normalize data
#'
#' @param data Input data.
#' @param normalize_method .
#' @param nfeatures .
#' @param metadata .
#' @param assay .
#' @param log_file log file.
#'
#' @return Normalized data
#'
#' @import dplyr
#' @export
normalize_data <- function(data, normalize_method, nfeatures = 2000, metadata = NULL, assay = NULL, log_file = NULL) {
  UseMethod("normalize_data")
}

normalize_data.default <- function(data, normalize_method, nfeatures = 2000, metadata = NULL, assay = NULL, log_file = NULL) {

  if (normalize_method == "log_norm") {

    normed <- log_normalize_data(data = data,
                                 log_file = log_file)

  } else if ( normalize_method == "sct") {

    normed <- sctransform_data(counts = data,
                               metadata = metadata,
                               nfeatures =  nfeatures,
                               log_file = log_file)
  }
  return(normed)
}

normalize_data.Seurat <- function(data, normalize_method, nfeatures = 2000, metadata = NULL, assay = "RNA", log_file = NULL) {

  if (normalize_method == "sct") {

    # get counts data for the specified assay
    assay.obj <- GetAssay(object = data, assay = assay)
    counts <- GetAssayData(object = assay.obj, slot = 'counts')

    # get metadata
    metadata <- data@meta.data

    # run normalization
    normed_data <- normalize_data(data = counts,
                                  normalize_method = normalize_method,
                                  nfeatures = nfeatures,
                                  metadata = metadata,
                                  log_file = log_file,
                                  assay = assay)

    # Create new SCT assay.

    # set the counts of the SCT assay to the output of the sct
    assay.out <- CreateAssayObject(counts = normed_data[["vst.out"]]$umi_corrected)

    # set the variable features to the top variable features from sct
    VariableFeatures(object = assay.out) <- normed_data[["top.features"]]

    # log norm the counts data from sct to make the data data
    assay.out <- SetAssayData(
      object = assay.out,
      slot = 'data',
      new.data = log1p(GetAssayData(object = assay.out, slot = 'counts')))

    # set the scale data to the scaled sct'ed data
    assay.out <- SetAssayData(
      object = assay.out,
      slot = 'scale.data',
      new.data = normed_data[["scale.data"]]
    )

    data[["SCT"]] <- assay.out

  } else {

    # get counts for specified assay
    assay.obj <- GetAssay(object = data, assay = assay)
    counts <- GetAssayData(object = assay.obj, slot = 'counts')

    # normalize data
    normed_data <- normalize_data(counts, normalize_method)

    # set data slot to normalized data
    data[[assay]]@data <- normed_data

    # calculate variance of tje assay and get most variable genes
    scaled.data <- calculate_variance(seurat_obj = data,
                                      assay = assay,
                                      nfeatures = nfeatures,
                                      log_file = NULL)

    # add scaled data to scale data slot
    data[[assay]]@scale.data <- scaled.data$scaled.data

    # add top features to variable features
    VariableFeatures(data[[assay]]) <- scaled.data$top.features
  }
  return(data)
}

#' Log normalize data.
#'
#' @param data A seurat object.
#' @param log_file log file.
#'
#' @return normalized data
#'
#' @importFrom Seurat NormalizeData
#' @export
log_normalize_data <- function(data, log_file = NULL) {
  # log normalize data

  message_str <- "\n\n ========== log normalize ========== \n\n"
  write_message(message_str, log_file)

  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  data = NormalizeData(data,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000,
                       verbose = FALSE)

  return(data)
}

#' SCT normalize data.
#'
#' @param counts Raw counts.
#' @param metadata Metadata for each cell.
#' @param nfeatures Number of variable features to output.
#' @param log_file A log file.
#'
#' @return A named list of the vst output, the final scaled data, and the top variable genes.
#'
#' @import dplyr
#' @importFrom sctransform vst
#' @export
sctransform_data <- function(counts, metadata, nfeatures, log_file = NULL){
  # sc transform data

  message_str = "\n\n ========== sc transform ========== \n\n"
  write_message(message_str, log_file)

  # define clip range. This is unchanged from Seurat
  clip.range <- c(-sqrt(ncol(counts)),sqrt(ncol(counts)))

  vst.out <- vst(counts,
                 cell_attr = metadata,
                 latent_var = "log_umi",
                 show_progress = TRUE,
                 return_cell_attr = TRUE,
                 return_gene_attr = TRUE,
                 return_corrected_umi = TRUE,
                 residual_type = "pearson",
                 res_clip_range = clip.range)

  # name the residual variance for each gene
  feature.variance <- setNames(object = vst.out$gene_attr$residual_variance,
                               nm = rownames(vst.out$gene_attr))

  # sort the residual variance for each gene
  feature.variance <- sort(feature.variance,
                           decreasing = TRUE)

  # get the genes with the top variance
  top.features <- names(feature.variance)[1:min(nfeatures,
                                                length(feature.variance))]
  # get and clean sc transformed data
  scale.data <- vst.out$y

  scale.data[scale.data < clip.range[1]] <- clip.range[1]
  scale.data[scale.data > clip.range[2]] <- clip.range[2]

  # regress out percent mitochondria and ncount_RNA
  scale.data <- ScaleData(
    scale.data,
    features = NULL,
    vars.to.regress = c("percent.mito", "nCount_RNA"),
    latent.data = metadata[, c("percent.mito", "nCount_RNA"), drop = FALSE],
    model.use = 'linear',
    use.umi = FALSE,
    do.scale = FALSE,
    do.center = TRUE,
    scale.max = Inf,
    block.size = 750,
    min.cells.to.block = 3000,
    verbose = TRUE
  )

  return(list(vst.out = vst.out,
              scale.data = scale.data,
              top.features = top.features))

}

#' Get variable genes and scale data.
#'
#' @param seurat_obj Seurat object.
#' @param assay Assay.
#' @param nfeatures Number of variable features to output.
#' @param log_file A log file.
#'
#' @return A named list of the top features, and the scaled data.
#'
#' @import dplyr
#' @importFrom Seurat FindVariableFeatures ScaleData
#' @export
calculate_variance <- function(seurat_obj, assay = "RNA", nfeatures = 2000, log_file = NULL){
  # calculate variance of genes in a seurat object

  s_obj = seurat_obj

  message_str <- "\n\n ========== Seurat::FindVariableGenes() ========== \n\n"
  write_message(message_str, log_file)

  # identify features that are outliers on a 'mean variability plot'
  var_features = FindVariableFeatures(s_obj[[assay]],
                                      selection.method = "vst",
                                      nfeatures = nfeatures,
                                      verbose = FALSE)

  message_str <- "\n\n ========== Seurat::ScaleData() ========== \n\n"
  write_message(message_str, log_file)

  # regress out unwanted sources of variation
  # regressing uninteresting sources of variation can improve dimensionality reduction and clustering
  # could include technical noise, batch effects, biological sources of variation (cell cycle stage)
  # scaled z-scored residuals of these models are stored in scale.data slot
  # used for dimensionality reduction and clustering
  # RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData
  # s_obj = ScaleData(s_obj, features = rownames(s_obj), vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)

  metadata <- s_obj@meta.data
  # scale data in the assay
  scaled.data = ScaleData(s_obj[[assay]],
                          vars.to.regress = c("percent.mito", "nCount_RNA"),
                          latent.data = metadata[,c("percent.mito", "nCount_RNA")],
                          verbose = FALSE)

  return(list(top.features = var_features@var.features,
              scaled.data = scaled.data@scale.data))
}

