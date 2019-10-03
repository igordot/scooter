#' Read in 10x data.
#'
#' @param data_path Path to directory that holds the 10x output files.
#' @param gene.column The column with the gene names.
#'
#' @return list of matrices.
#'
#' @importFrom Matrix readMM
#' @export
import_matrix <- function(data_path, gene.column = 2){
  # Heavily sourced from Seurat

  # check if the directory exists
  if (!dir.exists(paths = data_path)) {
    stop(glue("dir {data_path} does not exist"))
  }

  # read in features.tsv.gz
  feature.names <- read.delim(
    file = file.path(data_path, 'features.tsv.gz'),
    header = FALSE,
    stringsAsFactors = FALSE)

  # read in sparse matrix
  data <- readMM(file = file.path(data_path, 'matrix.mtx.gz'))

  # set the rownames to unique gene names
  rownames(data) <- make.unique(feature.names[, gene.column])

  # read in cell barcodes
  cell.names <- readLines(file.path(data_path, 'barcodes.tsv.gz'))

  # if the cell names all have -1 at the end, remove the -1
  if (all(str_detect(string = cell.names, pattern = "\\-1$"))) {
    cell.names <- as.vector(as.character(sapply(
      X = cell.names,
      function(x) str_sub(x, end = (nchar(x) -2))
    )))
  }

  # set the colnames of the gene expression data to the barcodes
  colnames(data) <- cell.names

  # if there is a third column, it indicates that there is more than 1 datatype
  if (ncol(feature.names) > 2) {

    # get the different datatypes
    data_types <- factor(feature.names$V3)
    data_levels <- levels(data_types)

    if (length(data_levels) > 1) {
      message("10X data contains more than one type and is being
              returned as a list containing matrices of each type.")
    }

    # split the data into a list with a matrix for each datatype
    data <- lapply(
      X = data_levels,
      FUN = function(l) {
        return(data[data_types == l, ])
      })

    # name the elements in the list by datatype
    names(data) <- data_levels

  } else {

    data <- list(`Gene Expression` = data)

  }

  return(data)
}

#' Reads in count data from 10x, or from a delimited file.
#'
#' @param sample_names A character vector of sample names.
#' @param data_path_10x A path to the data /<data_path>/outs.
#' @param AC_text_file A path to the antibody capture file.
#' @param delim A delimiter for the AC_test_file.
#' @param log_file A log filename.
#'
#' @return A named list of matrices.
#'
#' @import dplyr
#' @importFrom glue glue
#' @importFrom Seurat Read10X
#' @importFrom stringr str_subset str_c
#' @export
load_sample_counts_matrix = function(sample_names, data_path_10x = NULL, AC_text_file = NULL, delim = NULL, log_file = NULL) {
  # Reads in count data from 10x from one path or multiple paths
  # FROM IGOR DOLGALEV
  # Args:
  #   sample_names: Names to set each sample
  #   data_path: Paths to each sample /data_path/outs
  #   log_file: Name of log_file
  #
  # Returns:
  #   Counts matrix from 10x, merged from several samples or from one

  # Only does 10x or AC
  if(!is.null(data_path_10x) & !is.null(AC_text_file)) {
    stop("One file at a time please.")
  }

  if(!is.null(data_path_10x) & !is.null(delim)) {
    warning("delim will not be taken into account for reading 10x")
  }

  counts = list()

  if(!is.null(data_path_10x)) {

    message_str <- "\n\n ========== import cell ranger counts matrix ========== \n\n"
    # write message will output a message, and write to a log file if a log file is
    # supplied
    write_message(message_str, log_file)

    message_str <- glue("loading counts matrix for sample: {sample_names}")
    write_message(message_str, log_file)

    # check if sample dir is valid
    if (!dir.exists(data_path_10x)) {
      stop(glue("dir {data_path_10x} does not exist"))
    }

    # check if there is an outs directory
    data_dir = glue("{data_path_10x}/outs")

    if (!dir.exists(data_dir)) {
      stop(glue("dir {data_path_10x} does not contain outs directory"))
    }

    # check if the outs directory has the necesary elements
    data_dir = list.files(path = data_dir,
                          pattern = "matrix.mtx",
                          full.names = TRUE,
                          recursive = TRUE)

    data_dir = str_subset(data_dir, "filtered_.*_bc_matrix")[1]
    data_dir = dirname(data_dir)

    if (!dir.exists(data_dir)) {
      stop(glue("dir {data_path_10x} does not contain matrix.mtx"))
    }

    message_str <- glue("loading counts matrix dir: {data_dir}")
    write_message(message_str, log_file)

    # read in the 10x counts
    counts_matrix = import_matrix(data_dir)

  } else if(!is.null(AC_text_file)) {

    # read in AC data
    AC <- read.table(AC_text_file, sep = delim,
                     header = TRUE)

    # remove any row that is called unmapped
    unmpp <- str_detect(toupper(rownames(AC)), "UNMAPPED")
    if(sum(unmpp) > 0) {
      AC <- AC[!unmpp,]
    }

    counts_matrix <- list(`Antibody Capture` = AC)

  }

  # rename the column names by adding a sample name prefix to every barcode
  counts_out <- list()
  for(i in 1:length(counts_matrix)) {
    current_mat <- counts_matrix[[i]]
    colnames(current_mat) <- str_c(sample_names, ":", colnames(current_mat))
    counts_out[[i]] <- current_mat
  }

  return(counts_out)
}


#' Create a new Seurat object from a matrix.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param assay assay.
#' @param log_file log file.
#'
#' @return Seurat object.
#'
#' @import Matrix dplyr readr tibble
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject AddMetaData
#' @export
create_seurat_obj <- function(counts_matrix, assay = "RNA",
                               log_file = NULL) {

  # check that the size of the input matrix is reasonable
  if (ncol(counts_matrix) < 10) {
    stop(glue("matrix contains too few cells: {ncol(counts_matrix)}"))
  }

  # remove genes with very few counts
  counts_matrix <- counts_matrix[Matrix::rowSums(counts_matrix) > 0, ]

  message_str <- glue("\n\n ========== create seurat object ========== \n\n
                     input cells: {ncol(counts_matrix)}
                     input genes: {nrow(counts_matrix)}")
  write_message(message_str, log_file)

  # Create seurat object
  s_obj <- CreateSeuratObject(
    counts = counts_matrix,
    project = "proj",
    assay = assay,
    names.field = 1,
    names.delim = ":"
  )

  # Calculate mito pct
  if(assay == "RNA"){
    s_obj <- calculate_mito_pct(s_obj)
  }

  return(s_obj)
}

#' Calculate mitochondrial percentage from Seurat object.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return Seurat object.
#'
#' @import Matrix
#' @importFrom Seurat AddMetaData
#' @export
calculate_mito_pct <- function(seurat_obj){

  # nGene and nUMI are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  s_obj <- seurat_obj

  # get all mitochondrial genes
  mt_genes <- grep("^MT-", rownames(s_obj@assays$RNA@counts),
                   ignore.case = TRUE, value = TRUE)

  # calculate the percent mitochondrial reads
  percent_mt <- Matrix::colSums(s_obj@assays$RNA@counts[mt_genes, ]) / Matrix::colSums(s_obj@assays$RNA@counts)
  percent_mt <- round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj <- AddMetaData(s_obj, metadata = percent_mt, col.name = "percent.mito")

 return(s_obj)
}

#' Add assay to Seurat object. If Seurat object doesn't exist, create one.
#'
#' @param seurat_obj A Seurat object.
#' @param assay Assay slot.
#' @param counts_matrix counts matrix.
#' @param log_file log file.
#'
#' @return Seurat object.
#'
#' @importFrom Seurat CreateAssayObject
#' @export
add_seurat_assay <- function(seurat_obj, assay, counts_matrix, log_file = NULL){

  #If s_obj is not a seurat object, create one
  if(!is(seurat_obj, "Seurat")){

    s_obj <- create_seurat_obj(counts_matrix = counts_matrix,
                               assay = assay, log_file = log_file)

  } else if(is(seurat_obj, "Seurat")){

    s_obj <- seurat_obj
    cells_to_use <- intersect(colnames(s_obj), colnames(counts_matrix))

    if(length(s_obj) != length(cells_to_use)){
      message_str <- "some cells in scrna matrix not in counts matrix"
      write_message(message_str, log_file)
    }
    if(ncol(counts_matrix) != length(cells_to_use)){
      message_str <- "some cells in counts matrix not in scrna matrix"
      write_message(message_str, log_file)
    }

    # Subset  counts by joint cell barcodes
    counts_matrix <- as.matrix(counts_matrix[, cells_to_use])
    s_obj <- subset(s_obj, cells = cells_to_use)

    # add ADT slot
    s_obj[[assay]] <- CreateAssayObject(counts = counts_matrix)
  }
}


#' Filter out cells based on minimum and maximum number of genes and max mito percentage.
#'
#' @param metadata_tbl A tibble with metadata.
#' @param min_genes Minimum number of genes per cell.
#' @param max_genes Maximim number of genes per cell.
#' @param max_mt Maximum percentage of mitochondrial reads per cell.
#' @param log_file log file.
#'
#' @return cells to keep
#'
#' @import dplyr
#' @importFrom stats quantile
#' @export
filter_data <- function(metadata_tbl, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10) {
  UseMethod("filter_data")
}

filter_data.default <- function(metadata_tbl, log_file = NULL, min_genes, max_genes, max_mt) {

  cells_subset =
    metadata_tbl %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    filter(.data$nFeature_RNA > min_genes,
           .data$nFeature_RNA < max_genes,
           .data$percent.mito < max_mt) %>%
    pull(.data$cell)

  return(cells_subset)
}

#' @export
filter_data.Seurat <- function(metadata_tbl, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10) {
  # filter data by number of genes and mitochondrial percentage
  #
  # Args:
  #   seurat_obj: Seurat object
  #   out_dir: Output directory
  #   proj_name: Name or project and name of output files
  #   min_genes: Minimum number of genes
  #   max_genes: Maximum number of genes
  #   max_mt: Maximum mito pct
  #
  # Results:
  #   Filtered seurat object

  s_obj = metadata_tbl

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

  s_obj = subset(s_obj, cells = cells_subset)

  message_str <- glue("filtered cells: {ncol(s_obj)}
                      filtered genes: {nrow(s_obj)}")
  return(s_obj)
}


normalize_data <- function(data, normalize_method, metadata, assay, slot, log_file = NULL) {
  UseMethod("normalize_data")
}

normalize_data.default <- function(data, normalize_method, metadata, assay, slot, log_file = NULL) {
  if (normalize_method == "log_norm") {
    normed <- log_normalize_data(data = data, log_file = log_file)
  } else if ( normalize_method == "sct") {
    normed <- sctransform_data(umi = data, cell.attr = metadata, variable.features.n = 2000, log_file = log_file)
  }
  return(normed)
}

normalize_data.Seurat <- function(data, normalize_method, metadata, assay, slot, log_file = NULL) {

  if (normalize_method == "sct") {

    assay.obj <- GetAssay(object = data, assay = assay)
    umi <- GetAssayData(object = assay.obj, slot = 'counts')
    cell.attr <- slot(object = data, name = 'meta.data')

    normed_data <- normalize_data(data = umi, normalize_method = normalize_method,
                                  metadata = cell.attr,
                                  log_file = log_file,
                                  slot = slot,
                                  assay = assay)

    assay.out <- CreateAssayObject(counts = normed_data[["vst.out"]]$umi_corrected)


    VariableFeatures(object = assay.out) <- normed_data[["top.features"]]

    assay.out <- SetAssayData(
      object = assay.out,
      slot = 'data',
      new.data = log1p(GetAssayData(object = assay.out, slot = 'counts')))

    assay.out <- SetAssayData(
      object = assay.out,
      slot = 'scale.data',
      new.data = normed_data[["scale.data"]]
    )

  } else {
   data <- GetAssay(data)
   normed_data <- normalize_data(data, normalize_method)
   object[[assay]] <- normed_data
  }
}

#' Log normalize data.
#'
#' @param seurat_obj A seurat object.
#' @param assay assay.
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
#' @param seurat_obj A seurat object.
#' @param log_file log file.
#'
#' @return normalized data
#'
#' @import dplyr
#' @importFrom Seurat PercentageFeatureSet SCTransform
#' @export
sctransform_data <- function(umi, cell.attr, variable.features.n, log_file = NULL){
  # sc transform data

  message("\n\n ========== sc transform ========== \n\n")

  clip.range <- c(-sqrt(ncol(umi)),sqrt(ncol(umi)))

  vst.out <- sctransform::vst(umi, cell_attr = cell.attr,
                 latent_var = "log_umi",
                 show_progress = TRUE,
                 return_cell_attr = TRUE,
                 return_gene_attr = TRUE,
                 return_corrected_umi = TRUE,
                 residual_type = "pearson",
                 res_clip_range = clip.range)

  feature.variance <- setNames(object = vst.out$gene_attr$residual_variance,
                               nm = rownames(vst.out$gene_attr))

  feature.variance <- sort(feature.variance,
                           decreasing = TRUE)

  top.features <- names(feature.variance)[1:min(variable.features.n,
                                                    length(feature.variance))]

  scale.data <- vst.out$y

  scale.data[scale.data < clip.range[1]] <- clip.range[1]
  scale.data[scale.data > clip.range[2]] <- clip.range[2]

  scale.data <- ScaleData(
    scale.data,
    features = NULL,
    vars.to.regress = c("percent.mito", "nCount_RNA"),
    latent.data = cell.attr[, c("percent.mito", "nCount_RNA"), drop = FALSE],
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
