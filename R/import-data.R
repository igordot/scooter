#' Create a new Seurat object from a matrix.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param assay assay.
#' @param color_scheme A character vector of colors for plots (per sample/library).
#' @param log_file log file.
#'
#' @return Seurat object.
#'
#' @import Matrix dplyr readr tibble
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject AddMetaData
#' @export
create_seurat_obj <- function(counts_matrix, assay = "RNA",
                              color_scheme = NULL, log_file = NULL) {
  # color scheme
  if (is.null(color_scheme)) color_scheme <- get_color_scheme(type = "samples")

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

  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  raw_data <- counts_matrix %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(.data$gene)

  s_obj <- CreateSeuratObject(
    counts = counts_matrix,
    project = "proj",
    assay = assay,
    names.field = 1,
    names.delim = ":"
  )

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

  mt_genes <- grep("^MT-", rownames(s_obj@assays$RNA@counts),
                   ignore.case = TRUE, value = TRUE)

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

#' Read 10x data.
#'
#' @param data_path path to 10x.
#'
#' @return list of matrices.
#'
#' @importFrom Matrix readMM
#' @export
import_matrix <- function(data_path, gene.column = 2){

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

  rownames(data) <- make.unique(names = feature.names[, gene.column])

  # read in cell names
  cell.names <- readLines(file.path(data_path, 'barcodes.tsv.gz'))

  # if the cell names all have -1 at the end, remove the -1
  if (all(str_detect(string = cell.names, pattern = "\\-1$"))) {
    cell.names <- as.vector(as.character(sapply(
      X = cell.names,
      function(x) str_sub(x, end = (nchar(x) -2))
    )))
  }
  colnames(data) <- cell.names

  if (ncol(feature.names) > 2) {

    data_types <- factor(feature.names$V3)
    data_levels <- levels(data_types)

    if (length(data_levels) > 1 && length(full.data) == 0) {
      message("10X data contains more than one type and is being
              returned as a list containing matrices of each type.")
    }

    data <- lapply(
      X = data_levels,
      FUN = function(l) {
        return(data[data_types == l, ])
      })

    names(data) <- data_levels

  } else{
    data <- list(data)
  }
  return(data)
}


#' Reads in count data from 10x from one path or multiple paths.
#'
#' @param sample_names A character vector of sample names.
#' @param data_path A character vector of paths to each sample /<data_path>/outs.
#' @param log_file A log filename.
#'
#' @return Counts matrix from 10x, merged from several samples or from one.
#'
#' @import dplyr
#' @importFrom glue glue
#' @importFrom Seurat Read10X
#' @importFrom stringr str_subset str_c
#' @export
load_sample_counts_matrix = function(sample_names, data_path_10x, HTO_file, ADT_file, log_file = NULL) {
  # Reads in count data from 10x from one path or multiple paths
  # FROM IGOR DOLGALEV
  # Args:
  #   sample_names: Names to set each sample
  #   data_path: Paths to each sample /data_path/outs
  #   log_file: Name of log_file
  #
  # Returns:
  #   Counts matrix from 10x, merged from several samples or from one

  message_str <- "\n\n ========== import cell ranger counts matrix ========== \n\n"
  # write message will output a message, and write to a log file if a log file is
  # supplied
  write_message(message_str, log_file)

  counts_matrix_total = NULL
  for (i in 1:length(data_path)) {

    sample_name = sample_names[i]

    message_str <- glue("loading counts matrix for sample: {sample_name}")
    write_message(message_str, log_file)

    # check if sample dir is valid
    if (!dir.exists(data_path)) stop(glue("dir {data_path} does not exist"))

    data_dir = glue("{data_path}/outs")

    if (!dir.exists(data_dir)) {
      stop(glue("dir {data_path} does not contain outs directory"))
    }

    data_dir = list.files(path = data_dir,
                          pattern = "matrix.mtx",
                          full.names = TRUE,
                          recursive = TRUE)

    data_dir = str_subset(data_dir, "filtered_.*_bc_matri")[1]
    data_dir = dirname(data_dir)

    if (!dir.exists(data_dir)) stop(glue("dir {data_path} does not contain matrix.mtx"))

    message_str <- glue("loading counts matrix dir: {data_dir}")
    write_message(message_str, log_file)

    counts_matrix = import_matrix(data_dir)

    message_str <- glue("library {sample_name} cells: {ncol(counts_matrix)}
                        library {sample_name} genes: {nrow(counts_matrix)}")
    write_message(message_str, log_file)

    # clean up counts matrix to make it more readable
    counts_matrix = counts_matrix[sort(rownames(counts_matrix)), ]
    colnames(counts_matrix) = str_c(sample_name, ":", colnames(counts_matrix))

    # combine current matrix with previous
    if (i == 1) {

      # skip if there is no previous matrix
      counts_matrix_total = counts_matrix

    } else {

      # check if genes are the same for current and previous matrices
      if (!identical(rownames(counts_matrix_total), rownames(counts_matrix))) {

        # generate a warning, since this is probably a mistake
        warning("counts matrix genes are not the same for different libraries")
        write("counts matrix genes are not the same for different libraries",
              file = log_file,
              append = TRUE)
        Sys.sleep(1)

        # get common genes
        common_genes = intersect(rownames(counts_matrix_total), rownames(counts_matrix))
        common_genes = sort(common_genes)

        message_str <- glue("num genes for previous libraries: {length(rownames(counts_matrix_total))} num genes for current library: {length(rownames(counts_matrix))} num genes in common: {length(common_genes)}")
        write_message(message_str, log_file)

        # exit if the number of overlapping genes is too few
        if (length(common_genes) < (length(rownames(counts_matrix)) * 0.9)) stop("libraries have too few genes in common")

        # subset current and previous matrix to overlapping genes
        counts_matrix_total = counts_matrix_total[common_genes, ]
        counts_matrix = counts_matrix[common_genes, ]
      }

      # combine current matrix with previous
      counts_matrix_total = cbind(counts_matrix_total, counts_matrix)
      Sys.sleep(1)

    }

  }
  return(counts_matrix_total)
}

#' Reads in HTO data.
#'
#' @param HTO_file Path to HTO file.
#'
#' @return hto matrix.
#'
#' @importFrom utils read.table
#' @export
read_remove.unmapped <- function(read_file) {
  # remove unmapped row
  hto_matrix <- read.table(read_file)
  hto_matrix <- hto_matrix[-which(rownames(hto_matrix) == "unmapped"),]
  return(hto_matrix)
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

normalize_data <- function() {

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
log_normalize_data <- function(seurat_obj, assay = "RNA", log_file = NULL) {
  # log normalize data

  message_str <- "\n\n ========== log normalize ========== \n\n"
  write_message(message_str, log_file)

  s_obj <- seurat_obj
  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  s_obj = NormalizeData(s_obj,
                        normalization.method = "LogNormalize",
                        assay = "RNA",
                        scale.factor = 10000,
                        verbose = FALSE)
  # log to file
  message_str <- glue("filtered cells: {ncol(s_obj)}
                      filtered genes: {nrow(s_obj)}
                      filtered mean num genes: {round(mean(s_obj$num_genes), 3)}
                      filtered median num genes: {median(s_obj$num_genes)}")
  write_message(message_str, log_file)

  return(s_obj)
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
sctransform_data <- function(seurat_obj, log_file = NULL){
  # sc transform data

  s_obj <-seurat_obj

  message("\n\n ========== sc transform ========== \n\n")

  s_obj <- PercentageFeatureSet(s_obj,
                                pattern = "^MT-",
                                col.name = "percent.feature.set.mt")
  # run sctransform
  s_obj <- SCTransform(s_obj,
                       vars.to.regress = c("percent.feature.set.mt", "num_UMIs"),
                       verbose = FALSE)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts_norm = GetAssayData(s_obj) %>%
    as.matrix() %>%
    round(3)
  counts_norm = counts_norm %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(gene)

  # log to file
  message_str <- glue("filtered cells: {ncol(s_obj)}
                      filtered genes: {nrow(s_obj)}
                      filtered mean num genes: {round(mean(s_obj$num_genes), 3)}
                      filtered median num genes: {median(s_obj$num_genes)}")
  write_message(message_str, log_file)

  return(s_obj)

}
