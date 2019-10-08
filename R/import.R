#' Reads in data from 10x, or from a delimited file.
#'
#' @param sample_name A character for the sample name to be a prefix for every cell.
#' @param path A path to the 10x data /<data_path>/outs.
#' @param log_file A log filename.
#'
#' @return Named list of matrices. One matrix for each data type.
#'
#' @import dplyr
#' @importFrom glue glue
#' @importFrom Seurat Read10X
#' @importFrom stringr str_subset str_c
#' @importFrom data.table fread
#' @export
load_sample_counts_matrix = function(sample_name, path, log_file = NULL) {

  counts = list()

  message_str <- glue("loading counts matrix for sample: {sample_name}")
  write_message(message_str, log_file)

  # if path is a file, read it in
  if(file_test("-f", path)) {

    # read in AC data
    counts_df <- fread(path, stringsAsFactors = FALSE, data.table = FALSE)

    # remove any row that is called unmapped
    unmapped <- str_detect(toupper(rownames(counts_df)), "UNMAPPED")
    if(any(unmapped)) {
      counts_df <- counts_df[!unmapped,]
    }

    counts_matrix <- list("Antibody Capture" = counts_df)

  } else {
    # check if path is a directory
    if (!dir.exists(path)) {
      stop(glue("{path} is not a file and is not a directory"))
    }
    # check if the outs directory has the necesary elements
    data_dir = list.files(path = path,
                          pattern = "matrix.mtx",
                          full.names = TRUE,
                          recursive = TRUE)
    data_dir = str_subset(data_dir, "filtered_.*_bc_matrix")[1]
    data_dir = dirname(data_dir)
    if (!dir.exists(data_dir)) {
      stop(glue("dir {data_dir} does not contain matrix.mtx"))
    }
    message_str <- glue("loading counts matrix dir: {data_dir}")
    write_message(message_str, log_file)

    # read in the 10x counts
    counts_matrix = import_mtx(data_dir)
  }

  # rename the column names by adding a sample name prefix to every barcode
  counts_out <- list()
  for(i in 1:length(counts_matrix)) {
    current_mat <- counts_matrix[[i]]
    colnames(current_mat) <- str_c(sample_name, ":", colnames(current_mat))
    counts_out[[i]] <- current_mat
  }
  names(counts_out) <- names(counts_matrix)

  return(counts_out)
}

#' Read in 10x Matrix Market format data.
#'
#' @param data_path Path to directory that holds the files output from 10x.
#' @param gene.column The column with the gene names.
#'
#' @return Named list of matrices. One matrix for each data type as specified in
#' the third column of the features.tsv file. As of Oct 3rd 2019, the two options
#' are `Gene Expression` and `Antibody Capture`
#'
#' @importFrom Matrix readMM
#' @export
import_mtx <- function(data_path, gene_column = 2, log_file = NULL){
  # Heavily sourced from Seurat
  # caveat the name of the files have to be features.tsv.gz, matrix.mtx.gz, barcodes.tsv.gz

  message_str <- "\n\n ========== import cell ranger counts matrix ========== \n\n"
  # write message will output a message, and write to a log file if a log file is
  # supplied
  write_message(message_str, log_file)

  # check if the directory exists
  if (!dir.exists(paths = data_path)) {
    stop(glue("dir {data_path} does not exist"))
  }

  # read in features.tsv.gz
  feature_names <- read.delim(
    file = file.path(data_path, 'features.tsv.gz'),
    header = FALSE,
    stringsAsFactors = FALSE)

  # read in sparse gene expression matrix
  data <- readMM(file = file.path(data_path, 'matrix.mtx.gz'))

  # set the rownames to unique gene names
  rownames(data) <- make.unique(feature_names[, gene_column])

  # read in cell barcodes
  cell_names <- readLines(file.path(data_path, 'barcodes.tsv.gz'))

  # if the cell names all have -1 at the end, remove the -1
  if (all(str_detect(string = cell_names, pattern = "\\-1$"))) {
    cell_names <- as.vector(as.character(sapply(
      X = cell_names,
      function(x) str_sub(x, end = (nchar(x) -2))
    )))
  }

  # set the barcodes as colnames of the gene expression data
  colnames(data) <- cell_names

  # a third column indicates that there is more than 1 datatype
  if (ncol(feature_names) > 2) {

    # get the different datatypes
    data_types <- factor(feature_names$V3)
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

#' Create a new Seurat object from a matrix.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param assay Seurat assay to add the data to.
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
  s_obj <- AddMetaData(s_obj, metadata = percent_mt, col.name = "pct_mito")

 return(s_obj)
}

#' Adds assay to Seurat object.
#'
#' @param seurat_obj Existing seurat object.
#' @param assay Seurat assay to add the matrix to.
#' @param counts_matrix counts matrix.
#' @param log_file log file.
#'
#' @return Seurat object.
#'
#' @importFrom Seurat CreateAssayObject
#' @importFrom methods is
#' @export
add_seurat_assay <- function(seurat_obj, assay, counts_matrix, log_file = NULL){

  if(!is(seurat_obj, "Seurat")){

    stop(glue("{seurat_obj} is not a Seurat object. Cannot add Assay"))

  }

  if(assay %in% names(seurat_obj)) {
    stop(glue("{assay} already exists in this Seurat object. Overwriting {assay}"))
  }

  s_obj <- seurat_obj

  # use cells that are found in both antibody capture and RNA
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

  # add assay
  s_obj[[assay]] <- CreateAssayObject(counts = counts_matrix)

  return(s_obj)
}

