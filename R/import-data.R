#' Create a new Seurat object from a matrix.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param min_cells X.
#' @param min_genes X.
#' @param out_dir Directory where all the output files will be saved.
#' @param color_scheme A character vector of colors for plots (per sample/library).
#'
#' @return Seurat object.
#'
#' @import Matrix dplyr readr tibble
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject AddMetaData
#' @export
create_seurat_obj <- function(counts_matrix, min_cells = 10, min_genes = 100, out_dir = ".", color_scheme = NULL, log_file = NULL) {

  # create output directories
  if (!dir.exists(out_dir)) dir.create(out_dir)
  qc_dir = glue("{out_dir}/qc")
  if (!dir.exists(qc_dir)) dir.create(qc_dir)

  # color scheme
  if (is.null(color_scheme)) color_scheme <- get_color_scheme(type = "samples")

  # check that the size of the input matrix is reasonable
  if (ncol(counts_matrix) < 10) stop(glue("matrix contains too few cells: {ncol(counts_matrix)}"))
  if (nrow(counts_matrix) < 100) stop(glue("matrix contains too few genes: {nrow(counts_matrix)}"))

  # remove genes with very few counts
  counts_matrix <- counts_matrix[Matrix::rowSums(counts_matrix) > 0, ]

  message_str <- glue("\n\n ========== create seurat object ========== \n\n
                     input cells: {ncol(counts_matrix)}
                     input genes: {nrow(counts_matrix)}")
  write_message(message_str, log_file)

  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  raw_data <- counts_matrix %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(raw_data, path = "counts.raw.csv.gz")

  s_obj <- CreateSeuratObject(
    counts = counts_matrix,
    min.cells = min_cells,
    min.features = min_genes,
    project = "proj",
    names.field = 1,
    names.delim = ":"
  )

  # nGene and nUMI are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  mt_genes <- grep("^MT-", rownames(s_obj@assays$RNA@counts), ignore.case = TRUE, value = TRUE)
  percent_mt <- Matrix::colSums(s_obj@assays$RNA@counts[mt_genes, ]) / Matrix::colSums(s_obj@assays$RNA@counts)
  percent_mt <- round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj <- AddMetaData(s_obj, metadata = percent_mt, col.name = "percent.mito")

  plot_distribution(
    so = s_obj, features = c("nGene", "nUMI", "percent.mito"), grouping = "orig.ident",
    color_scheme = color_scheme,
    filename = glue("{qc_dir}/qc.distribution.unfiltered.png"), width = 10, height = 5
  )

  return(s_obj)
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
#' @export
load_sample_counts_matrix = function(sample_names, data_path, log_file = NULL) {
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

    # determine counts matrix directory (HDF5 is not the preferred option)
    # "filtered_gene_bc_matrices" for single library
    # "filtered_gene_bc_matrices_mex" for aggregated
    # Cell Ranger 3.0: "genes" has been replaced by "features" to account for feature barcoding
    # Cell Ranger 3.0: the matrix and barcode files are now gzipped
    data_dir = glue("{data_path}/outs")
    if (!dir.exists(data_dir)) stop(glue("dir {data_path} does not contain outs directory"))
    data_dir = list.files(path = data_dir, pattern = "matrix.mtx", full.names = TRUE, recursive = TRUE)
    data_dir = str_subset(data_dir, "filtered_.*_bc_matri")[1]
    data_dir = dirname(data_dir)
    if (!dir.exists(data_dir)) stop(glue("dir {data_path} does not contain matrix.mtx"))

    message_str <- glue("loading counts matrix dir: {data_dir}")
    write_message(message_str, log_file)


    counts_matrix = Read10X(data_dir)

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

        message_str <- glue("num genes for previous libraries: {length(rownames(counts_matrix_total))}
                            num genes for current library: {length(rownames(counts_matrix))}
                            num genes in common: {length(common_genes)}")
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

#' Filter out cells based on minimum and maximum number of genes and max mito percentage.
#'
#' @param metadata_tbl A tibble with metadata.
#' @param min_genes Minimum number of genes per cell.
#' @param max_genes Maximim number of genes per cell.
#' @param max_mt Maximum percentage of mitochondrial reads per cell.
#'
#' @return cells to keep
#'
#' @import dplyr
#' @export
filter_data <- function(metadata, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10) {
  UseMethod("filter_data")
}

filter_data.default <- function(metadata_tbl, log_file = NULL, min_genes, max_genes, max_mt) {

  cells_subset =
    metadata_tbl %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    filter(nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mito < max_mt) %>%
    pull(cell)

  return(cells_subset)
}

#' @return A filtered Seurat object
#' @export
filter_data.Seurat <- function(metadata, log_file = NULL, min_genes = NULL, max_genes = NULL, max_mt = 10) {
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

  s_obj = metadata

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
