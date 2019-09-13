#' Small function to write to message and to log file.
#'
#' @param message_str A string to write as a message.
#' @param log_file A log filename.
#
#' @return A message and writes the message to the specified log file.
#'
#' @examples
#' write_message(message_str = "Finished Step 1", log_file = "log.file.txt")
#'
#' @export
write_message <- function(message_str, log_file = NULL) {
  # Small function to write to message and to log file if log file is not null
  message(message_str)
  if(!is.null(log_file)){
    write(message_str,
          file = log_file,
          append = TRUE)
  }
}
#' Function to write Seurat counts matrix to csv.
#'
#' @param seurat_obj A Seurat object.
#' @param proj_name Name of the project that will be the prefix of the file name.
#' @param label An optional label for the file.
#' @param out_dir Directory in which to save csv.
#' @param assay The assay within the Seurat object to retrieve data from.
#' @param slot The slot within the Seurat object to retrieve data from.
#' @param log_file A log filename.
#'
#' @return A csv file in the out_dir.
#'
#' @import dplyr
#' @importFrom glue glue
#' @importFrom Seurat GetAssayData
#' @export
save_counts_matrix <- function(seurat_obj, proj_name = "", label = "", out_dir = ".", assay = "RNA", slot = "data", log_file = NULL) {
  # save counts matrix as a csv file (to be consistent with the rest of the tables)

  message_str <- "\n\n ========== saving Seurat counts ========== \n\n"
  write_message(message_str, log_file)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts = GetAssayData(seurat_obj, assay = assay) %>%
    as.matrix() %>%
    round(3)

  counts = counts %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(gene)

  write_csv(counts,
            path = glue("{out_dir}/{proj_name}.{label}.counts.csv.gz"))

  rm(counts)
}
#' Function to merge two metadata tables together.
#'
#' @param metadata_one A Seurat object or a tibble containing metadata with
#' either a column called "cell" with cell IDs or rownames with cell IDs.
#' @param metadata_two A tibble containing metadata with
#' either a column called "cell" with cell IDs or rownames with cell IDs.
#' @param log_file A log filename.
#' @param proj_name Name of the project that will be the prefix of the file name.
#' @param label An optional label for the file.
#' @param write Boolean. To write out the metadata or not.
#' @param out_dir Directory in which to write metadata.
#'
#' @return A metadata file merged on cell identifiers.
#'
#' @import dplyr
#' @importFrom glue glue
#' @export
merge_metadata <- function(metadata_one, metadata_two, log_file = NULL, write = TRUE, proj_name = "", label = "", out_dir = ".") {
  # save metadata from seurat object

  message_str <- "\n\n ========== saving metadata ========== \n\n"
  write_message(message_str, log_file)

  metadata_one <- switch(class(metadata_one),
                 Seurat = metadata_one@meta.data,
                 data.frame = metadata_one)

  # check that there is a cell column, if not, make the rownames the cell column
  if(sum(grepl("^cell$", colnames(metadata_one)))) {
    metadata_one = metadata_one %>% as.tibble()
  } else {
    metadata_one = metadata_one %>% as.data.frame %>% rownames_to_column("cell") %>%  as.tibble()
  }

  # check that there is a cell column, if not, make the rownames the cell column
  if(sum(grepl("^cell$", colnames(metadata_two)))){
    metadata_two = metadata_two %>%  as.tibble()
  } else {
    metadata_two = metadata_two %>% as.data.frame %>% rownames_to_column("cell") %>%  as.tibble()
  }

  # compile all cell metadata into a single table
  cells_metadata = metadata_one %>%
    full_join(metadata_two,by = "cell") %>%
    arrange(cell) %>%
    as.data.frame()

  if(write) {
    write_excel_csv(cells_metadata, path = glue("{out_dir}/{proj_name}.{label}.metadata.csv"))
  }

  return(cells_metadata)
}