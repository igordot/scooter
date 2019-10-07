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
save_seurat_counts_matrix <- function(seurat_obj, proj_name = "", label = "", out_dir = ".", assay = "RNA", slot = "data", log_file = NULL) {
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
#' @param metadata1 A Seurat object or a tibble containing metadata with
#' either a column called "cell" with cell IDs or rownames with cell IDs.
#' @param metadata2 A tibble containing metadata with
#' either a column called "cell" with cell IDs or rownames with cell IDs.
#' @param log_file A log filename.
#'
#' @return A metadata file merged on cell identifiers.
#'
#' @import dplyr
#' @importFrom glue glue
#' @importFrom stringr str_detect
#' @export
merge_metadata <- function(metadata1, metadata2, log_file = NULL) {
  UseMethod("merge_metadata")
}

#' @export
merge_metadata.default <- function(metadata1, metadata2, log_file = NULL) {

  message_str <- "\n\n ========== saving metadata ========== \n\n"
  write_message(message_str, log_file)

  # check that there is a cell column, if not, make the rownames the cell column
  if(any(str_detect("^cell$", colnames(metadata1)))) {
    metadata1 = metadata1 %>%
      as_tibble()
  } else {
    metadata1 = metadata1 %>%
      as.data.frame %>%
      rownames_to_column("cell") %>%
      as_tibble()
  }

  # check that there is a cell column, if not, make the rownames the cell column
  if(any(str_detect("^cell$", colnames(metadata2)))){
    metadata2 = metadata2 %>%
      as_tibble()
  } else {
    metadata2 = metadata2 %>%
      as.data.frame %>%
      rownames_to_column("cell") %>%
      as_tibble()
  }

  # compile all cell metadata into a single table
  cells_metadata = metadata1 %>%
    full_join(metadata2, by = "cell") %>%
    arrange(.data$cell) %>%
    as.data.frame()


  return(cells_metadata)
}

#' @export
merge_metadata.Seurat <- function(metadata1, metadata2, log_file = NULL) {

  message_str <- "\n\n ========== saving metadata ========== \n\n"
  write_message(message_str, log_file)

  metadata1 = metadata1@meta.data
  merge_metadata(metadata1, metadata2, log_file = log_file)
}

#' Function to extract data from Seurat object.
#'
#' @param seurat_obj A Seurat object.
#' @param assay Assay such as RNA.
#' @param slot Slot such as counts.
#' @param reduction Character vector of reduction types.
#' @param metadata Boolean. To grab metadata or not
#'
#' @return A metadata file merged on cell identifiers.
#'
#' @import dplyr
#' @importFrom purrr reduce
#' @export
seurat_to_matrix <- function(seurat_obj, assay = NULL, slot = NULL, reduction = NULL, metadata = TRUE) {
  #TODO THIS FUNCTION IS NOT DONE
  s_obj <- seurat_obj

  if(metadata == TRUE) {
    s_obj_metadata <- s_obj@meta.data %>%
      rownames_to_column("cell")
  }

  if(!is.null(reduction)) {
   # get index of reductions in the list
   reductions_to_save <- lapply(reduction, function(x){
     as.data.frame(s_obj@reductions[[x]]@cell.embeddings)
   })

   # set rownames to columns for easier joining
   reductions_to_save <- lapply(reductions_to_save,
                                function(x) rownames_to_column(.data = x,
                                                               "cell"))
   # join lists by the new column
   s_obj_reduction <- reduce(reductions_to_save,
                             .f = left_join,
                             by = "cell")
  }

  if(!is.null(assay)) {
    if(!is.null(slot)) {
      s_obj_assay <- s_obj@assays$assay$slot
    } else {
      s_obj_assay <- s_obj@assays$assay
    }
  }

  if(!is.null(slot)) {
    if(!is.null(assay)) {
      stop("No assay specified for slot")
    }
  }

  out <- full_join(s_obj_metadata, s_obj_reduction, by = "cell")
  return(out)
}

#' Function to create a color vector.
#'
#' @param seurat_obj A Seurat object.
#' @param group Assay such as RNA.
#'
#' @return A vector of colors.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggsci pal_igv
#' @export
create_color_vect <- function(seurat_obj, group = "orig.ident") {
  # create a vector of colors for the Idents of the s_obj
  sample_names <- switch(class(seurat_obj),
                         Seurat = seurat_obj[[group]] %>% unique() %>% arrange(get(group)),
                         data.frame = unique(seurat_obj))

  colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  # create a named color scheme to ensure names and colors are in the proper order
  sample_names[] <- lapply(sample_names, as.character)
  colors_samples_named = colors_samples[1:nrow(sample_names)]
  names(colors_samples_named) = sample_names[,1]
  return(colors_samples_named)
}