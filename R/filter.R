#' Filter out cells based on minimum and maximum number of genes and maximum percentage mitochondrial reads. If cutoffs are not provided, the min_genes will be the 0.02 quantile, and the max genes will be 0.98 quantile and the mitochondrial percentage will be 10%.
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

filter_data.default <- function(data, log_file = NULL, min_genes, max_genes, max_mt) {

  cells_subset =
    data %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    filter(.data$nFeature_RNA > min_genes,
           .data$nFeature_RNA < max_genes,
           .data$pct_mito < max_mt) %>%
    pull(.data$cell)

  return(cells_subset)
}

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

