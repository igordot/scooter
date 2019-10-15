#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param pcs Data.
#' @param num_dim Number of PCs to use for tsne and umap.
#' @param num_neighbors Number of neighbors to use for umap.
#' @param log_file log file.
#' @param res .
#'
#' @return .
#'
#' @import dplyr
#' @import Seurat
#' @export
calculate_clusters <- function(pcs, num_dim, log_file, num_neighbors = 30, res = NULL){
  # TODO: allow UMAP graphs to be used

  message_str <- "========== Seurat::FindNeighbors() =========="
  write_message(message_str, log_file)

  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)

  snn_graph <- Seurat:::FindNeighbors.default(
    pcs[,1:num_dim],
    distance.matrix = FALSE,
    k.param = num_neighbors,
    compute.SNN = TRUE,
    prune.SNN = 1/15,
    nn.eps = 0,
    force.recalc = TRUE)

  snn_graph <- as.Graph(snn_graph[["snn"]])

  message_str <- "\n\n ========== Seurat::FindClusters() ========== \n\n"
  write_message(message_str, log_file)

  if(is.null(res)){
    res_range <- seq(0.1, 2.5, 0.1)
    if (nrow(pcs) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8, 9)
  }else{
    res_range <- res
  }

  clusters <-  Seurat:::FindClusters.default(snn_graph,
                                             algorithm = 3,
                                             resolution = res_range,
                                             verbose = FALSE)
  return(clusters)
}

#' Calculate differential expression for
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param metadata_column Column in metadata.
#' @param list_groups list of groups to compare in the metadata column.
#' @param proj_name Project name as a prefix for the output file.
#' @param log_fc_threshold Log fc threshold.
#' @param min.pct Minimum percentage of cells a gene must be expressed in to be tested.
#' @param test.use Test to use.
#' @param out_path output path.
#' @param write Boolean to write or not.
#' @param log_file log file.
#'
#' @return .
#'
#' @import dplyr
#' @import Seurat
#' @export
differential_expression <- function(data, metadata, metadata_column, list_groups = NULL, log_fc_threshold = 0.5, min.pct = 0.1, test.use = "wilcox", out_path = ".", write = FALSE, log_file = NULL){


  if(!is.null(list_groups)) {
    # get unique combinations
    unique_ids <- unique(unique(metadata[,metadata_column]))[[1]]
    list_groups <- expand.grid(unique_ids,
                               unique_ids,
                               stringsAsFactors = FALSE)
    list_groups <- list_groups %>%
      filter(Var1 != Var2)
    indx <- !duplicated(t(apply(list_groups, 1, sort)))
    list_groups <- list_groups[indx,]
  } else {
    list_groups <- list_groups
  }

  # If there is a column called cell leave it, if not, make rownames cell
  if("cell" %in% colnames (metadata)) {
      metadata = metadata %>%  as.tibble()
  } else {
    metadata = metadata %>% as_tibble(rownames = "cell")
  }

  # initialize a tibble to save de in
  diff_exp_stats <- tibble(
    gene = character(),
    p_val = numeric(),
    avg_logFC = numeric(),
    p_val_adj = numeric(),
    group_1 = character(),
    group_2 = character()
  )

  for(current_group in list_groups){

    cell_group1 <-  metadata %>%
      filter(get(metadata_column) == current_group[1]) %>%
      select("cell")

    cell_group2 <-  metadata %>%
      filter(get(metadata_column) == current_group[2]) %>%
      select("cell")

    message_str <- glue("{current_group[1]} versus {current_group[2]}")
    write_message(message_str, log_file)

    current_comparison <- Seurat:::FindMarkers.default(
      object = as.matrix(data),
      reduction = NULL,
      slot = "data",
      cells.1 = cell_group1$cell,
      cells.2 = cell_group2$cell,
      logfc.threshold = log_fc_threshold,
      test.use = test.use,
      min.pct =  min.pct)

    current_comparison_filt <- current_comparison %>%
      select(p_val, avg_logFC, p_val_adj) %>%
      rownames_to_column("gene") %>%
      mutate(group_1 = rep(current_group[1])) %>%
      mutate(group_2 = rep(current_group[2]))

    if(write == TRUE) {
      write_excel_csv(current_comparison_filt,
                      path = glue("{out_path}/{proj_name}.diff_exp.{metadata_column}{current_group[1]}.{current_group[2]}.csv"))
    }

    diff_exp_stats <- rbind(diff_exp_stats, current_comparison_filt)
  }

  if(write == TRUE) {
    write_excel_csv(diff_exp_stats, path = glue("{out_path}/{proj_name}.diff_exp.{metadata_column}.csv"))
  }

  return(diff_exp_stats)
}

#' Get cluster averages.
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param group Column in metadata.
#'
#' @return .
#'
#' @import dplyr
#' @import data.table
#' @export
calc_clust_averages <- function(metadata, data, group){

  # get relevant metadata
  metadata <- metadata %>%
    select("cell", group)

  # merge metadata and data on cell
  if(nrow(metadata) != (ncol(data) - 1)) {
    stop("the number of cells in metadata is not the same as the number of cells in data")
  }

  # manitpulate data to merge with metadata
  data <- data %>%
    column_to_rownames("gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell")

  current_data <- full_join(metadata, data, by = "cell") %>%
    column_to_rownames("cell")

  current_data <- setDT(current_data)
  current_mean <- current_data[, lapply(.SD, mean), by = .(get(group)), .SDcols = 2:ncol(current_data)]

  colnames(current_mean)[which(colnames(current_mean) == "get")] <- group

  current_mean <- as.data.frame(current_mean)

  return(current_mean)
}

