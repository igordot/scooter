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

loop_plot_scatter <- function(metadata, out_path, proj_name, log_file, X, Y, colors_vect){

  for(i in colors_vect){
    current_plot <- plot_scatter(metadata,
                                 out_path,
                                 proj_name,
                                 log_file,
                                 X,
                                 Y,
                                 i,
                                 write = TRUE)
  }

}

diff_exp <- function(data, metadata, metadata_column, list_groups, test.use = "wilcox", out_path, proj_name){

  if(sum(grepl("^cell$", colnames(metadata)))){
    metadata = metadata %>%  as.tibble()

  } else {
    metadata = metadata %>% as.data.frame %>% rownames_to_column("cell") %>%  as.tibble()
  }

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
    print(dim(cell_group2))
    print(dim(cell_group1))

    current_comparison <- Seurat:::FindMarkers.default(
      object = as.matrix(data),
      reduction = NULL,
      slot = "data",
      cells.1 = cell_group1$cell,
      cells.2 = cell_group2$cell,
      logfc.threshold = 0.00001,
      test.use = test.use,
      min.pct =  0.0001)

    current_comparison_filt <- current_comparison %>%
      select(p_val, avg_logFC, p_val_adj) %>%
      rownames_to_column("gene") %>%
      mutate(group_1 = rep(current_group[1])) %>%
      mutate(group_2 = rep(current_group[2]))

    write_excel_csv(current_comparison_filt, path = glue("{out_path}/{proj_name}.diff_exp.{metadata_column}{current_group[1]}.{current_group[2]}.csv"))

    diff_exp_stats <- rbind(diff_exp_stats, current_comparison_filt)
  }
  write_excel_csv(diff_exp_stats, path = glue("{out_path}/{proj_name}.diff_exp.{metadata_column}.csv"))

  return(diff_exp_stats)
}

cluster_stats_bar <- function(metadata, group1, group2, write = FALSE, g1_col = NULL, g2_col = NULL, cluster = TRUE){
  # TODO: pull out plots into new function
  # make barplots and output cluster stats
  summary_metadata <- metadata %>%
    group_by(!!!syms(c(group1, group2))) %>%
    summarize(n = n()) %>%
    group_by(!!sym(group1)) %>%
    mutate(pct_g2_in_g1 = n / sum(n)) %>%
    group_by(!!sym(group2)) %>%
    mutate(pct_g1_in_g2 = n / sum(n)) %>%
    ungroup()

  write_excel_csv(summary_metadata, path = glue("{out_path}/{proj_name}.summary.{group1}{group2}.csv"))

  if(write == TRUE){
    # make both grouping variables factors
    summary_metadata %<>% mutate(!!group1 := as.factor(!!sym(group1)))
    summary_metadata %<>% mutate(!!group2 := as.factor(!!sym(group2)))
    if(cluster){
      mat_g1 = summary_metadata %>%
        select(!!c(group1, group2, "pct_g1_in_g2")) %>%
        spread(group2, 'pct_g1_in_g2', fill = 0) %>%
        as.data.frame %>%
        column_to_rownames(group1) %>%
        as.matrix()

      hc_g1 = hclust(dist(mat_g1), method = 'ward.D2')  # clusters rows of mat
      levels_g1 = rownames(mat_g1)[order.dendrogram(as.dendrogram(hc_g1))]

      summary_metadata <- summary_metadata %>%
        mutate(!!group1 := fct_relevel(!!sym(group1), levels_g1))

      mat_g2 = summary_metadata %>%
        select(!!c(group1, group2, "pct_g2_in_g1")) %>%
        spread(group1, 'pct_g2_in_g1', fill = 0) %>%
        as.data.frame %>%
        column_to_rownames(group2) %>%
        as.matrix()

      hc_g2 = hclust(dist(mat_g2), method = 'ward.D2')  # clusters rows of mat
      levels_g2 = rownames(mat_g2)[order.dendrogram(as.dendrogram(hc_g2))]

      summary_metadata <- summary_metadata %>%
        mutate(!!group2 := fct_relevel(!!sym(group2), levels_g2))
    }
    # use levels to re-order factor
    if(is.null(g1_col)){
      group1_col <- create_color_vect(as.data.frame(summary_metadata[group1]))
    } else{
      group1_col <- g1_col
    }
    if(is.null(g2_col)){
      group2_col <- create_color_vect(as.data.frame(summary_metadata[group2]))
    } else{
      group2_col <- g2_col
    }

    summary_plots_g2 <- ggplot(summary_metadata) +
      geom_col(aes_string(x = group2, y = "pct_g1_in_g2", fill = group1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = group1_col,
                        name = group1) +
      ylab(glue("percent {group1} in {group2}"))

    summary_plots_g2_legend <- get_legend(summary_plots_g2)


    summary_plots_g1 <- ggplot(summary_metadata) +
      geom_col(aes_string(x = group1, y = "pct_g2_in_g1", fill = group2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = group2_col,
                        name = group2) +
      ylab(glue("percent {group2} in {group1}"))

    summary_plots_g1_legend <- get_legend(summary_plots_g1)


    summary_plots <- plot_grid(summary_plots_g2 + theme(legend.position = "none"),
                               summary_plots_g2_legend,
                               summary_plots_g1 + theme(legend.position = "none"),
                               summary_plots_g1_legend)

    ggsave(summary_plots,
           file = glue("{out_path}/{proj_name}.{group1}{group2}.bar.png"),
           height = 10,
           width = 10)
  }
  return(summary_metadata)
}

calc_clust_averages <- function(metadata, data, group){
  # merge relevant metadata and data and row avg in group

  # get relevant metadata
  metadata <- metadata %>%
    select("cell", group)

  # merge metadata and data on cell
  if(nrow(metadata) != (ncol(data) - 1)) {stop("the number of cells in metadata is not the same as the number of cells in data")}
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

  write_excel_csv(current_mean, path = glue("{out_path}/{proj_name}.{group}.means.csv"))

  return(current_mean)
}

