#' Plot the distribution of specified features/variables.
#'
#' @param data Seurat object or metadata.
#' @param features Vector of features to plot (such as "nGene", "nUMI", "percent.mito").
#' @param grouping X.
#' @param color_scheme (optional) Vector of colors.
#'
#' @return A vector of colors.
#'
#' @import ggplot2
#' @export
plot_distribution <- function(data, features, grouping, color_scheme = NULL) {
  UseMethod("plot_distribution")
}

#' @export
plot_distribution.default <- function(data, features, grouping, color_scheme = NULL) {
  if (is.null(color_scheme)) color_scheme <- get_color_scheme()

  dist_tbl <- data %>%
    select(c(features, grouping)) %>%
    tidyr::gather(key = "var", value = "val", -grouping)

  plot_dist <-
    ggplot(dist_tbl, aes(x = !!sym(grouping), y = val)) +
    cowplot::theme_cowplot() +
    geom_violin(aes(fill = !!sym(grouping))) +
    geom_jitter(size = 1, color = "black", alpha = 0.3, width = 0.3) +
    scale_fill_manual(values = color_scheme) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank()
    ) +
    facet_wrap(. ~ var, scales = "free")
  return(plot_dist)
}

#' @export
plot_distribution.Seurat <- function(data, features, grouping, color_scheme = NULL) {
  # compile the data table
  dist_tbl <- Seurat::FetchData(object = data, vars = c(features, grouping))
  plot_distribution(dist_tbl, features, grouping, color_scheme)
}

plot_scatter<- function(metadata, out_path = NULL, proj_name = NULL, log_file = NULL, X, Y, color, write = FALSE, color_vect = NULL){

  if(is.null(color_vect)){
    colors_samples_named <- create_color_vect(as.data.frame(metadata[color]))
  } else {
    colors_samples_named <- color_vect
  }

  current_plot <- ggplot(sample_frac(metadata),
                         aes(x = eval(as.name(X)),
                             y = eval(as.name(Y)),
                             color = eval(as.name(color)))) +
    geom_point(size = 1, alpha = 0.5) +
    coord_fixed(ratio = (max(metadata[X]) - min(metadata[X]))/(max(metadata[Y]) - min(metadata[Y]))) +
    xlab(X) +
    ylab(Y) +
    scale_color_manual(values = colors_samples_named,
                       name = color)

  if(write){
    ggsave(glue("{out_path}/{proj_name}.{X}.{Y}.{color}.png"),
           plot = current_plot,
           width = 8,
           height = 6,
           units = "in")
  }

  return(current_plot)
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

  if(write == TRUE) {

    write_excel_csv(summary_metadata, path = glue("{out_path}/{proj_name}.summary.{group1}{group2}.csv"))

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

    if(write == TRUE) {
      ggsave(summary_plots,
             file = glue("{out_path}/{proj_name}.{group1}{group2}.bar.png"),
             height = 10,
             width = 10)
    }
  }
  return(list(summary_metadata = summary_metadata,
              summary_plots = summary_plots))
}

