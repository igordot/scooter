#' Plot the distribution of specified features/variables.
#'
#' @param so Seurat object.
#' @param features Vector of features to plot (such as "nGene", "nUMI", "percent.mito").
#' @param grouping X.
#' @param color_scheme (optional) Vector of colors.
#'
#' @return A vector of colors.
#'
#' @examples
#' plot_distribution(so, features, grouping)
#'
#' @import ggplot2
#' @export
plot_distribution <- function(so, features, grouping, filename, color_scheme = NULL, width = NA, height = NA) {

  # color scheme
  if (is.null(color_scheme)) color_scheme <- get_color_scheme()

  # compile the data table
  dist_tbl <- Seurat::FetchData(object = so, vars.all = c(features, grouping))
  dist_tbl <- dist_tbl %>% tidyr::gather(key = "var", value = "val", -grouping)

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
  cowplot::save_plot(filename = filename, plot = plot_dist, base_width = width, base_height = height)

}
