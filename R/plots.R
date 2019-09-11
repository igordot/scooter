#' Determine the color scheme.
#'
#' Determine the color scheme. Can be specified for samples or clusters to avoid confusion.
#'
#' @param type Type of scheme ("samples" or "clusters").
#'
#' @return A vector of colors.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggsci pal_d3 pal_igv
get_color_scheme <- function(type = "clusters") {
  if (type == "samples") {
    color_scheme <- c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  }
  if (type == "clusters") {
    color_scheme <- c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))
  }

  return(color_scheme)
}

#' Determine the point size for tSNE plots (smaller for larger datasets).
#'
#' @param num_cells Number of cells (points on the plot).
#'
#' @return Numeric point size.
get_tsne_point_size <- function(num_cells) {
  pt_size <- 1

  # bigger for smaller datasets
  if (num_cells < 5000) pt_size <- 1.2
  if (num_cells < 3000) pt_size <- 1.5
  if (num_cells < 1000) pt_size <- 2

  # smaller for larger datasets
  if (num_cells > 10000) pt_size <- 0.8
  if (num_cells > 25000) pt_size <- 0.6

  return(pt_size)
}

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
