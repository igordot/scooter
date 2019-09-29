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


