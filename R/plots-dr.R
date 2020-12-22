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
#' @export
get_color_scheme <- function(type = "clusters") {
  if (type == "samples") {
    color_scheme <- c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  }
  if (type == "clusters") {
    color_scheme <- c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51), pal_igv(alpha = 0.6)(51))
  }

  return(color_scheme)
}

#' Determine the point size for reduced dimensions scatter plots (smaller for larger datasets).
#'
#' @param num_cells Number of cells (points on the plot).
#'
#' @return Numeric point size.
#' @export
get_dr_point_size <- function(num_cells) {

  pt_size <- 1.8
  if (num_cells > 1000) pt_size <- 1.2
  if (num_cells > 5000) pt_size <- 1.0
  if (num_cells > 10000) pt_size <- 0.8
  if (num_cells > 25000) pt_size <- 0.6
  if (num_cells > 50000) pt_size <- 0.4
  return(pt_size)

}

# plot colored by specified variable
plot_scatter_group = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", aspect_ratio = 1, color_var, color_scheme) {

  ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(aes(color = !!sym(color_var)), size = get_dr_pt_size(metadata_tbl)) +
    theme(
      aspect.ratio = aspect_ratio,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_scheme)

}

# plot split by specified variable
plot_scatter_split = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", aspect_ratio = 1, rows_var = NULL, cols_var = NULL, color_var, color_scheme) {

  gp =
    ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(aes(color = !!sym(color_var)), size = get_dr_pt_size(metadata_tbl)) +
    theme(
      aspect.ratio = aspect_ratio,
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      strip.background = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_scheme)

  if (is.null(rows_var)) {
    gp + facet_grid(cols = vars(!!sym(cols_var)))
  } else if (is.null(cols_var)) {
    gp + facet_grid(rows = vars(!!sym(rows_var)))
  } else {
    gp + facet_grid(rows = vars(!!sym(rows_var)), cols = vars(!!sym(cols_var)))
  }

}

# density plot split by specified variable
# calculate density normalized to 1, independently for each facet variable
plot_density_split = function(metadata_tbl, x_var, y_var, split_var, num_bins) {

  # ran into some issues with merging split geom_hex
  ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    # geom_hex(aes(fill = stat(ndensity)), bins = num_bins) +
    stat_bin_2d(aes(fill = stat(ndensity)), bins = num_bins) +
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      strip.background = element_blank()
    ) +
    scale_fill_gradient2(low = "white", high = "darkred") +
    facet_wrap(vars(!!sym(split_var)))

}

# get table for density plot, split by stage
get_density_diff_table = function(metadata_tbl, x_var, y_var, split_var, num_bins) {

  # generate a density plot split by stage
  density_plot = plot_density_split(metadata_tbl = metadata_tbl, x_var = x_var, y_var = y_var, split_var = split_var, num_bins = num_bins)

  # produce an object that can be rendered
  density_plot_tbl = ggplot_build(density_plot)

  # panel labels
  panels_tbl =
    tibble(
      PANEL = density_plot_tbl$layout$layout$PANEL,
      stage = density_plot_tbl$layout$layout$stage
    )

  # merge panel contents and panel names
  density_tbl = density_plot_tbl$data[[1]]
  density_tbl = density_tbl %>% full_join(panels_tbl, by = "PANEL")

  return(density_tbl)

}

# density plot split by specified variable
# split normalization (adding norm_split_var) may not work
plot_density_diff = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", split_var, num_bins, group_pos, group_neg, interpolate = FALSE) {

  density_tbl = get_density_diff_table(metadata_tbl = metadata_tbl, x_var = x_var, y_var = y_var, split_var = split_var, num_bins = num_bins)

  min_density = quantile(density_tbl$density, 0)

  density_pos_tbl =
    density_tbl %>%
    filter(stage == group_pos) %>%
    select(x, y, cells_pos = count, density_pos = density)
  density_neg_tbl =
    density_tbl %>%
    filter(stage == group_neg) %>%
    select(x, y, cells_neg = count, density_neg = density)
  density_split_tbl = full_join(density_pos_tbl, density_neg_tbl, by = c("x", "y"))
  density_split_tbl[is.na(density_split_tbl)] = min_density
  density_split_tbl = density_split_tbl %>% mutate(density_diff = density_pos - density_neg)
  density_split_tbl = density_split_tbl %>% mutate(density_ratio = log(density_pos/density_neg))

  min_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.01)
  max_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.99)
  min_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.01)
  max_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.99)

  density_split_tbl =
    density_split_tbl %>%
    mutate(
      cells = cells_pos + cells_neg,
      log_density = log(density_pos + density_neg),
      density_ratio = if_else(density_ratio < min_density_ratio, min_density_ratio, density_ratio),
      density_ratio = if_else(density_ratio > max_density_ratio, max_density_ratio, density_ratio)
    ) %>%
    filter(cells > 0)

  ggplot(density_split_tbl, aes(x = x, y = y)) +
    # geom_tile(aes(fill = density_ratio)) +
    geom_raster(aes(fill = density_ratio), interpolate = interpolate) +
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    labs(title = glue("{group_pos} vs {group_neg}"), x = x_var, y = y_var) +
    scale_fill_gradient2(low = "#053061", mid = "gray80", high = "#E41A1C")

}

# violin plot split by specified group.
# default group is orig.ident
plot_violin <- function(metadata_tbl, color_scheme, y_var, x_var = "orig.ident") {

  violin_plot <- ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(x_var))) +
    geom_violin() +
    xlab(x_var) +
    ylab(y_var) +
    scale_fill_manual(values = color_scheme,
                      name = x_var) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(violin_plot)
}
