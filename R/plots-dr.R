#' Plot a single feature overlaid on a dimensionality reduction.
#'
#' One point per cell, colored by expression/abundance of `feature` (a gene,
#' antibody, or any other value `FeaturePlot()` can fetch). Unlike passing a
#' multi-color ramp straight to `FeaturePlot(cols = )`, the ramp is applied
#' afterward with `scale_color_gradientn()` so the color scale stays
#' continuous - see the "cols" note below.
#'
#' @param x Seurat object.
#' @param feature Feature name (gene, antibody, or metadata column) to plot.
#' @param assay (optional) Assay to fetch `feature` from. `NULL` uses
#'   `DefaultAssay(x)`.
#' @param reduction Dimensionality reduction to plot on (e.g. "umap", "tsne",
#'   "pca").
#' @param cells (optional) Cells to plot, in the given order (controls z-order
#'   of overlapping points). Defaults to a random shuffle of all cells so no
#'   single cell/group dominates due to plot order.
#' @param color_scheme (optional) Vector of colors for the low-to-high
#'   expression gradient.
#' @param pt_size (optional) Point size. Defaults to `get_dr_point_size(x)`.
#'
#' @return A ggplot object.
#'
#' @section Seurat's `FeaturePlot(cols = )` quirk:
#' As of Seurat 5, passing a `cols` vector of anything other than exactly 2
#' colors makes `FeaturePlot()` bin the continuous values into 2 groups before
#' coloring (`cut(x, breaks = 2)`), which turns a smooth expression gradient
#' into a flat low/high split. Seurat 4 only did this when `cols` had exactly 2
#' colors, so a long custom ramp used to pass through untouched. This function
#' sidesteps the quirk (and the version difference) entirely by never passing a
#' multi-color `cols` to `FeaturePlot()` - the ramp is layered on afterward
#' instead.
#'
#' @export
plot_dr_feature <- function(
  x,
  feature,
  assay = NULL,
  reduction = "umap",
  cells = NULL,
  color_scheme = NULL,
  pt_size = NULL
) {
  if (is.null(cells)) {
    cells <- sample(colnames(x))
  }
  if (is.null(color_scheme)) {
    color_scheme <- colorRampPalette(c(
      "#d9cfcb",
      "#d49070",
      "#ca5528",
      "#b72600",
      "#981000",
      "#730000"
    ))(100)
  }
  if (is.null(pt_size)) {
    pt_size <- get_dr_point_size(x)
  }

  feat_plot <- FeaturePlot(
    x,
    features = feature,
    assay = assay,
    reduction = reduction,
    cells = cells,
    pt.size = pt_size,
    raster = FALSE
  )
  # FeaturePlot() always attaches its own color scale, so replacing it with
  # color_scheme (see the "cols" quirk note above) always triggers ggplot2's
  # scale-collision message
  suppressMessages(feat_plot + scale_color_gradientn(colors = color_scheme)) +
    labs(title = feature) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
}

#' Scatter plot of a reduction, colored by a grouping variable.
#'
#' The categorical counterpart of [plot_dr_feature()]: cells are colored by a
#' metadata column (or by the current identities) rather than by a continuous
#' feature. Cells are shuffled so no one group is drawn on top of the others.
#'
#' @param x Seurat object.
#' @param reduction Dimensionality reduction to plot on (e.g. "umap", "tsne",
#'   "pca").
#' @param group_by (optional) Metadata column to color by. `NULL` uses the
#'   current identities.
#' @param color_scheme (optional) Vector of colors. Defaults to the "clusters"
#'   scheme.
#' @param pt_size (optional) Point size. Defaults to `get_dr_point_size(x)`.
#' @param na_value Color for cells with no value in `group_by`.
#' @param file_prefix (optional) Path to save to.
#' @param file_format File formats to write when `file_prefix` is given.
#' @param width (optional) Saved width in inches. `NULL` measures the plot's
#'   own legend and widens to fit it - see [save_dr_plot()].
#' @param height Saved height in inches.
#'
#' @return A ggplot object.
#'
#' @export
plot_dr_group <- function(
  x,
  reduction = "umap",
  group_by = NULL,
  color_scheme = NULL,
  pt_size = NULL,
  na_value = "grey50",
  file_prefix = NULL,
  file_format = c("png", "pdf"),
  width = NULL,
  height = 6
) {
  if (is.null(color_scheme)) {
    color_scheme <- get_color_scheme("clusters")
  }
  if (is.null(pt_size)) {
    pt_size <- get_dr_point_size(x)
  }

  plot_obj <-
    DimPlot(
      x,
      reduction = reduction,
      group.by = group_by,
      pt.size = pt_size,
      cols = color_scheme,
      shuffle = TRUE,
      na.value = na_value,
      raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  if (!is.null(file_prefix)) {
    save_dr_plot(
      plot_obj,
      file_prefix = file_prefix,
      file_format = file_format,
      width = width,
      height = height
    )
  }

  plot_obj
}

#' Save a reduction plot, widened to fit its own legend.
#'
#' Unlike [dr_plot_width()] (a label-count guess made before the plot exists),
#' this measures the legend actually rendered on `plot` - so long label text
#' widens the saved image even when the label count is small. Requires an
#' open graphics device to get accurate text metrics, so it opens and closes
#' a throwaway one around the measurement.
#'
#' @param plot A ggplot object, such as from [plot_dr_group()]. Forwarded
#'   straight to `ggsave(plot = )`.
#' @param file_prefix Path to save to without an extension.
#' @param file_format File formats to write.
#' @param width (optional) Saved width in inches. `NULL` measures `plot`'s
#'   legend and adds it to `panel_width`.
#' @param height Saved height in inches.
#' @param panel_width Width of the plot panel itself, excluding the legend, in
#'   inches. Only used when `width` is `NULL`.
#'
#' @return `plot`, invisibly.
#'
#' @importFrom cowplot get_legend
#' @importFrom grDevices dev.off pdf
#' @export
save_dr_plot <- function(
  plot,
  file_prefix,
  file_format = c("png", "pdf"),
  width = NULL,
  height = 6,
  panel_width = height
) {
  if (is.null(width)) {
    # a plot with no legend (e.g. a single-color scatter) has no "guide-box"
    # grob to measure at all, not a zero-width one
    legend_width <- tryCatch(
      {
        legend_grob <- get_legend(plot)
        pdf(NULL)
        on.exit(dev.off())
        grid::convertWidth(grid::grobWidth(legend_grob), "in", valueOnly = TRUE)
      },
      error = function(e) 0
    )
    width <- panel_width + legend_width + 0.5
  }

  for (format in file_format) {
    ggsave(
      glue("{file_prefix}.{format}"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      create.dir = TRUE
    )
    # deliberate: one second between writes keeps the modification times
    # distinct, so the output files sort in the order they were produced. Do
    # not remove as an optimization.
    Sys.sleep(1)
  }

  invisible(plot)
}

# plot colored by specified variable
plot_scatter_group <- function(
  metadata_tbl,
  x_var = "UMAP_1",
  y_var = "UMAP_2",
  aspect_ratio = 1,
  color_var,
  color_scheme
) {
  ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(
      aes(color = !!sym(color_var)),
      size = get_dr_point_size(metadata_tbl)
    ) +
    theme(
      aspect.ratio = aspect_ratio,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_scheme)
}

# plot split by specified variable
plot_scatter_split <- function(
  metadata_tbl,
  x_var = "UMAP_1",
  y_var = "UMAP_2",
  aspect_ratio = 1,
  rows_var = NULL,
  cols_var = NULL,
  color_var,
  color_scheme
) {
  gp <-
    ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(
      aes(color = !!sym(color_var)),
      size = get_dr_point_size(metadata_tbl)
    ) +
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

# density plot split by specified variable calculate density normalized to 1,
# independently for each facet variable
plot_density_split <- function(
  metadata_tbl,
  x_var,
  y_var,
  split_var,
  num_bins
) {
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
get_density_diff_table <- function(
  metadata_tbl,
  x_var,
  y_var,
  split_var,
  num_bins
) {
  # generate a density plot split by stage
  density_plot <- plot_density_split(
    metadata_tbl = metadata_tbl,
    x_var = x_var,
    y_var = y_var,
    split_var = split_var,
    num_bins = num_bins
  )

  # produce an object that can be rendered
  density_plot_tbl <- ggplot_build(density_plot)

  # panel labels
  panels_tbl <-
    tibble(
      PANEL = density_plot_tbl$layout$layout$PANEL,
      stage = density_plot_tbl$layout$layout$stage
    )

  # merge panel contents and panel names
  density_tbl <- density_plot_tbl$data[[1]]
  density_tbl <- density_tbl |> full_join(panels_tbl, by = "PANEL")

  return(density_tbl)
}

# density plot split by specified variable split normalization (adding
# norm_split_var) may not work
plot_density_diff <- function(
  metadata_tbl,
  x_var = "UMAP_1",
  y_var = "UMAP_2",
  split_var,
  num_bins,
  group_pos,
  group_neg,
  interpolate = FALSE
) {
  density_tbl <- get_density_diff_table(
    metadata_tbl = metadata_tbl,
    x_var = x_var,
    y_var = y_var,
    split_var = split_var,
    num_bins = num_bins
  )

  min_density <- quantile(density_tbl$density, 0)

  density_pos_tbl <-
    density_tbl |>
    filter(stage == group_pos) |>
    select(x, y, cells_pos = count, density_pos = density)
  density_neg_tbl <-
    density_tbl |>
    filter(stage == group_neg) |>
    select(x, y, cells_neg = count, density_neg = density)
  density_split_tbl <- full_join(
    density_pos_tbl,
    density_neg_tbl,
    by = c("x", "y")
  )
  density_split_tbl[is.na(density_split_tbl)] <- min_density
  density_split_tbl <- density_split_tbl |>
    mutate(density_diff = density_pos - density_neg)
  density_split_tbl <- density_split_tbl |>
    mutate(density_ratio = log(density_pos / density_neg))

  min_density_diff <- density_split_tbl |> pull(density_diff) |> quantile(0.01)
  max_density_diff <- density_split_tbl |> pull(density_diff) |> quantile(0.99)
  min_density_ratio <- density_split_tbl |>
    pull(density_ratio) |>
    quantile(0.01)
  max_density_ratio <- density_split_tbl |>
    pull(density_ratio) |>
    quantile(0.99)

  density_split_tbl <-
    density_split_tbl |>
    mutate(
      cells = cells_pos + cells_neg,
      log_density = log(density_pos + density_neg),
      density_ratio = if_else(
        density_ratio < min_density_ratio,
        min_density_ratio,
        density_ratio
      ),
      density_ratio = if_else(
        density_ratio > max_density_ratio,
        max_density_ratio,
        density_ratio
      )
    ) |>
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

# violin plot split by specified group. default group is orig.ident
plot_violin <- function(
  metadata_tbl,
  color_scheme,
  y_var,
  x_var = "orig.ident"
) {
  violin_plot <- ggplot(
    metadata_tbl,
    aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(x_var))
  ) +
    geom_violin() +
    xlab(x_var) +
    ylab(y_var) +
    scale_fill_manual(values = color_scheme, name = x_var) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(violin_plot)
}

#' Determine the color scheme.
#'
#' Determine the color scheme. Can be specified for samples or clusters to
#' avoid confusion.
#'
#' @param type Type of scheme ("samples" or "clusters").
#'
#' @return A vector of colors.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggsci pal_d3 pal_igv
#' @export
get_color_scheme <- function(type = "clusters") {
  color_scheme <- c(
    pal_igv("default")(51),
    pal_igv(alpha = 0.6)(51),
    pal_igv(alpha = 0.3)(51)
  )
  if (type == "samples") {
    color_scheme <- c(
      brewer.pal(5, "Set1"),
      brewer.pal(8, "Dark2"),
      color_scheme
    )
  }
  if (type == "clusters") {
    color_scheme <- c(
      pal_d3("category10")(10),
      pal_d3("category20b")(20),
      color_scheme
    )
  }

  return(color_scheme)
}

#' Width for a reduction plot, widened so many legend labels stay legible.
#'
#' A label-count estimate made without building the plot. [plot_dr_group()]
#' now sizes itself from the actual rendered legend instead (see
#' [save_dr_plot()]), so this has no callers left in the package - kept as a
#' cheap standalone estimate for callers that don't have a plot object yet.
#'
#' @param x Seurat object, or the number of legend labels.
#' @param group_by Metadata column whose distinct values become the plot's
#'   legend. Ignored if `x` is already a number.
#'
#' @return Numeric width in inches.
#'
#' @export
dr_plot_width <- function(x, group_by = "orig.ident") {
  num_labels <- if (is(x, "Seurat")) {
    group_by <- check_identity_column(x, group_by)
    x@meta.data[[group_by]] |> as.character() |> n_distinct()
  } else {
    x
  }

  plot_width <- 12
  if (num_labels > 15) {
    plot_width <- 15
  }
  if (num_labels > 30) {
    plot_width <- 20
  }
  if (num_labels > 50) {
    plot_width <- 30
  }
  plot_width
}

#' Determine the point size for reduced dimensions scatter plots (smaller for
#' larger datasets).
#'
#' @param x Number of cells (points on the plot), a Seurat object, or a data
#'   frame with one row per cell.
#'
#' @return Numeric point size.
#' @export
get_dr_point_size <- function(x) {
  num_cells <- if (is(x, "Seurat")) {
    ncol(x)
  } else if (is(x, "data.frame")) {
    nrow(x)
  } else {
    x
  }

  pt_size <- 1.8
  if (num_cells > 1000) {
    pt_size <- 1.4
  }
  if (num_cells > 5000) {
    pt_size <- 1.0
  }
  if (num_cells > 10000) {
    pt_size <- 0.6
  }
  if (num_cells > 50000) {
    pt_size <- 0.2
  }
  return(pt_size)
}

#' Plot a UMAP colored by one or more cluster/grouping variables.
#'
#' @param x Seurat object with a "umap" reduction.
#' @param cluster_var Metadata column(s) to color by, as a single string
#'   (comma-separated for more than one).
#' @param color_scheme (optional) Named vector of colors. Defaults to
#'   `get_color_scheme("clusters")`.
#'
#' @return Nothing; called for its file output, one `dr-umap-<cluster_var>.png`
#'   per variable.
#'
#' @importFrom cowplot save_plot
#' @importFrom stringr str_c
#' @export
plot_dr_umap_clusters <- function(
  x,
  cluster_var,
  color_scheme = get_color_scheme("clusters")
) {
  if (is.null(x@reductions$umap)) {
    stop("UMAP not computed yet")
  }

  cluster_vars <- unlist(strsplit(cluster_var, ",", fixed = TRUE))
  dr_point_size <- get_dr_point_size(x)
  set.seed(99)
  shuffled_cells <- sample(colnames(x))

  for (c_var in cluster_vars) {
    message("")
    x <- set_identity(x = x, c_var)
    grouping_label <- check_identity_column(x, c_var)
    grouping_label <- gsub("\\.", "", grouping_label)
    message(grouping_label, " : ", str_c(levels(x), collapse = ", "))

    num_clusters <- Idents(x) |> as.character() |> n_distinct()
    plot_width <- 12
    if (num_clusters > 15) {
      plot_width <- 15
    }
    if (num_clusters > 30) {
      plot_width <- 20
    }
    if (num_clusters > 50) {
      plot_width <- 30
    }

    cluster_umap <-
      DimPlot(
        x,
        reduction = "umap",
        cells = shuffled_cells,
        pt.size = dr_point_size,
        cols = color_scheme,
        na.value = "grey90",
        raster = FALSE
      ) +
      theme(
        plot.background = element_rect(fill = "white"),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )
    save_plot(
      glue("dr-umap-{grouping_label}.png"),
      plot = cluster_umap,
      base_height = 6,
      base_width = plot_width
    )
    Sys.sleep(1)
  }
}
