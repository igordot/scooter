#' Scatter plots of the per-cell QC metrics against each other.
#'
#' One panel per pair in `pairs`, colored by sample. No legend, since this only
#' ever runs on a single sample - see [filter_cells()]. Pairs are a curated
#' list rather than every combination of a metric vector, since that would make
#' the panel count quadratic in the number of metrics (5 metrics -> 10 panels
#' via `combn()`) and most pairs are not biologically meaningful together.
#'
#' @param x Seurat object.
#' @param pairs A list of length-2 character vectors, `c(feature1, feature2)`,
#'   one per panel.
#' @param group_by Metadata column to color the points by.
#' @param color_scheme (optional) Named vector of colors. Defaults to the
#'   "samples" scheme, named by the sorted levels of `group_by`.
#' @param file Path to save the plot to. `NULL` returns it without writing.
#'   Parent directories are created as needed.
#' @param width,height Saved size in inches. The default is wide and short: one
#'   square panel per pair of metrics, side by side.
#'
#' @return A combined plot, one panel per pair.
#'
#' @importFrom cowplot plot_grid
#' @export
plot_metrics_correlations <- function(
  x,
  pairs = list(
    c("total_counts", "detected_genes"),
    c("pct_mito", "detected_genes"),
    c("pct_mito", "pct_ribo")
  ),
  group_by = "orig.ident",
  color_scheme = NULL,
  file = NULL,
  width = 16,
  height = 5
) {
  # validates that group_by/the pair features exist; the returned scheme itself
  # is unused below - points are always solid black regardless of sample, since
  # this only ever runs on one sample
  metrics_color_scheme(x, unique(unlist(pairs)), group_by, color_scheme)

  scatter_theme <-
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )

  scatter_plots <- lapply(pairs, function(feature_pair) {
    scatter_plot <- FeatureScatter(
      x,
      feature1 = feature_pair[1],
      feature2 = feature_pair[2],
      group.by = group_by,
      cols = "black"
    )
    scatter_plot$layers[[1]]$aes_params$alpha <- 0.1
    # red to match VariableFeaturePlot()'s own highlight color
    # (normalize_counts()'s variance-features.png), no confidence ribbon
    scatter_plot <- scatter_plot +
      geom_smooth(method = "gam", se = FALSE, color = "red", linewidth = 0.8)
    scatter_plot + scatter_theme
  })

  plot_obj <- plot_grid(plotlist = scatter_plots, ncol = length(scatter_plots))

  if (!is.null(file)) {
    ggsave(
      file,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      create.dir = TRUE
    )
  }

  plot_obj
}

#' Violin plots of the per-cell QC metrics.
#'
#' One violin panel per metric, grouped by sample. The metrics are metadata
#' columns, so `layer` is set explicitly: `VlnPlot()` otherwise searches for
#' the "data" layer, which does not exist before `NormalizeData()` has been
#' run.
#'
#' @param x Seurat object.
#' @param metrics Metadata columns to plot, one panel each. No default - the
#'   caller decides, since the right set depends on context (e.g. a single
#'   sample can afford `pct_ribo`/`pct_hb` alongside
#'   `detected_genes`/`total_counts`/`pct_mito`; a multi-sample plot usually
#'   can't). A `pct_*` metric that is under 1% in every cell is skipped (a flat
#'   near-zero violin says nothing) - its range is reported with a message
#'   instead.
#' @param group_by Metadata column to group the violins by.
#' @param color_scheme (optional) Named vector of colors. Defaults to the
#'   "samples" scheme, named by the sorted levels of `group_by`.
#' @param file Path to save the plot to. `NULL` returns it without writing.
#'   Parent directories are created as needed.
#' @param width,height Saved size in inches. Widen it when there are many
#'   samples.
#'
#' @return A combined plot, one panel per plotted metric.
#'
#' @importFrom cowplot plot_grid
#' @importFrom scales comma
#' @export
plot_metrics_distribution <- function(
  x,
  metrics,
  group_by = "orig.ident",
  color_scheme = NULL,
  file = NULL,
  width = 12,
  height = 6
) {
  color_scheme <- metrics_color_scheme(x, metrics, group_by, color_scheme)

  # a pct_* metric that never clears 1% plots as a flat, uninformative violin -
  # skip it and report the range instead
  is_negligible <- vapply(
    metrics,
    function(metric) {
      startsWith(metric, "pct_") && max(x@meta.data[[metric]]) < 1
    },
    logical(1)
  )
  for (metric in metrics[is_negligible]) {
    metric_range <- range(x@meta.data[[metric]])
    message(sprintf(
      "%s: %.3g-%.3g%% in every cell, skipping plot",
      metric,
      metric_range[1],
      metric_range[2]
    ))
  }
  metrics <- metrics[!is_negligible]

  vln_theme <-
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )

  vln_plot_grid <- suppressMessages({
    vln_plots <- lapply(metrics, function(qc_metric) {
      vln_plot <-
        VlnPlot(
          x,
          features = qc_metric,
          group.by = group_by,
          layer = "counts",
          pt.size = 0.1,
          sort = TRUE,
          combine = TRUE,
          cols = color_scheme
        ) +
        scale_y_continuous(labels = comma) +
        vln_theme
      vln_plot$layers[[2]]$aes_params$alpha <- 0.1
      vln_plot
    })
    plot_grid(plotlist = vln_plots, ncol = length(metrics))
  })

  if (!is.null(file)) {
    ggsave(
      file,
      plot = vln_plot_grid,
      width = width,
      height = height,
      units = "in",
      create.dir = TRUE
    )
  }

  vln_plot_grid
}

#' Euler diagram of the variable genes shared between batches.
#'
#' Becomes unreadable and can take a very long time to fit for many sets, so
#' returns `NULL` when there are too many batches.
#'
#' @param var_genes_list Named list of variable gene vectors, one element per
#'   batch.
#' @param color_scheme (optional) Vector of colors.
#' @param max_sets Return `NULL` when there are more sets than this.
#'
#' @return A plot, or `NULL` if there are too many sets to draw.
#'
#' @importFrom eulerr euler
#' @export
plot_var_genes_euler <- function(
  var_genes_list,
  color_scheme = NULL,
  max_sets = 8
) {
  if (length(var_genes_list) >= max_sets) {
    return(NULL)
  }

  if (is.null(color_scheme)) {
    color_scheme <- get_color_scheme("samples")[seq_along(var_genes_list)]
  }

  euler_fit <- euler(var_genes_list, shape = "ellipse")
  plot(
    euler_fit,
    fills = list(fill = color_scheme, alpha = 0.7),
    edges = list(col = color_scheme)
  )
}

#' UpSet plot of the variable genes shared between batches.
#'
#' @param var_genes_list Named list of variable gene vectors, one element per
#'   batch.
#' @param nsets Maximum number of sets to include.
#' @param nintersects Maximum number of intersections to show.
#'
#' @return An UpSet plot.
#'
#' @importFrom UpSetR upset fromList
#' @export
plot_var_genes_upset <- function(var_genes_list, nsets = 50, nintersects = 15) {
  # UpSetR's own internals still call the ggplot2 3.0-era aes_string() and the
  # pre-3.4 `size` aesthetic/argument - both deprecated upstream in ggplot2,
  # neither fixable from here
  suppressWarnings(upset(
    fromList(var_genes_list),
    nsets = nsets,
    nintersects = nintersects,
    order.by = "freq",
    mb.ratio = c(0.5, 0.5)
  ))
}

plot_qc <- function(
  x,
  proj_name = "",
  label = "",
  features = "detected_genes",
  grouping = "orig.ident",
  out_dir = "."
) {
  # create output directories
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  qc_dir <- glue("{out_dir}/qc")
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir)
  }

  qc_plot <- plot_distribution(
    x = x,
    features = c("detected_genes"),
    grouping = grouping,
    color_scheme = NULL
  )

  ggsave(filename = glue("{qc_dir}/{proj_name}{label}-qc.png"))
}

#' Check the QC metadata columns and build the default sample palette.
#'
#' Shared by the QC plots so they agree on both the error and the colors.
#'
#' @param x Seurat object.
#' @param features Metadata columns being plotted.
#' @param group_by Metadata column the plot groups by.
#' @param color_scheme (optional) Named vector of colors, returned unchanged
#'   when given.
#'
#' @return A named vector of colors, one per group.
#'
#' @noRd
metrics_color_scheme <- function(x, features, group_by, color_scheme = NULL) {
  missing_features <- setdiff(c(features, group_by), colnames(x@meta.data))
  if (length(missing_features)) {
    stop(
      "metadata columns not found: ",
      paste(missing_features, collapse = ", ")
    )
  }

  if (!is.null(color_scheme)) {
    return(color_scheme)
  }

  # named, to keep colors matched to groups regardless of the order they are
  # plotted in
  group_names <- x@meta.data[[group_by]] |>
    as.character() |>
    sort() |>
    unique()
  color_scheme <- get_color_scheme("samples")[seq_along(group_names)]
  names(color_scheme) <- group_names
  color_scheme
}

#' Width for a multi-sample QC violin plot, wide enough that individual violins
#' stay readable.
#'
#' Unlike [dr_plot_width()] (sized for a single scatter plot's legend, which is
#' fine to cap at a fixed width), a violin panel holds one violin per sample
#' side by side - the space each violin gets keeps shrinking as more samples
#' are added, so this scales in wider, uncapped steps than [dr_plot_width()]
#' rather than sharing its tiers.
#'
#' @param x A Seurat object, or a cell/sample count.
#'
#' @return A plot width in inches.
#'
#' @noRd
metrics_plot_width <- function(x) {
  num_samples <- if (is(x, "Seurat")) length(levels(x$orig.ident)) else x

  plot_width <- 12
  if (num_samples > 5) {
    plot_width <- 15
  }
  if (num_samples > 20) {
    plot_width <- 20
  }
  if (num_samples > 40) {
    plot_width <- 30
  }
  plot_width
}
