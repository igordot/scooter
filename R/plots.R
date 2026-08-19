#' Plot the distribution of specified features/variables.
#'
#' @param x Seurat object or metadata.
#' @param features Vector of features to plot (such as "nGene", "nUMI",
#'   "percent.mito").
#' @param grouping X.
#' @param color_scheme (optional) Vector of colors.
#'
#' @return A vector of colors.
#'
#' @importFrom data.table :=
#' @importFrom rlang sym syms
#' @export
plot_distribution <- function(x, features, grouping, color_scheme = NULL) {
  UseMethod("plot_distribution")
}

#' @export
plot_distribution.default <- function(
  x,
  features,
  grouping,
  color_scheme = NULL
) {
  if (is.null(color_scheme)) {
    color_scheme <- get_color_scheme()
  }

  dist_tbl <- x |>
    select(c(features, grouping)) |>
    tidyr::gather(key = "var", value = "val", -grouping)

  plot_dist <-
    ggplot(dist_tbl, aes(x = !!sym(grouping), y = val)) +
    cowplot::theme_cowplot() +
    geom_jitter(size = 1, color = "black", alpha = 0.1, width = 0.3) +
    geom_violin(aes(fill = !!sym(grouping))) +
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
plot_distribution.Seurat <- function(
  x,
  features,
  grouping,
  color_scheme = NULL
) {
  # compile the data table
  dist_tbl <- Seurat::FetchData(object = x, vars = c(features, grouping))
  plot_distribution(dist_tbl, features, grouping, color_scheme)
}

plot_scatter <- function(
  metadata,
  out_path = NULL,
  proj_name = NULL,
  log_file = NULL,
  X,
  Y,
  color,
  write = FALSE,
  color_vect = NULL
) {
  if (is.factor(metadata[color][[1]]) | is.character(metadata[color][[1]])) {
    # Create color vector if not supplied
    if (is.null(color_vect)) {
      colors_samples_named <- create_color_vect(as.data.frame(metadata[color]))
    } else {
      colors_samples_named <- color_vect
    }
    current_plot <- ggplot(
      sample_frac(metadata),
      aes(
        x = !!sym(X),
        y = !!sym(Y),
        color = !!sym(color)
      )
    ) +
      geom_point(size = 1, alpha = 0.7) +
      # coord_fixed(ratio = (max(metadata[,which(colnames(metadata) == X)]) - min(metadata[,which(colnames(metadata) == X)]))/(max(metadata[,which(colnames(metadata) == Y)]) - min(metadata[,which(colnames(metadata) == Y)]))) +
      xlab(X) +
      ylab(Y) +
      scale_color_manual(
        values = colors_samples_named,
        name = color
      )
  } else {
    current_plot <- ggplot(
      sample_frac(metadata),
      aes(
        x = !!sym(X),
        y = !!sym(Y),
        color = !!sym(color)
      )
    ) +
      geom_point(size = 1, alpha = 0.7) +
      # coord_fixed(ratio = (max(metadata[,which(colnames(metadata) == X)]) - min(metadata[,which(colnames(metadata) == X)]))/(max(metadata[,which(colnames(metadata) == Y)]) - min(metadata[,which(colnames(metadata) == Y)]))) +
      xlab(X) +
      ylab(Y) +
      ggsci::scale_color_material("purple")
  }

  #
  if (is.null(proj_name)) {
    proj_name <- ""
  } else {
    proj_name <- paste0(proj_name, ".")
  }

  if (write) {
    ggsave(
      filename = glue("{out_path}/{proj_name}{X}-{Y}-{color}.png"),
      plot = current_plot,
      width = 8,
      height = 6,
      units = "in"
    )
  }

  return(current_plot)
}

cluster_stats_bar <- function(
  metadata,
  group1,
  group2,
  out_path = ".",
  write = FALSE,
  g1_col = NULL,
  g2_col = NULL,
  cluster = TRUE
) {
  # make barplots and output cluster stats
  summary_metadata <- metadata |>
    group_by(!!!syms(c(group1, group2))) |>
    summarize(n = n()) |>
    group_by(!!sym(group1)) |>
    mutate(pct_g2_in_g1 = n / sum(n)) |>
    group_by(!!sym(group2)) |>
    mutate(pct_g1_in_g2 = n / sum(n)) |>
    ungroup()

  if (write) {
    write_excel_csv(
      summary_metadata,
      file = glue("{out_path}/summary-{group1}-{group2}.csv")
    )

    # make both grouping variables factors
    summary_metadata <- summary_metadata |>
      mutate(
        !!group1 := as.factor(!!sym(group1)),
        !!group2 := as.factor(!!sym(group2))
      )
    if (cluster) {
      mat_g1 <- summary_metadata |>
        select(!!c(group1, group2, "pct_g1_in_g2")) |>
        spread(group2, "pct_g1_in_g2", fill = 0) |>
        as.data.frame() |>
        column_to_rownames(group1) |>
        as.matrix()

      hc_g1 <- hclust(dist(mat_g1), method = "ward.D2") # clusters rows of mat
      levels_g1 <- rownames(mat_g1)[order.dendrogram(as.dendrogram(hc_g1))]

      summary_metadata <- summary_metadata |>
        mutate(!!group1 := fct_relevel(!!sym(group1), levels_g1))

      mat_g2 <- summary_metadata |>
        select(!!c(group1, group2, "pct_g2_in_g1")) |>
        spread(group1, "pct_g2_in_g1", fill = 0) |>
        as.data.frame() |>
        column_to_rownames(group2) |>
        as.matrix()

      hc_g2 <- hclust(dist(mat_g2), method = "ward.D2") # clusters rows of mat
      levels_g2 <- rownames(mat_g2)[order.dendrogram(as.dendrogram(hc_g2))]

      summary_metadata <- summary_metadata |>
        mutate(!!group2 := fct_relevel(!!sym(group2), levels_g2))
    }
    # use levels to re-order factor
    if (is.null(g1_col)) {
      group1_col <- create_color_vect(as.data.frame(summary_metadata[group1]))
    } else {
      group1_col <- g1_col
    }
    if (is.null(g2_col)) {
      group2_col <- create_color_vect(as.data.frame(summary_metadata[group2]))
    } else {
      group2_col <- g2_col
    }

    summary_plots_g2 <- ggplot(summary_metadata) +
      geom_col(aes_string(x = group2, y = "pct_g1_in_g2", fill = group1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(
        values = group1_col,
        name = group1
      ) +
      ylab(glue("percent {group1} in {group2}"))

    summary_plots_g2_legend <- get_legend(summary_plots_g2)

    summary_plots_g1 <- ggplot(summary_metadata) +
      geom_col(aes_string(x = group1, y = "pct_g2_in_g1", fill = group2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(
        values = group2_col,
        name = group2
      ) +
      ylab(glue("percent {group2} in {group1}"))

    summary_plots_g1_legend <- get_legend(summary_plots_g1)

    summary_plots <- plot_grid(
      summary_plots_g2 + theme(legend.position = "none"),
      summary_plots_g2_legend,
      summary_plots_g1 + theme(legend.position = "none"),
      summary_plots_g1_legend
    )

    if (write) {
      ggsave(
        summary_plots,
        filename = glue("{out_path}/{group1}-{group2}-bar.png"),
        height = 10,
        width = 10
      )
    }
  }
  return(list(
    summary_metadata = summary_metadata,
    summary_plots = summary_plots
  ))
}

#' Plot a gene list in multiple ways: UMAP, dot plot, violin, and per-cluster
#' bar plot.
#'
#' @param x Seurat object with an identity set.
#' @param genes Genes to plot.
#' @param filename_base Path prefix for the output files.
#' @param color_scheme (optional) Named vector of cluster colors, used by the
#'   violin and bar plots. Defaults to `get_color_scheme("clusters")`.
#'
#' @return Nothing; called for its file output.
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom tidyr gather
#' @export
plot_cluster_markers_top <- function(
  x,
  genes,
  filename_base,
  color_scheme = get_color_scheme("clusters")
) {
  gradient_colors <- colorRampPalette(c(
    "#d9cfcb",
    "#d49070",
    "#ca5528",
    "#b72600",
    "#981000",
    "#730000"
  ))(100)

  # switch to "RNA" assay from potentially "integrated"
  DefaultAssay(x) <- "RNA"

  # UMAP plots color-coded by expression level (should be square to match the
  # original tSNE plots)
  feat_plot <- FeaturePlot(
    x,
    features = genes,
    reduction = "umap",
    cells = sample(colnames(x)),
    pt.size = 0.5,
    ncol = 4
  )
  # FeaturePlot() always attaches its own color scale, so replacing it with
  # gradient_colors always triggers ggplot2's scale-collision message
  feat_plot <- suppressMessages(
    feat_plot & scale_color_gradientn(colors = gradient_colors)
  )
  ggsave(
    glue("{filename_base}-umap.png"),
    plot = feat_plot,
    width = 16,
    height = 10,
    units = "in",
    dpi = 150
  )

  # dot plot visualization
  dot_plot <- DotPlot(
    x,
    features = genes,
    dot.scale = 12,
    cols = gradient_colors
  )
  ggsave(
    glue("{filename_base}-dotplot.png"),
    plot = dot_plot,
    width = 20,
    height = 8,
    units = "in"
  )
  ggsave(
    glue("{filename_base}-dotplot.pdf"),
    plot = dot_plot,
    width = 20,
    height = 8,
    units = "in"
  )

  # gene violin plots (size.use below 0.2 doesn't seem to make a difference)
  # skip PDF since every cell has to be plotted and they become too big
  vln_plot <- VlnPlot(
    x,
    features = genes,
    pt.size = 0.1,
    combine = TRUE,
    cols = color_scheme,
    ncol = 4
  )
  ggsave(
    glue("{filename_base}-violin.png"),
    plot = vln_plot,
    width = 16,
    height = 10,
    units = "in",
    dpi = 150
  )

  # expression levels per cluster for bar plots (averaging and output are in
  # non-log space)
  cluster_avg_exp <- AverageExpression(
    x,
    assay = "RNA",
    features = genes,
    verbose = FALSE
  )[["RNA"]]
  cluster_avg_exp_long <- cluster_avg_exp |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    gather(cluster, avg_exp, -gene)

  # bar plots - a named color scheme keeps names and colors in the proper order
  clust_names <- levels(x)
  color_scheme_named <- color_scheme[seq_along(clust_names)]
  names(color_scheme_named) <- clust_names
  barplot_plot <- ggplot(
    cluster_avg_exp_long,
    aes(x = cluster, y = avg_exp, fill = cluster)
  ) +
    geom_col(color = "black") +
    theme(
      plot.background = element_rect(fill = "white"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = color_scheme_named) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot() +
    facet_wrap(~gene, ncol = 4, scales = "free")
  ggsave(
    glue("{filename_base}-barplot.png"),
    plot = barplot_plot,
    width = 20,
    height = 10,
    units = "in",
    dpi = 150
  )
}

#' Heatmap of the top cluster markers, at several top-N cutoffs.
#'
#' @param x Seurat object with an identity set.
#' @param markers_tbl Marker table, such as from [calculate_cluster_markers()] -
#'   either the standard `log2FC`/`cluster` columns, or the pairwise
#'   `log2FC_min`/`cluster` columns.
#' @param num_genes Top-N genes per cluster to include, one heatmap per value.
#' @param filename_base Path prefix for the output files.
#'
#' @return Nothing; called for its file output.
#'
#' @importFrom grDevices dev.off pdf png
#' @importFrom pheatmap pheatmap
#' @importFrom readr write_csv
#' @importFrom stats quantile
#' @export
plot_cluster_markers_heatmap <- function(
  x,
  markers_tbl,
  num_genes,
  filename_base
) {
  # adjust pairwise clusters to match the standard format
  if ("log2FC_min" %in% colnames(markers_tbl)) {
    markers_tbl <- markers_tbl |> mutate(log2FC = log2FC_min)
  }

  # keep only the top cluster for each gene so each gene appears once
  markers_tbl <- markers_tbl |> filter(log2FC > 0)
  markers_tbl <- markers_tbl |>
    group_by(gene) |>
    top_n(1, log2FC) |>
    slice(1) |>
    ungroup()

  num_clusters <- Idents(x) |> as.character() |> n_distinct()
  marker_genes <- markers_tbl |> pull(gene) |> unique() |> sort()
  cluster_avg_exp <- AverageExpression(
    x,
    assay = "RNA",
    features = marker_genes,
    verbose = FALSE
  )[["RNA"]]
  cluster_avg_exp <- cluster_avg_exp |> as.matrix() |> log1p()
  cluster_avg_exp <- cluster_avg_exp[rowSums(cluster_avg_exp) > 0, ]

  # heatmap settings
  hm_colors <- colorRampPalette(c("#053061", "#FFFFFF", "#E41A1C"))(51)
  hm_width <- (num_clusters / 2) + 2

  for (ng in num_genes) {
    hm_base <- glue("{filename_base}-heatmap-top{ng}")

    markers_top_tbl <- markers_tbl |>
      group_by(cluster) |>
      top_n(ng, log2FC) |>
      ungroup()
    markers_top_tbl <- markers_top_tbl |> arrange(cluster, -log2FC)

    # generate the scaled expression matrix and save the text version
    hm_mat <- cluster_avg_exp[markers_top_tbl$gene, ]
    hm_mat <- hm_mat |> t() |> scale() |> t()
    hm_mat |>
      round(3) |>
      as_tibble(rownames = "gene") |>
      write_csv(glue("{hm_base}.csv"))
    Sys.sleep(1)

    # set outliers to 95th percentile to yield a more balanced color scale
    scale_cutoff <- as.numeric(quantile(abs(hm_mat), 0.95))
    hm_mat[hm_mat > scale_cutoff] <- scale_cutoff
    hm_mat[hm_mat < -scale_cutoff] <- -scale_cutoff

    # generate the heatmap
    ph_obj <- pheatmap(
      hm_mat,
      scale = "none",
      color = hm_colors,
      border_color = NA,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 12,
      show_colnames = TRUE,
      main = glue("Cluster Markers: Top {ng}")
    )

    png(
      glue("{hm_base}.png"),
      width = hm_width,
      height = 10,
      units = "in",
      res = 300
    )
    grid::grid.newpage()
    grid::grid.draw(ph_obj$gtable)
    dev.off()
    Sys.sleep(1)
    pdf(glue("{hm_base}.pdf"), width = hm_width, height = 10)
    grid::grid.newpage()
    grid::grid.draw(ph_obj$gtable)
    dev.off()
    Sys.sleep(1)
  }
}

#' Plot a cluster resolution on tSNE and UMAP.
#'
#' A single cluster says nothing, and past the color scheme's length the colors
#' start repeating, so nothing is plotted outside `2:length(color_scheme)`
#' clusters.
#'
#' @param x Seurat object.
#' @param resolution Metadata column (or bare resolution value, resolved via
#'   [check_identity_column()]) to set as the identity before plotting.
#' @param filename_base Path prefix for the two plots (`-tsne`/`-umap` are
#'   appended by [plot_dr_group()]).
#' @param color_scheme (optional) Named vector of cluster colors. Defaults to
#'   `get_color_scheme("clusters")`.
#'
#' @return `x` with `resolution` set as its identity.
#'
#' @export
plot_clusters <- function(
  x,
  resolution,
  filename_base,
  color_scheme = get_color_scheme("clusters")
) {
  x <- set_identity(x = x, identity_column = resolution)

  num_clusters <- Idents(x) |> as.character() |> n_distinct()
  message("resolution: ", resolution)
  message("num clusters: ", num_clusters)

  if (num_clusters > 1 && num_clusters < length(color_scheme)) {
    for (reduction in c("tsne", "umap")) {
      plot_dr_group(
        x,
        reduction = reduction,
        color_scheme = color_scheme,
        na_value = "grey90",
        file_prefix = glue("{filename_base}-{reduction}")
      )
    }
    if (file.exists("Rplots.pdf")) {
      file.remove("Rplots.pdf")
    }
  }

  x
}
