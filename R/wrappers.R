# Wrappers for pipeline

# load data
# create seurat obj
# plot qc

plot_qc <- function(seurat_obj, proj_name = "", label = "",features = "nFeature_RNA", grouping = "orig.ident", out_dir = ".") {

  # create output directories
  if (!dir.exists(out_dir)) dir.create(out_dir)
  qc_dir = glue("{out_dir}/qc")
  if (!dir.exists(qc_dir)) dir.create(qc_dir)

  qc_plot <- plot_distribution(data = seurat_obj,
                    features = c("nFeature_RNA"),
                    grouping = grouping,
                    color_scheme = NULL)

  ggsave(filename = glue("{qc_dir}/{proj_name}{label}.qc.png"))
}

plot_HTO <- function(seurat_obj, proj_name = "", label = "", out_dir = ".") {
  # create output directories
  if (!dir.exists(out_dir)) dir.create(out_dir)
  hto_dir = glue("{out_dir}/hto")
  if (!dir.exists(hto_dir)) dir.create(hto_dir)

  # plot standard seurat plots for HTO analysis
  Idents(seurat_obj) <- "HTO_maxID"

  seurat_obj_ridge_plot <- RidgePlot(seurat_obj,
                                     assay = "HTO",
                                     features = rownames(seurat_obj[["HTO"]]),
                                     ncol = 2)
  ggsave(seurat_obj_ridge_plot,
         filename = glue("{hto_dir}/{proj_name}{label}.htodemux_ridge_plot.png"),
         height = 7,
         width =7 ,
         units = "in")

  # heatmap
  seurat_obj_heatmap <- HTOHeatmap(seurat_obj,
                                   assay = "HTO",
                                   ncells = 800)
  ggsave(seurat_obj_heatmap,
         filename =  glue("{hto_dir}/{proj_name}{label}.htodemux_heatmap.png"),
         height = 7,
         width =7 ,
         units = "in")

  # scatter plots
  HTO_seurat <- as.data.frame(t(as.data.frame(seurat_obj@assays$HTO@data)))
  hto_pairs <- ggpairs(HTO_seurat)
  ggsave(hto_pairs,
         filename =  glue("{hto_dir}/{proj_name}{label}.htodemux_pairs.png"),
         height = 7,
         width =7 ,
         units = "in")

  # plot the doublet trends
  doublet_trend <- ggplot(as.data.frame(seurat_obj$HTO_classification),
                          aes(seurat_obj$HTO_classification)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(doublet_trend,
         filename = glue("{hto_dir}/{proj_name}{label}.htodemux_doubet.png"),
         height = 7,
         width =7 ,
         units = "in")

  doub_col <- c("red", "blue", "green")
  names(doub_col) <- c("Singlet", "Negative", "Doublet")
  doublet_bar <- data.frame(doub = seurat_obj$HTO_classification.global)

  doublet_bar <- doublet_bar %>%
    group_by(doub) %>%
    summarise (n = n()) %>%
    unique() %>%
    mutate(percentage = n /sum(n))

  doublet_bar <- doublet_bar %>%
    mutate(sample = rep(proj_name, nrow(doublet_bar)))

  plot_num_doublet <- ggplot(doublet_bar, aes(x = sample, y =percentage, fill = doub))+
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = doub_col) +
    theme_bw()

    ggsave(plot_num_doublet,
           filename =  glue("{hto_dir}/{proj_name}{label}.htodemux_doublet_count.png"),
           height = 7,
           width =7 ,
           units = "in")

    # ridge plots
    Idents(seurat_obj) <- "MULTI_ID"
    seurat_obj_ridge_plot <- RidgePlot(seurat_obj,
                                       assay = "HTO",
                                       features = rownames(seurat_obj[["HTO"]]),
                                       ncol = 2)

    ggsave(seurat_obj_ridge_plot,
           filename = glue("{hto_dir}/{proj_name}{label}.htomulti_ridge_plot.png"),
           height = 7,
           width =7 ,
           units = "in")

    # heatmap
    seurat_obj_heatmap <- HTOHeatmap(seurat_obj, assay = "HTO", ncells = 800)
    ggsave(seurat_obj_heatmap,
           filename =  glue("{hto_dir}/{proj_name}{label}.htomulti_heatmap.png"),
           height = 7,
           width =7 ,
           units = "in")

    # scatter plots
    HTO_seurat <- as.data.frame(t(as.data.frame(seurat_obj@assays$HTO@data)))
    hto_pairs <- ggpairs(HTO_seurat)

    ggsave(hto_pairs, filename =  glue("{hto_dir}/{proj_name}{label}.htomulti_pairs.png"),
           height = 7,
           width =7 ,
           units = "in")

    # plot the doublet trends
    doublet_trend <- ggplot(as.data.frame(seurat_obj$MULTI_classification),
                            aes(seurat_obj$MULTI_classification)) +
      geom_bar() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(doublet_trend,
           filename = glue("{hto_dir}/{proj_name}{label}.htomulti_doubet.png"),
           height = 7,
           width =7 ,
           units = "in")
}




