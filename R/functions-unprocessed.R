# # functions from the old cli script not yet converted to the package format
#
#

#

# # plot a set of genes
# plot_genes = function(seurat_obj, genes, filename_base) {
#
#   # color gradient for FeaturePlot-based plots
#   gradient_colors = c("gray85", "red2")
#
#   # switch to "RNA" assay from potentially "integrated"
#   DefaultAssay(seurat_obj) = "RNA"
#
#   # tSNE plots color-coded by expression level (should be square to match the original tSNE plots)
#   feat_plot =
#     FeaturePlot(
#       seurat_obj, features = genes, reduction = "tsne", cells = sample(colnames(seurat_obj)),
#       pt.size = 0.5, cols = gradient_colors, ncol = 4
#     )
#   ggsave(glue("{filename_base}.tsne.png"), plot = feat_plot, width = 16, height = 10, units = "in")
#   ggsave(glue("{filename_base}.tsne.pdf"), plot = feat_plot, width = 16, height = 10, units = "in")
#
#   # UMAP plots color-coded by expression level (should be square to match the original tSNE plots)
#   feat_plot =
#     FeaturePlot(
#       seurat_obj, features = genes, reduction = "umap", cells = sample(colnames(seurat_obj)),
#       pt.size = 0.5, cols = gradient_colors, ncol = 4
#     )
#   ggsave(glue("{filename_base}.umap.png"), plot = feat_plot, width = 16, height = 10, units = "in")
#   ggsave(glue("{filename_base}.umap.pdf"), plot = feat_plot, width = 16, height = 10, units = "in")
#
#   # dot plot visualization
#   dot_plot = DotPlot(seurat_obj, features = genes, dot.scale = 12, cols = gradient_colors)
#   ggsave(glue("{filename_base}.dotplot.png"), plot = dot_plot, width = 20, height = 8, units = "in")
#   ggsave(glue("{filename_base}.dotplot.pdf"), plot = dot_plot, width = 20, height = 8, units = "in")
#
#   # gene violin plots (size.use below 0.2 doesn't seem to make a difference)
#   # skip PDF since every cell has to be plotted and they become too big
#   vln_plot = VlnPlot(seurat_obj, features = genes, pt.size = 0.1, combine = TRUE, cols = colors_clusters, ncol = 4)
#   ggsave(glue("{filename_base}.violin.png"), plot = vln_plot, width = 16, height = 10, units = "in")
#
#   # expression levels per cluster for bar plots (averaging and output are in non-log space)
#   cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", features = genes, verbose = FALSE)[["RNA"]]
#   cluster_avg_exp_long = cluster_avg_exp %>% rownames_to_column("gene") %>% gather(cluster, avg_exp, -gene)
#
#   # bar plots
#   # create a named color scheme to ensure names and colors are in the proper order
#   clust_names = levels(seurat_obj)
#   color_scheme_named = colors_clusters[1:length(clust_names)]
#   names(color_scheme_named) = clust_names
#   barplot_plot = ggplot(cluster_avg_exp_long, aes(x = cluster, y = avg_exp, fill = cluster)) +
#     geom_col(color = "black") +
#     theme(legend.position = "none") +
#     scale_fill_manual(values = color_scheme_named) +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme_cowplot() +
#     facet_wrap(~ gene, ncol = 4, scales = "free")
#   ggsave(glue("{filename_base}.barplot.png"), plot = barplot_plot, width = 16, height = 10, units = "in")
#   ggsave(glue("{filename_base}.barplot.pdf"), plot = barplot_plot, width = 16, height = 10, units = "in")
#
# }
