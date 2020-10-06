# # functions from the old cli script not yet converted to the package format
#
#

#

# # plot ADT metrics
# plot_adt_qc = function(seurat_obj) {
#
#   # create a named color scheme to ensure names and colors are in the proper order
#   sample_names = seurat_obj$orig.ident %>% as.character() %>% sort() %>% unique()
#   colors_samples_named = colors_samples[1:length(sample_names)]
#   names(colors_samples_named) = sample_names
#
#   vln_theme =
#     theme(
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "none"
#     )
#   suppressMessages({
#     dist_adt_g_plot =
#       VlnPlot(
#         seurat_obj, features = "num_ADT_genes", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_adt_u_plot =
#       VlnPlot(
#         seurat_obj, features = "num_ADT_UMIs", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_pct_m_plot =
#       VlnPlot(
#         seurat_obj, features = "pct_mito", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_plot = plot_grid(dist_adt_g_plot, dist_adt_u_plot, dist_pct_m_plot, ncol = 3)
#     ggsave("qc.adt.distribution.png", plot = dist_plot, width = 20, height = 6, units = "in")
#   })
#   Sys.sleep(1)
#
#   cor_umis_plot =
#     FeatureScatter(
#       seurat_obj, feature1 = "num_genes", feature2 = "num_ADT_genes",
#       group.by = "orig.ident", cols = colors_samples_named
#     ) +
#     theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
#   cor_genes_plot =
#     FeatureScatter(
#       seurat_obj, feature1 = "num_UMIs", feature2 = "num_ADT_UMIs",
#       group.by = "orig.ident", cols = colors_samples_named
#     ) +
#     theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
#   cor_mito_plot =
#     FeatureScatter(
#       seurat_obj, feature1 = "pct_mito", feature2 = "num_ADT_UMIs",
#       group.by = "orig.ident", cols = colors_samples_named
#     ) +
#     theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
#   cor_adt_plot = plot_grid(cor_umis_plot, cor_genes_plot, cor_mito_plot, ncol = 3)
#   ggsave("qc.adt.correlations.png", plot = cor_adt_plot, width = 18, height = 5, units = "in")
#   Sys.sleep(1)
#
# }
#
#
# # plot HTO metrics
# plot_hto_qc = function(seurat_obj) {
#
#   # normalized HTO signal combined with metadata table
#   hto_tbl = GetAssayData(seurat_obj, assay = "HTO") %>% t()
#   hto_tbl = hto_tbl[, sort(colnames(hto_tbl))] %>% as_tibble(rownames = "cell")
#   id_tbl = seurat_obj@meta.data %>% as_tibble(rownames = "cell") %>% select(cell, hash.ID, HTO_classification.global)
#   hto_tbl = full_join(hto_tbl, id_tbl, by = "cell")
#   hto_tbl
#
#   # HTO color scheme
#   colors_hto_names = c(levels(hto_tbl$HTO_classification.global), levels(hto_tbl$hash.ID)) %>% unique()
#   colors_hto = colors_clusters[1:length(colors_hto_names)]
#   names(colors_hto) = colors_hto_names
#
#   # visualize pairs of HTO signals
#   hto_facet_plot =
#     ggplot(sample_frac(hto_tbl), aes(x = .panel_x, y = .panel_y, fill = hash.ID, color = hash.ID)) +
#     geom_point(shape = 16, size = 0.2) +
#     geom_autodensity(color = NA, fill = "gray20") +
#     geom_density2d(color = "black", alpha = 0.5) +
#     scale_color_manual(values = colors_hto) +
#     scale_fill_manual(values = colors_hto) +
#     facet_matrix(vars(-cell, -hash.ID, -HTO_classification.global), layer.diag = 2, layer.upper = 3) +
#     guides(color = guide_legend(override.aes = list(size = 5))) +
#     theme(aspect.ratio = 1, legend.title = element_blank(), strip.background = element_blank())
#   save_plot(filename = "qc.hto.correlation.png", plot = hto_facet_plot, base_height = 8, base_width = 10)
#   Sys.sleep(1)
#   save_plot(filename = "qc.hto.correlation.pdf", plot = hto_facet_plot, base_height = 8, base_width = 10)
#   Sys.sleep(1)
#
#   # number of UMIs for singlets, doublets and negative cells
#   hto_umi_plot =
#     VlnPlot(seurat_obj, features = "num_UMIs", group.by = "HTO_classification.global", pt.size = 0.1) +
#     theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
#     scale_fill_manual(values = colors_hto) +
#     scale_y_continuous(labels = comma)
#   save_plot("qc.hto.umis.png", plot = hto_umi_plot, base_height = 6, base_width = 6)
#   Sys.sleep(1)
#
#   # number of genes for singlets, doublets and negative cells
#   hto_gene_plot =
#     VlnPlot(seurat_obj, features = "num_genes", group.by = "HTO_classification.global", pt.size = 0.1) +
#     theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
#     scale_fill_manual(values = colors_hto) +
#     scale_y_continuous(labels = comma)
#   save_plot("qc.hto.genes.png", plot = hto_gene_plot, base_height = 6, base_width = 6)
#   Sys.sleep(1)
#
#   # number of genes for singlets, doublets and negative cells
#   hto_mito_plot =
#     VlnPlot(seurat_obj, features = "pct_mito", group.by = "HTO_classification.global", pt.size = 0.1) +
#     theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
#     scale_fill_manual(values = colors_hto) +
#     scale_y_continuous(labels = comma)
#   save_plot("qc.hto.mito.png", plot = hto_mito_plot, base_height = 6, base_width = 6)
#   Sys.sleep(1)
#
#   group_var = "HTO_classification.global"
#   Idents(seurat_obj) = group_var
#   plot_umap =
#     DimPlot(
#       seurat_obj, reduction = "umap",
#       cells = sample(colnames(seurat_obj)), pt.size = get_dr_point_size(seurat_obj), cols = colors_hto
#     ) +
#     theme(aspect.ratio = 1)
#   save_plot(glue("dr.umap.{group_var}.png"), plot = plot_umap, base_height = 6, base_width = 8)
#   Sys.sleep(1)
#   save_plot(glue("dr.umap.{group_var}.pdf"), plot = plot_umap, base_height = 6, base_width = 8)
#   Sys.sleep(1)
#
#   group_var = "hash.ID"
#   Idents(seurat_obj) = group_var
#   plot_umap =
#     DimPlot(
#       seurat_obj, reduction = "umap",
#       cells = sample(colnames(seurat_obj)), pt.size = get_dr_point_size(seurat_obj), cols = colors_hto
#     ) +
#     theme(aspect.ratio = 1)
#   save_plot(glue("dr.umap.{group_var}.png"), plot = plot_umap, base_height = 6, base_width = 8)
#   Sys.sleep(1)
#   save_plot(glue("dr.umap.{group_var}.pdf"), plot = plot_umap, base_height = 6, base_width = 8)
#   Sys.sleep(1)
#
#   if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
#
# }
#
# # split Seurat object (by sample by default)
# split_seurat_obj = function(seurat_obj, original_wd, split_var = "orig.ident") {
#
#   # set identity to the column used for splitting
#   s_obj = seurat_obj
#   Idents(s_obj) = split_var
#
#   # clean up metadata
#   s_obj@assays$RNA@var.features = vector()
#   s_obj@assays$RNA@scale.data = matrix()
#   s_obj@reductions = list()
#   s_obj@meta.data = s_obj@meta.data %>% select(-starts_with("snn_res"))
#
#   setwd(original_wd)
#
#   # process each split group
#   split_names = Idents(s_obj) %>% as.character() %>% unique() %>% sort()
#   for (s in split_names) {
#
#     message(glue("split: {s}"))
#     split_obj = subset(s_obj, idents = s)
#     message(glue("split {s} cells: {ncol(split_obj)}"))
#     if (ncol(split_obj) > 10) {
#       split_dir = glue("split-{s}")
#       if (dir.exists(split_dir)) {
#         stop(glue("output dir {split_dir} already exists"))
#       } else {
#         dir.create(split_dir)
#         setwd(split_dir)
#         split_obj = calculate_variance(seurat_obj = split_obj, jackstraw_max_cells = 100)
#         Idents(split_obj) = "orig.ident"
#         saveRDS(split_obj, file = "seurat_obj.rds")
#       }
#     }
#
#     # clean up and return to the main dir before processing the next split
#     if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
#     setwd(original_wd)
#
#   }
#
# }
#
# # merge multiple Seurat objects
# combine_seurat_obj = function(original_wd, sample_analysis_dirs) {
#
#   if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")
#
#   message("\n\n ========== combine samples ========== \n\n")
#
#   seurat_obj_list = list()
#   for (i in 1:length(sample_analysis_dirs)) {
#
#     sample_analysis_dir = sample_analysis_dirs[i]
#     sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
#     sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")
#
#     # check if analysis dir is valid
#     if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
#     # check if seurat object exists
#     if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))
#
#     # load seurat object
#     seurat_obj_list[[i]] = readRDS(sample_seurat_rds)
#
#     # clean up object
#     seurat_obj_list[[i]]@assays$RNA@var.features = vector()
#     seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
#     seurat_obj_list[[i]]@reductions = list()
#     seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))
#
#     # print single sample sample stats
#     sample_name = seurat_obj_list[[i]]$orig.ident %>% as.character() %>% sort()
#     sample_name = sample_name[1]
#     message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
#     write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
#     message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
#     write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
#     message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
#     write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
#     message(" ")
#
#   }
#
#   # merge
#   merged_obj = merge(seurat_obj_list[[1]], seurat_obj_list[2:length(seurat_obj_list)])
#   rm(seurat_obj_list)
#
#   # print combined sample stats
#   message(glue("combined unfiltered cells: {ncol(merged_obj)}"))
#   write(glue("combined unfiltered cells: {ncol(merged_obj)}"), file = "create.log", append = TRUE)
#   message(glue("combined unfiltered genes: {nrow(merged_obj)}"))
#   write(glue("combined unfiltered genes: {nrow(merged_obj)}"), file = "create.log", append = TRUE)
#
#   # filter poorly expressed genes (detected in less than 10 cells)
#   filtered_genes = Matrix::rowSums(GetAssayData(merged_obj, assay = "RNA", slot = "counts") > 0)
#   min_cells = 10
#   if (ncol(merged_obj) > 100000) { min_cells = 50 }
#   filtered_genes = filtered_genes[filtered_genes >= min_cells] %>% names() %>% sort()
#   # keep all HTO and ADT features if present
#   if ("HTO" %in% names(merged_obj@assays)) { filtered_genes = c(filtered_genes, rownames(merged_obj@assays$HTO)) }
#   if ("ADT" %in% names(merged_obj@assays)) { filtered_genes = c(filtered_genes, rownames(merged_obj@assays$ADT)) }
#   merged_obj = subset(merged_obj, features = filtered_genes)
#
#   # encode sample name as factor (also sets alphabetical sample order)
#   merged_obj@meta.data$orig.ident = factor(merged_obj@meta.data$orig.ident)
#
#   # print combined sample stats
#   message(glue("combined cells: {ncol(merged_obj)}"))
#   write(glue("combined cells: {ncol(merged_obj)}"), file = "create.log", append = TRUE)
#   message(glue("combined genes: {nrow(merged_obj)}"))
#   write(glue("combined genes: {nrow(merged_obj)}"), file = "create.log", append = TRUE)
#
#   # print gene/cell minimum cutoffs
#   min_cells = Matrix::rowSums(GetAssayData(merged_obj, assay = "RNA", slot = "counts") > 0) %>% min()
#   min_genes = Matrix::colSums(GetAssayData(merged_obj, assay = "RNA", slot = "counts") > 0) %>% min()
#   message(glue("min cells per gene: {min_cells}"))
#   write(glue("min cells per gene: {min_cells}"), file = "create.log", append = TRUE)
#   message(glue("min genes per cell: {min_genes}"))
#   write(glue("min genes per cell: {min_genes}"), file = "create.log", append = TRUE)
#
#   # check that the full counts table is small enough to fit into an R matrix (max around 100k x 21k)
#   num_matrix_elements = GetAssayData(merged_obj, assay = "RNA", slot = "counts") %>% length()
#   if (num_matrix_elements < 2^31) {
#
#     # save raw counts matrix as a text file
#     counts_raw = GetAssayData(merged_obj, assay = "RNA", slot = "counts") %>% as.matrix()
#     counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
#     # write_csv(counts_raw, path = "counts.raw.csv.gz")
#     fwrite(counts_raw, file = "counts.raw.csv", sep = ",")
#     R.utils::gzip("counts.raw.csv")
#     rm(counts_raw)
#
#     # save normalized counts matrix as a text file
#     counts_norm = GetAssayData(merged_obj, assay = "RNA") %>% as.matrix() %>% round(3)
#     counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
#     # write_csv(counts_norm, path = "counts.normalized.csv.gz")
#     fwrite(counts_norm, file = "counts.normalized.csv", sep = ",")
#     R.utils::gzip("counts.normalized.csv")
#     rm(counts_norm)
#
#   }
#
#   # create a named color scheme to ensure names and colors are in the proper order
#   sample_names = merged_obj$orig.ident %>% as.character() %>% sort() %>% unique()
#   colors_samples_named = colors_samples[1:length(sample_names)]
#   names(colors_samples_named) = sample_names
#
#   vln_theme =
#     theme(
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "none"
#     )
#   suppressMessages({
#     dist_nft_plot =
#       VlnPlot(
#         merged_obj, features = "num_genes", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_nct_plot =
#       VlnPlot(
#         merged_obj, features = "num_UMIs", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_pmt_plot =
#       VlnPlot(
#         merged_obj, features = "pct_mito", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_plot = plot_grid(dist_nft_plot, dist_nct_plot, dist_pmt_plot, ncol = 3)
#     ggsave("qc.distribution.png", plot = dist_plot, width = 20, height = 6, units = "in")
#   })
#   Sys.sleep(1)
#
#   return(merged_obj)
#
# }
#
# # integrate multiple Seurat objects
# integrate_seurat_obj = function(original_wd, sample_analysis_dirs, num_dim) {
#
#   # check if the inputs seems reasonable
#   if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")
#   num_dim = as.integer(num_dim)
#   if (num_dim < 5) { stop("too few dims: ", num_dim) }
#   if (num_dim > 50) { stop("too many dims: ", num_dim) }
#
#   message("\n\n ========== integrate samples ========== \n\n")
#
#   seurat_obj_list = list()
#   var_genes_list = list()
#   exp_genes = c()
#   for (i in 1:length(sample_analysis_dirs)) {
#
#     sample_analysis_dir = sample_analysis_dirs[i]
#     sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
#     sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")
#
#     # check if analysis dir is valid
#     if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
#     # check if seurat object exists
#     if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))
#
#     # load seurat object
#     seurat_obj_list[[i]] = readRDS(sample_seurat_rds)
#     sample_name = seurat_obj_list[[i]]$orig.ident %>% as.character() %>% sort()
#     sample_name = sample_name[1]
#
#     # clean up object
#     seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
#     seurat_obj_list[[i]]@reductions = list()
#     seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))
#
#     # save expressed genes keeping only genes present in all the datasets (for genes to integrate in IntegrateData)
#     if (length(exp_genes) > 0) {
#       exp_genes = intersect(exp_genes, rownames(seurat_obj_list[[i]])) %>% sort()
#     } else {
#       exp_genes = rownames(seurat_obj_list[[i]])
#     }
#
#     # save variable genes
#     var_genes_list[[sample_name]] = VariableFeatures(seurat_obj_list[[i]])
#
#     # print single sample sample stats
#     message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
#     write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
#     message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
#     write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
#     message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
#     write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
#     message(" ")
#
#   }
#
#   # euler plot of variable gene overlaps (becomes unreadable and can take days for many overlaps)
#   if (length(var_genes_list) < 8) {
#     colors_euler = colors_samples[1:length(var_genes_list)]
#     euler_fit = euler(var_genes_list, shape = "ellipse")
#     euler_plot = plot(euler_fit,
#                       fills = list(fill = colors_euler, alpha = 0.7),
#                       edges = list(col = colors_euler))
#     png("variance.vargenes.euler.png", res = 200, width = 5, height = 5, units = "in")
#       print(euler_plot)
#     dev.off()
#   }
#
#   # upset plot of variable gene overlaps
#   png("variance.vargenes.upset.png", res = 200, width = 8, height = 5, units = "in")
#     upset(fromList(var_genes_list), nsets = 50, nintersects = 15, order.by = "freq", mb.ratio = c(0.5, 0.5))
#   dev.off()
#
#   message("\n\n ========== Seurat::FindIntegrationAnchors() ========== \n\n")
#
#   # find the integration anchors
#   anchors = FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = 2000, dims = 1:num_dim)
#   rm(seurat_obj_list)
#
#   message("\n\n ========== Seurat::IntegrateData() ========== \n\n")
#
#   # integrating all genes may cause issues and may not add any relevant information
#   # integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim, features.to.integrate = exp_genes)
#   integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim)
#   rm(anchors)
#
#   # after running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix
#   # the original (uncorrected values) are still stored in the object in the “RNA” assay
#
#   # switch to integrated assay
#   DefaultAssay(integrated_obj) = "integrated"
#
#   # print integrated sample stats
#   message(glue("integrated unfiltered cells: {ncol(integrated_obj)}"))
#   write(glue("integrated unfiltered cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
#   message(glue("integrated unfiltered genes: {nrow(integrated_obj)}"))
#   write(glue("integrated unfiltered genes: {nrow(integrated_obj)}"), file = "create.log", append = TRUE)
#
#   # filter poorly expressed genes (detected in less than 10 cells)
#   filtered_genes = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0)
#   min_cells = 10
#   if (ncol(integrated_obj) > 100000) { min_cells = 50 }
#   filtered_genes = filtered_genes[filtered_genes >= min_cells] %>% names() %>% sort()
#   # keep all HTO and ADT features if present
#   if ("HTO" %in% names(integrated_obj@assays)) { filtered_genes = c(filtered_genes, rownames(integrated_obj@assays$HTO)) }
#   if ("ADT" %in% names(integrated_obj@assays)) { filtered_genes = c(filtered_genes, rownames(integrated_obj@assays$ADT)) }
#   integrated_obj = subset(integrated_obj, features = filtered_genes)
#
#   # encode sample name as factor (also sets alphabetical sample order)
#   integrated_obj@meta.data$orig.ident = factor(integrated_obj@meta.data$orig.ident)
#
#   # print integrated sample stats
#   message(glue("integrated cells: {ncol(integrated_obj)}"))
#   write(glue("integrated cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
#   message(glue("integrated genes: {nrow(GetAssayData(integrated_obj, assay = 'RNA'))}"))
#   write(glue("integrated genes: {nrow(GetAssayData(integrated_obj, assay = 'RNA'))}"), file = "create.log", append = TRUE)
#
#   # print gene/cell minumum cutoffs
#   min_cells = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0) %>% min()
#   min_genes = Matrix::colSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0) %>% min()
#   message(glue("min cells per gene: {min_cells}"))
#   write(glue("min cells per gene: {min_cells}"), file = "create.log", append = TRUE)
#   message(glue("min genes per cell: {min_genes}"))
#   write(glue("min genes per cell: {min_genes}"), file = "create.log", append = TRUE)
#
#   # check that the full counts table is small enough to fit into an R matrix (max around 100k x 21k)
#   num_matrix_elements = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% length()
#   if (num_matrix_elements < 2^31) {
#
#     # save raw counts matrix as a text file
#     counts_raw = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% as.matrix()
#     counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
#     # write_csv(counts_raw, path = "counts.raw.csv.gz")
#     fwrite(counts_raw, file = "counts.raw.csv", sep = ",")
#     R.utils::gzip("counts.raw.csv")
#     rm(counts_raw)
#
#     # save normalized counts matrix as a text file
#     counts_norm = GetAssayData(integrated_obj, assay = "RNA") %>% as.matrix() %>% round(3)
#     counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
#     # write_csv(counts_norm, path = "counts.normalized.csv.gz")
#     fwrite(counts_norm, file = "counts.normalized.csv", sep = ",")
#     R.utils::gzip("counts.normalized.csv")
#     rm(counts_norm)
#
#   }
#
#   vln_theme =
#     theme(
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "none"
#     )
#   suppressMessages({
#     dist_nft_plot =
#       VlnPlot(
#         integrated_obj, features = "num_genes", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_nct_plot =
#       VlnPlot(
#         integrated_obj, features = "num_UMIs", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_pmt_plot =
#       VlnPlot(
#         integrated_obj, features = "pct_mito", group.by = "orig.ident",
#         pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
#       ) +
#       scale_y_continuous(labels = comma) +
#       vln_theme
#     dist_plot = plot_grid(dist_nft_plot, dist_nct_plot, dist_pmt_plot, ncol = 3)
#     ggsave("qc.distribution.png", plot = dist_plot, width = 20, height = 6, units = "in")
#   })
#   Sys.sleep(1)
#
#   return(integrated_obj)
#
# }
#
# # calculate various variance metrics and perform basic analysis
# # PC selection approaches:
# # - PCHeatmap - more supervised, exploring PCs to determine relevant sources of heterogeneity
# # - PCElbowPlot - heuristic that is commonly used and can be calculated instantly
# # - JackStrawPlot - implements a statistical test based on a random null model, but is time-consuming
# # jackStraw procedure is very slow, so skip for large projects (>10,000 cells)
# calculate_variance = function(seurat_obj, jackstraw_max_cells = 10000) {
#
#   s_obj = seurat_obj
#
#   message("\n\n ========== Seurat::FindVariableGenes() ========== \n\n")
#
#   # identify features that are outliers on a 'mean variability plot'
#   # Seurat v3 implements an improved method based on a variance stabilizing transformation ("vst")
#   s_obj = FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#
#   # export highly variable feature information (mean, variance, variance standardized)
#   hvf_tbl = HVFInfo(s_obj) %>% round(3) %>% rownames_to_column("gene") %>% arrange(-variance.standardized)
#   write_excel_csv(hvf_tbl, path = "variance.csv")
#
#   # plot variance
#   var_plot = VariableFeaturePlot(s_obj, pt.size = 0.5)
#   var_plot = LabelPoints(var_plot, points = head(hvf_tbl$gene, 30), repel = TRUE, xnudge = 0, ynudge = 0)
#   ggsave("variance.features.png", plot = var_plot, width = 12, height = 5, units = "in")
#
#   message("\n\n ========== Seurat::ScaleData() ========== \n\n")
#
#   # regress out unwanted sources of variation
#   # regressing uninteresting sources of variation can improve dimensionality reduction and clustering
#   # could include technical noise, batch effects, biological sources of variation (cell cycle stage)
#   # scaled z-scored residuals of these models are stored in scale.data slot
#   # used for dimensionality reduction and clustering
#   # RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData
#   # s_obj = ScaleData(s_obj, features = rownames(s_obj), vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)
#   s_obj = ScaleData(s_obj, vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)
#
#   message("\n\n ========== Seurat::PCA() ========== \n\n")
#
#   # use fewer PCs for small datasets
#   num_pcs = 50
#   if (ncol(s_obj) < 100) num_pcs = 20
#   if (ncol(s_obj) < 25) num_pcs = 5
#
#   # PCA on the scaled data
#   # PCA calculation stored in object[["pca"]]
#   s_obj = RunPCA(s_obj, assay = "RNA", features = VariableFeatures(s_obj), npcs = num_pcs, verbose = FALSE)
#
#   # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
#   pca_plot =
#     DimPlot(
#       s_obj, cells = sample(colnames(s_obj)), group.by = "orig.ident", reduction = "pca",
#       pt.size = 0.5, cols = colors_samples
#     ) +
#     theme(aspect.ratio = 1)
#   ggsave("variance.pca.png", plot = pca_plot, width = 8, height = 6, units = "in")
#
#   message("\n\n ========== Seurat::DimHeatmap() ========== \n\n")
#
#   # PCHeatmap (former) allows for easy exploration of the primary sources of heterogeneity in a dataset
#   if (num_pcs > 15) {
#     png("variance.pca.heatmap.png", res = 300, width = 10, height = 16, units = "in")
#       DimHeatmap(s_obj, reduction = "pca", dims = 1:15, nfeatures = 20, cells = 250, fast = TRUE)
#     dev.off()
#   }
#
#   message("\n\n ========== Seurat::PCElbowPlot() ========== \n\n")
#
#   # a more ad hoc method for determining PCs to use, draw cutoff where there is a clear elbow in the graph
#   elbow_plot = ElbowPlot(s_obj, reduction = "pca", ndims = num_pcs)
#   ggsave("variance.pca.elbow.png", plot = elbow_plot, width = 8, height = 5, units = "in")
#
#   # resampling test inspired by the jackStraw procedure - very slow, so skip for large projects (>10,000 cells)
#   if (ncol(s_obj) < jackstraw_max_cells) {
#
#     message("\n\n ========== Seurat::JackStraw() ========== \n\n")
#
#     # determine statistical significance of PCA scores
#     s_obj = JackStraw(s_obj, assay = "RNA", reduction = "pca", dims = num_pcs, verbose = FALSE)
#
#     # compute Jackstraw scores significance
#     s_obj = ScoreJackStraw(s_obj, reduction = "pca", dims = 1:num_pcs, do.plot = FALSE)
#
#     # plot the results of the JackStraw analysis for PCA significance
#     # significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line)
#     jackstraw_plot =
#       JackStrawPlot(s_obj, reduction = "pca", dims = 1:num_pcs) +
#       guides(col = guide_legend(ncol = 2))
#     ggsave("variance.pca.jackstraw.png", plot = jackstraw_plot, width = 12, height = 6, units = "in")
#
#   }
#
#   return(s_obj)
#
# }
#
# # calculate various variance metrics and perform basic analysis (integrated analysis workflow)
# # specify neighbors for UMAP (default is 30 in Seurat 2 and 3 pre-release)
# calculate_variance_integrated = function(seurat_obj, num_dim, num_neighbors = 30) {
#
#   s_obj = seurat_obj
#
#   num_dim = as.integer(num_dim)
#   if (num_dim < 5) { stop("too few dims: ", num_dim) }
#   if (num_dim > 50) { stop("too many dims: ", num_dim) }
#
#   message("\n\n ========== Seurat::ScaleData() ========== \n\n")
#
#   # s_obj = ScaleData(s_obj, features = rownames(s_obj), verbose = FALSE)
#   s_obj = ScaleData(s_obj, verbose = FALSE)
#
#   message("\n\n ========== Seurat::PCA() ========== \n\n")
#
#   # PCA on the scaled data
#   s_obj = RunPCA(s_obj, npcs = num_dim, verbose = FALSE)
#
#   # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
#   pca_plot =
#     DimPlot(
#       s_obj, reduction = "pca", cells = sample(colnames(s_obj)), group.by = "orig.ident",
#       pt.size = 0.5, cols = colors_samples
#     ) +
#     theme(aspect.ratio = 1)
#   ggsave("variance.pca.png", plot = pca_plot, width = 10, height = 6, units = "in")
#
#   message("\n\n ========== Seurat::RunTSNE() ========== \n\n")
#
#   # use tSNE as a tool to visualize, not for clustering directly on tSNE components
#   # cells within the graph-based clusters determined above should co-localize on the tSNE plot
#   s_obj = RunTSNE(s_obj, reduction = "pca", dims = 1:num_dim, dim.embed = 2)
#
#   # reduce point size for larger datasets
#   dr_pt_size = get_dr_point_size(s_obj)
#
#   # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
#   s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
#   plot_tsne =
#     DimPlot(s_obj, reduction = "tsne", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
#     theme(aspect.ratio = 1)
#   ggsave(glue("dr.tsne.{num_dim}.sample.png"), plot = plot_tsne, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#   ggsave(glue("dr.tsne.{num_dim}.sample.pdf"), plot = plot_tsne, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#
#   message("\n\n ========== Seurat::RunUMAP() ========== \n\n")
#
#   # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
#   s_obj = RunUMAP(s_obj, reduction = "pca", dims = 1:num_dim, n.neighbors = num_neighbors, verbose = FALSE)
#
#   # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
#   s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
#   plot_umap =
#     DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
#     theme(aspect.ratio = 1)
#   ggsave(glue("dr.umap.{num_dim}.sample.png"), plot = plot_umap, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#   ggsave(glue("dr.umap.{num_dim}.sample.pdf"), plot = plot_umap, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#
#   save_metadata(seurat_obj = s_obj)
#
#   return(s_obj)
#
# }
#
# # determine point size for tSNE/UMAP plots (smaller for larger datasets)
# get_dr_point_size = function(seurat_obj) {
#
#   pt_size = 1.8
#   if (ncol(seurat_obj) > 1000) pt_size = 1.2
#   if (ncol(seurat_obj) > 5000) pt_size = 1.0
#   if (ncol(seurat_obj) > 10000) pt_size = 0.8
#   if (ncol(seurat_obj) > 25000) pt_size = 0.6
#
#   return(pt_size)
#
# }
#
# # perform graph-based clustering and tSNE
# # specify neighbors for UMAP and FindNeighbors (default is 30 in Seurat 2 and 3 pre-release)
# calculate_clusters = function(seurat_obj, num_dim, num_neighbors = 30) {
#
#   # check if number of dimensions seems reasonable
#   if (num_dim < 5) { stop("too few dims: ", num_dim) }
#   if (num_dim > 50) { stop("too many dims: ", num_dim) }
#
#   s_obj = seurat_obj
#
#   message("\n\n ========== Seurat::RunTSNE() ========== \n\n")
#
#   # use tSNE as a tool to visualize, not for clustering directly on tSNE components
#   # cells within the graph-based clusters determined above should co-localize on the tSNE plot
#   s_obj = RunTSNE(s_obj, reduction = "pca", dims = 1:num_dim, dim.embed = 2)
#
#   # reduce point size for larger datasets
#   dr_pt_size = get_dr_point_size(s_obj)
#
#   # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
#   s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
#   plot_tsne =
#     DimPlot(s_obj, reduction = "tsne", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
#     theme(aspect.ratio = 1)
#   ggsave(glue("dr.tsne.{num_dim}.sample.png"), plot = plot_tsne, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#   ggsave(glue("dr.tsne.{num_dim}.sample.pdf"), plot = plot_tsne, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#
#   message("\n\n ========== Seurat::RunUMAP() ========== \n\n")
#
#   # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
#   s_obj = RunUMAP(s_obj, reduction = "pca", dims = 1:num_dim, n.neighbors = num_neighbors, verbose = FALSE)
#
#   # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
#   s_obj = set_identity(seurat_obj = s_obj, group_var = "orig.ident")
#   plot_umap =
#     DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size = dr_pt_size, cols = colors_samples) +
#     theme(aspect.ratio = 1)
#   ggsave(glue("dr.umap.{num_dim}.sample.png"), plot = plot_umap, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#   ggsave(glue("dr.umap.{num_dim}.sample.pdf"), plot = plot_umap, width = 10, height = 6, units = "in")
#   Sys.sleep(1)
#
#   message("\n\n ========== Seurat::FindNeighbors() ========== \n\n")
#
#   message("assay: ", DefaultAssay(s_obj))
#   message("num dims: ", num_dim)
#
#   # construct a Shared Nearest Neighbor (SNN) Graph for a given dataset
#   s_obj =
#     FindNeighbors(
#       s_obj, dims = 1:num_dim, k.param = num_neighbors,
#       graph.name = "snn", compute.SNN = TRUE, force.recalc = TRUE
#     )
#
#   message("\n\n ========== Seurat::FindClusters() ========== \n\n")
#
#   message("initial metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))
#
#   # resolutions for graph-based clustering
#   # increased resolution values lead to more clusters (recommendation: 0.6-1.2 for 3K cells, 2-4 for 33K cells)
#   res_range = seq(0.1, 2.0, 0.1)
#   if (ncol(s_obj) > 1000) res_range = c(res_range, 3, 4, 5, 6, 7, 8, 9, 10)
#
#   # algorithm: 1 = original Louvain; 2 = Louvain with multilevel refinement; 3 = SLM
#   # identify clusters of cells by SNN modularity optimization based clustering algorithm
#   s_obj = FindClusters(s_obj, algorithm = 3, resolution = res_range, graph.name = "snn", verbose = FALSE)
#
#   # remove "seurat_clusters" column that is added automatically (added in v3 late dev version)
#   s_obj@meta.data = s_obj@meta.data %>% select(-seurat_clusters)
#
#   message("new metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))
#
#   # create a separate sub-directory for cluster resolution plots
#   clusters_dir = "clusters-resolutions"
#   if (!dir.exists(clusters_dir)) { dir.create(clusters_dir) }
#
#   # for calculated cluster resolutions: remove redundant (same number of clusters), rename, and plot
#   res_cols = str_subset(colnames(s_obj@meta.data), "snn_res")
#   res_cols = sort(res_cols)
#   res_num_clusters_prev = 1
#   for (res in res_cols) {
#
#     # proceed if current resolution has more clusters than previous and less than the color scheme length
#     res_vector = s_obj@meta.data[, res] %>% as.character()
#     res_num_clusters_cur = res_vector %>% n_distinct()
#     if (res_num_clusters_cur > res_num_clusters_prev && res_num_clusters_cur < length(colors_clusters)) {
#
#       # check if the resolution still has original labels (characters starting with 0)
#       if (min(res_vector) == "0") {
#
#         # convert to character vector
#         s_obj@meta.data[, res] = as.character(s_obj@meta.data[, res])
#         # relabel identities so they start with 1 and not 0
#         s_obj@meta.data[, res] = as.numeric(s_obj@meta.data[, res]) + 1
#         # pad with 0s to avoid sorting issues
#         s_obj@meta.data[, res] = str_pad(s_obj@meta.data[, res], width = 2, side = "left", pad = "0")
#         # pad with "C" to avoid downstream numeric conversions
#         s_obj@meta.data[, res] = str_c("C", s_obj@meta.data[, res])
#         # encode as a factor
#         s_obj@meta.data[, res] = factor(s_obj@meta.data[, res])
#
#       }
#
#       # resolution value based on resolution column name
#       res_val = sub("snn_res\\.", "", res)
#
#       # plot file name
#       res_str = gsub("\\.", "", res)
#       dr_filename = glue("{clusters_dir}/dr.{DefaultAssay(s_obj)}.{num_dim}.{res_str}.clust{res_num_clusters_cur}")
#
#       s_obj = plot_clusters(seurat_obj = s_obj, resolution = res_val, filename_base = dr_filename)
#
#       # add blank line to make output easier to read
#       message(" ")
#
#     } else {
#
#       # remove resolution if the number of clusters is same as previous
#       s_obj@meta.data = s_obj@meta.data %>% select(-one_of(res))
#
#     }
#
#     # update resolution cluster count for next iteration
#     res_num_clusters_prev = res_num_clusters_cur
#
#   }
#
#   message("updated metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))
#
#   save_metadata(seurat_obj = s_obj)
#
#   return(s_obj)
#
# }
#
# # compile all cell metadata into a single table
# save_metadata = function(seurat_obj) {
#
#   s_obj = seurat_obj
#   metadata_tbl = s_obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>% mutate(sample_name = orig.ident)
#   tsne_tbl = s_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
#   umap_tbl = s_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
#   cells_metadata = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
#   cells_metadata = cells_metadata %>% arrange(cell)
#   write_excel_csv(cells_metadata, path = "metadata.csv")
#
# }
#
# # plot tSNE with color-coded clusters at specified resolution
# plot_clusters = function(seurat_obj, resolution, filename_base) {
#
#   s_obj = seurat_obj
#
#   # set identities based on specified resolution
#   s_obj = set_identity(seurat_obj = s_obj, group_var = resolution)
#
#   # print stats
#   num_clusters = Idents(s_obj) %>% as.character() %>% n_distinct()
#   message("resolution: ", resolution)
#   message("num clusters: ", num_clusters)
#
#   # generate plot if there is a reasonable number of clusters
#   if (num_clusters > 1 && num_clusters < length(colors_clusters)) {
#
#     # shuffle cells so they appear randomly and one group does not show up on top
#     plot_tsne =
#       DimPlot(
#         s_obj, reduction = "tsne", cells = sample(colnames(s_obj)),
#         pt.size = get_dr_point_size(s_obj), cols = colors_clusters
#       ) +
#       theme(aspect.ratio = 1)
#     ggsave(glue("{filename_base}.tsne.png"), plot = plot_tsne, width = 9, height = 6, units = "in")
#     Sys.sleep(1)
#     ggsave(glue("{filename_base}.tsne.pdf"), plot = plot_tsne, width = 9, height = 6, units = "in")
#     Sys.sleep(1)
#
#     plot_umap =
#       DimPlot(
#         s_obj, reduction = "umap", cells = sample(colnames(s_obj)),
#         pt.size = get_dr_point_size(s_obj), cols = colors_clusters
#       ) +
#       theme(aspect.ratio = 1)
#     ggsave(glue("{filename_base}.umap.png"), plot = plot_umap, width = 9, height = 6, units = "in")
#     Sys.sleep(1)
#     ggsave(glue("{filename_base}.umap.pdf"), plot = plot_umap, width = 9, height = 6, units = "in")
#     Sys.sleep(1)
#
#     if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
#
#   }
#
#   return(s_obj)
#
# }
#
# # check grouping variable/resolution against existing meta data columns
# check_group_var = function(seurat_obj, group_var) {
#
#   s_obj = seurat_obj
#
#   # check if the grouping variable is one of meta data columns
#   if (!(group_var %in% colnames(s_obj@meta.data))) {
#
#     # check if grouping variable is the resolution value (X.X instead of res.X.X)
#     res_column = str_c("snn_res.", group_var)
#     if (res_column %in% colnames(s_obj@meta.data)) {
#       group_var = res_column
#     } else {
#       stop("unknown grouping variable: ", group_var)
#     }
#
#   }
#
#   return(group_var)
#
# }
#
# # set identity based on a specified variable/resolution
# set_identity = function(seurat_obj, group_var) {
#
#   s_obj = seurat_obj
#
#   group_var = check_group_var(seurat_obj = s_obj, group_var = group_var)
#
#   # set identities based on selected grouping variable
#   message("setting grouping variable: ", group_var)
#   Idents(s_obj) = group_var
#
#   return(s_obj)
#
# }
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
#
# # gather metadata and calculate cluster stats (number of cells)
# calculate_cluster_stats = function(seurat_obj, label) {
#
#   message("\n\n ========== calculate cluster stats ========== \n\n")
#
#   message("cluster names: ", str_c(levels(seurat_obj), collapse = ", "))
#
#   # compile relevant cell metadata into a single table
#   seurat_obj$cluster = Idents(seurat_obj)
#   metadata_tbl = seurat_obj@meta.data %>% rownames_to_column("cell") %>% as_tibble() %>%
#     select(cell, num_UMIs, num_genes, pct_mito, sample_name = orig.ident, cluster)
#   tsne_tbl = seurat_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
#   umap_tbl = seurat_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
#   cells_metadata = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
#   cells_metadata = cells_metadata %>% arrange(cell)
#   write_excel_csv(cells_metadata, path = glue("metadata.{label}.csv"))
#   Sys.sleep(1)
#
#   # get number of cells split by cluster and by sample (orig.ident)
#   summary_cluster_sample =
#     cells_metadata %>%
#     select(cluster, sample_name) %>%
#     mutate(num_cells_total = n()) %>%
#     group_by(sample_name) %>%
#     mutate(num_cells_sample = n()) %>%
#     group_by(cluster) %>%
#     mutate(num_cells_cluster = n()) %>%
#     group_by(cluster, sample_name) %>%
#     mutate(num_cells_cluster_sample = n()) %>%
#     ungroup() %>%
#     distinct() %>%
#     mutate(
#       pct_cells_cluster = num_cells_cluster / num_cells_total,
#       pct_cells_cluster_sample = num_cells_cluster_sample / num_cells_sample
#     ) %>%
#     mutate(
#       pct_cells_cluster = round(pct_cells_cluster * 100, 1),
#       pct_cells_cluster_sample = round(pct_cells_cluster_sample * 100, 1)
#     ) %>%
#     arrange(cluster, sample_name)
#
#   # get number of cells split by cluster (ignore samples)
#   summary_cluster = summary_cluster_sample %>% select(-contains("sample")) %>% distinct()
#   write_excel_csv(summary_cluster, path = glue("summary.{label}.csv"))
#   Sys.sleep(1)
#
#   # export results split by sample if multiple samples are present
#   num_samples = cells_metadata %>% pull(sample_name) %>% n_distinct()
#   if (num_samples > 1) {
#     write_excel_csv(summary_cluster_sample, path = glue("summary.{label}.per-sample.csv"))
#     Sys.sleep(1)
#   }
#
# }
#
# # calculate cluster average expression (non-log space)
# calculate_cluster_expression = function(seurat_obj, label) {
#
#   message("\n\n ========== calculate cluster average expression ========== \n\n")
#
#   message("cluster names: ", str_c(levels(seurat_obj), collapse = ", "))
#
#   # gene expression for an "average" cell in each identity class (averaging and output are in non-log space)
#   cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", verbose = FALSE)[["RNA"]]
#   cluster_avg_exp = cluster_avg_exp %>% round(3) %>% rownames_to_column("gene") %>% arrange(gene)
#   write_excel_csv(cluster_avg_exp, path = glue("expression.mean.{label}.csv"))
#   Sys.sleep(1)
#
#   # cluster averages split by sample (orig.ident)
#   num_samples = n_distinct(seurat_obj@meta.data$orig.ident)
#   if (num_samples > 1) {
#     sample_avg_exp = AverageExpression(seurat_obj, assay = "RNA", add.ident = "orig.ident", verbose = FALSE)[["RNA"]]
#     sample_avg_exp = sample_avg_exp %>% round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
#     write_excel_csv(sample_avg_exp, path = glue("expression.mean.{label}.per-sample.csv"))
#     Sys.sleep(1)
#   }
#
# }
#
# # calculate cluster markers (compared to all other cells) and plot top ones
# # tests:
# # - roc: ROC test returns the classification power (ranging from 0 - random, to 1 - perfect)
# # - wilcox: Wilcoxon rank sum test (default in Seurat 2)
# # - bimod: Likelihood-ratio test for single cell gene expression (McDavid, Bioinformatics, 2013) (default in Seurat 1)
# # - tobit: Tobit-test for differential gene expression (Trapnell, Nature Biotech, 2014)
# # - MAST: GLM-framework that treates cellular detection rate as a covariate (Finak, Genome Biology, 2015)
# # pairwise option compares each cluster to each of the other clusters to yield markers that are both local and global
# calculate_cluster_markers = function(seurat_obj, label, test, pairwise = FALSE) {
#
#   message("\n\n ========== calculate cluster markers ========== \n\n")
#
#   message("cluster set: ", label)
#   message("marker test: ", test)
#
#   # get cluster names
#   clusters = Idents(seurat_obj) %>% as.character() %>% unique() %>% sort()
#
#   # use only clusters with more than 10 cells
#   clusters = clusters[table(Idents(seurat_obj)) > 10]
#
#   if (!pairwise) {
#
#     # standard cluster markers calculation
#
#     markers_dir = "markers-global"
#
#     # capture output to avoid excessive warnings
#     markers_log =
#       capture.output({
#         all_markers =
#           FindAllMarkers(
#             seurat_obj, assay = "RNA", test.use = test, logfc.threshold = log(1.2), min.pct = 0.2,
#             only.pos = FALSE, min.diff.pct = -Inf, verbose = FALSE
#           )
#       }, type = "message")
#
#     # do some light filtering and clean up (ROC test returns slightly different output)
#     if (test == "roc") {
#
#       all_markers =
#         all_markers %>%
#         select(cluster, gene, logFC = avg_diff, myAUC, power) %>%
#         filter(power > 0.3) %>%
#         mutate(logFC = round(logFC, 5), myAUC = round(myAUC, 5), power = round(power, 5)) %>%
#         arrange(cluster, -power)
#       top_markers = all_markers %>% filter(logFC > 0)
#       top_markers = top_markers %>% group_by(cluster) %>% top_n(50, power) %>% ungroup()
#
#     } else {
#
#       all_markers =
#         all_markers %>%
#         select(cluster, gene, logFC = avg_logFC, p_val, p_val_adj) %>%
#         filter(p_val_adj < 0.001) %>%
#         mutate(logFC = round(logFC, 5)) %>%
#         arrange(cluster, p_val_adj, p_val)
#       top_markers = all_markers %>% filter(logFC > 0)
#       top_markers = top_markers %>% group_by(cluster) %>% top_n(50, logFC) %>% ungroup()
#
#     }
#
#   } else {
#
#     # pairwise (each cluster versus each other cluster) cluster markers calculation
#
#     markers_dir = "markers-pairwise"
#
#     # initialize empty results tibble
#     unfiltered_markers = tibble(
#       cluster = character(),
#       cluster2 = character(),
#       gene = character(),
#       logFC = numeric(),
#       p_val = numeric(),
#       p_val_adj = numeric()
#     )
#
#     # check each cluster combination
#     for (cluster1 in clusters) {
#       for (cluster2 in setdiff(clusters, cluster1)) {
#
#         # find differentially expressed genes between two specific clusters
#         # low fold change cutoff to maximize chance of appearing in all comparisons
#         # capture output to avoid excessive warnings
#         markers_log =
#           capture.output({
#             cur_markers =
#               FindMarkers(
#                 seurat_obj, assay = "RNA", ident.1 = cluster1, ident.2 = cluster2, test.use = test,
#                 logfc.threshold = log(1.1), min.pct = 0.1,
#                 only.pos = TRUE, min.diff.pct = -Inf, verbose = FALSE
#               )
#           }, type = "message")
#
#         # clean up markers table (would need to be modified for "roc" test)
#         cur_markers =
#           cur_markers %>%
#           rownames_to_column("gene") %>%
#           mutate(cluster = cluster1) %>%
#           mutate(cluster2 = cluster2) %>%
#           filter(p_val_adj < 0.01) %>%
#           mutate(logFC = round(avg_logFC, 5)) %>%
#           select(one_of(colnames(unfiltered_markers)))
#
#         # add current cluster combination genes to the table of all markers
#         unfiltered_markers = bind_rows(unfiltered_markers, cur_markers)
#
#       }
#     }
#
#     # adjust test name for output
#     test = glue("pairwise.{test}")
#
#     # sort the markers to make the table more readable
#     unfiltered_markers =
#       unfiltered_markers %>%
#       distinct() %>%
#       add_count(cluster, gene) %>%
#       rename(cluster_gene_n = n) %>%
#       arrange(cluster, gene, cluster2)
#
#     # filter for genes that are significant compared to all other clusters
#     all_markers =
#       unfiltered_markers %>%
#       filter(cluster_gene_n == (length(clusters) - 1)) %>%
#       select(-cluster_gene_n)
#
#     # extract the lowest and highest fold changes and p-values
#     all_markers =
#       all_markers %>%
#       group_by(cluster, gene) %>%
#       summarize_at(
#         c("logFC", "p_val", "p_val_adj"),
#         list(min = min, max = max)
#       ) %>%
#       ungroup() %>%
#       arrange(cluster, -logFC_min)
#
#     top_markers = all_markers %>% group_by(cluster) %>% top_n(50, logFC_min) %>% ungroup()
#
#   }
#
#   # create a separate sub-directory for all markers
#   if (!dir.exists(markers_dir)) { dir.create(markers_dir) }
#
#   # filename prefix
#   filename_base = glue("{markers_dir}/markers.{label}.{test}")
#
#   # save unfiltered markers for pairwise comparisons
#   if (pairwise) {
#     unfiltered_markers_csv = glue("{filename_base}.unfiltered.csv")
#     message("unfiltered markers: ", unfiltered_markers_csv)
#     write_excel_csv(unfiltered_markers, path = unfiltered_markers_csv)
#     Sys.sleep(1)
#   }
#
#   all_markers_csv = glue("{filename_base}.all.csv")
#   message("all markers: ", all_markers_csv)
#   write_excel_csv(all_markers, path = all_markers_csv)
#   Sys.sleep(1)
#
#   top_markers_csv = glue("{filename_base}.top.csv")
#   message("top markers: ", top_markers_csv)
#   write_excel_csv(top_markers, path = top_markers_csv)
#   Sys.sleep(1)
#
#   # plot cluster markers heatmap
#   plot_cluster_markers(seurat_obj, markers_tbl = all_markers, num_genes = c(5, 10, 20), filename_base = filename_base)
#
#   # plot top cluster markers for each cluster
#   for (cluster_name in clusters) {
#     filename_cluster_base = glue("{markers_dir}/markers.{label}-{cluster_name}.{test}")
#     cluster_markers = top_markers %>% filter(cluster == cluster_name)
#     if (nrow(cluster_markers) > 9) {
#       Sys.sleep(1)
#       top_cluster_markers = cluster_markers %>% head(12) %>% pull(gene)
#       plot_genes(seurat_obj, genes = top_cluster_markers, filename_base = filename_cluster_base)
#     }
#
#   }
#
# }
#
# # generate cluster markers heatmap
# plot_cluster_markers = function(seurat_obj, markers_tbl, num_genes, filename_base) {
#
#   # adjust pairwise clusters to match the standard format
#   if ("logFC_min" %in% colnames(markers_tbl)) {
#     markers_tbl = markers_tbl %>% mutate(logFC = logFC_min)
#   }
#
#   # keep only the top cluster for each gene so each gene appears once
#   markers_tbl = markers_tbl %>% filter(logFC > 0)
#   markers_tbl = markers_tbl %>% group_by(gene) %>% top_n(1, logFC) %>% slice(1) %>% ungroup()
#
#   num_clusters = Idents(seurat_obj) %>% as.character() %>% n_distinct()
#   marker_genes = markers_tbl %>% pull(gene) %>% unique() %>% sort()
#   cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", features = marker_genes, verbose = FALSE)[["RNA"]]
#   cluster_avg_exp = cluster_avg_exp %>% as.matrix() %>% log1p()
#   cluster_avg_exp = cluster_avg_exp[rowSums(cluster_avg_exp) > 0, ]
#
#   # heatmap settings
#   hm_colors = colorRampPalette(c("#053061", "#FFFFFF", "#E41A1C"))(51)
#   hm_width = ( num_clusters / 2 ) + 2
#
#   for (ng in num_genes) {
#
#     hm_base = glue("{filename_base}.heatmap.top{ng}")
#
#     markers_top_tbl = markers_tbl %>% group_by(cluster) %>% top_n(ng, logFC) %>% ungroup()
#     markers_top_tbl = markers_top_tbl %>% arrange(cluster, -logFC)
#
#     # generate the scaled expression matrix and save the text version
#     hm_mat = cluster_avg_exp[markers_top_tbl$gene, ]
#     hm_mat = hm_mat %>% t() %>% scale() %>% t()
#     hm_mat %>% round(3) %>% as_tibble(rownames = "gene") %>% write_excel_csv(path = glue("{hm_base}.csv"))
#     Sys.sleep(1)
#
#     # set outliers to 95th percentile to yield a more balanced color scale
#     scale_cutoff = as.numeric(quantile(abs(hm_mat), 0.95))
#     hm_mat[hm_mat > scale_cutoff] = scale_cutoff
#     hm_mat[hm_mat < -scale_cutoff] = -scale_cutoff
#
#     # generate the heatmap
#     ph_obj = pheatmap(
#       hm_mat, scale = "none", color = hm_colors, border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE,
#       fontsize = 10, fontsize_row = 8, fontsize_col = 12, show_colnames = TRUE,
#       main = glue("Cluster Markers: Top {ng}")
#     )
#
#     png(glue("{hm_base}.png"), width = hm_width, height = 10, units = "in", res = 300)
#       grid::grid.newpage()
#       grid::grid.draw(ph_obj$gtable)
#     dev.off()
#     Sys.sleep(1)
#     pdf(glue("{hm_base}.pdf"), width = hm_width, height = 10)
#       grid::grid.newpage()
#       grid::grid.draw(ph_obj$gtable)
#     dev.off()
#     Sys.sleep(1)
#
#   }
#
# }
#
# # calculate differentially expressed genes within each cluster
# calculate_cluster_de_genes = function(seurat_obj, label, test, group_var = "orig.ident") {
#
#   message("\n\n ========== calculate cluster DE genes ========== \n\n")
#
#   # create a separate sub-directory for differential expression results
#   de_dir = glue("diff-expression-{group_var}")
#   if (!dir.exists(de_dir)) { dir.create(de_dir) }
#
#   # common settings
#   num_de_genes = 50
#
#   # results table
#   de_all_genes_tbl = tibble()
#
#   # get DE genes for each cluster
#   clusters = levels(seurat_obj)
#   for (clust_name in clusters) {
#
#     message(glue("calculating DE genes for cluster {clust_name}"))
#
#     # subset to the specific cluster
#     clust_obj = subset(seurat_obj, idents = clust_name)
#
#     # revert back to the grouping variable sample/library labels
#     Idents(clust_obj) = group_var
#
#     message("cluster cells: ", ncol(clust_obj))
#     message("cluster groups: ", paste(levels(clust_obj), collapse = ", "))
#
#     # continue if cluster has multiple groups and more than 10 cells in each group
#     if (n_distinct(Idents(clust_obj)) > 1 && min(table(Idents(clust_obj))) > 10) {
#
#       # scale data for heatmap
#       clust_obj = ScaleData(clust_obj, assay = "RNA", vars.to.regress = c("num_UMIs", "pct_mito"))
#
#       # iterate through sample/library combinations (relevant if more than two)
#       group_combinations = combn(levels(clust_obj), m = 2, simplify = TRUE)
#       for (combination_num in 1:ncol(group_combinations)) {
#
#         # determine combination
#         g1 = group_combinations[1, combination_num]
#         g2 = group_combinations[2, combination_num]
#         comparison_label = glue("{g1}-vs-{g2}")
#         message(glue("comparison: {clust_name} {g1} vs {g2}"))
#
#         filename_label = glue("{de_dir}/de.{label}-{clust_name}.{comparison_label}.{test}")
#
#         # find differentially expressed genes (default Wilcoxon rank sum test)
#         de_genes = FindMarkers(clust_obj, ident.1 = g1, ident.2 = g2, assay = "RNA",
#                                test.use = test, logfc.threshold = log(1), min.pct = 0.1, only.pos = FALSE,
#                                print.bar = FALSE)
#
#         # perform some light filtering and clean up
#         de_genes =
#           de_genes %>%
#           rownames_to_column("gene") %>%
#           mutate(cluster = clust_name, group1 = g1, group2 = g2, de_test = test) %>%
#           select(cluster, group1, group2, de_test, gene, logFC = avg_logFC, p_val, p_val_adj) %>%
#           mutate(
#             logFC = round(logFC, 3),
#             p_val = if_else(p_val < 0.00001, p_val, round(p_val, 5)),
#             p_val_adj = if_else(p_val_adj < 0.00001, p_val_adj, round(p_val_adj, 5))
#           ) %>%
#           arrange(p_val_adj, p_val)
#
#         message(glue("{comparison_label} num genes: {nrow(de_genes)}"))
#
#         # save stats table
#         write_excel_csv(de_genes, path = glue("{filename_label}.csv"))
#
#         # add cluster genes to all genes
#         de_all_genes_tbl = bind_rows(de_all_genes_tbl, de_genes)
#
#         # heatmap of top genes
#         if (nrow(de_genes) > 5) {
#           top_de_genes = de_genes %>% top_n(num_de_genes, -p_val_adj) %>% arrange(logFC) %>% pull(gene)
#           plot_hm = DoHeatmap(clust_obj, features = top_de_genes, assay = "RNA", slot = "scale.data")
#           heatmap_prefix = glue("{filename_label}.heatmap.top{num_de_genes}")
#           ggsave(glue("{heatmap_prefix}.png"), plot = plot_hm, width = 15, height = 10, units = "in")
#           Sys.sleep(1)
#           ggsave(glue("{heatmap_prefix}.pdf"), plot = plot_hm, width = 15, height = 10, units = "in")
#           Sys.sleep(1)
#         }
#
#       }
#
#     } else {
#
#       message("skip cluster: ", clust_name)
#
#     }
#
#     message(" ")
#
#   }
#
#   # save stats table
#   write_excel_csv(de_all_genes_tbl, path = glue("{de_dir}/de.{label}.{group_var}.{test}.all.csv"))
#   de_all_genes_tbl = de_all_genes_tbl %>% filter(p_val_adj < 0.01)
#   write_excel_csv(de_all_genes_tbl, path = glue("{de_dir}/de.{label}.{group_var}.{test}.sig.csv"))
#
# }


