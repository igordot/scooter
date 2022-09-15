#!/usr/bin/env Rscript

"
Analysis of 10x Genomics Chromium single cell RNA-seq data using Seurat (version >=4.0.2) starting with Cell Ranger output.

Basic workflow steps:
  1 - create - import counts matrix, perform initial QC, and calculate various variance metrics (slowest step)
  2 - cluster - perform clustering based on number of PCs
  3 - identify - identify clusters based on specified clustering/resolution (higher resolution for more clusters)

Additional optional steps:
  combine - merge multiple samples/libraries (no batch correction)
  integrate - perform integration (batch correction) across multiple samples or sample batches
  de - differential expression between samples/libraries within clusters
  adt - add antibody-derived tag data [retired]
  hto - add hashtag oligo data and split by hashtag
  plot umap - generate gene expression overlaid on a UMAP based on a table of genes and groups

Usage:
  scrna-10x-scooter.R create <analysis_dir> <sample_name> <sample_dir> [--min_genes=<n> --max_genes=<n> --mt=<n>]
  scrna-10x-scooter.R cluster <analysis_dir> <num_dim>
  scrna-10x-scooter.R identify <analysis_dir> <resolution>
  scrna-10x-scooter.R combine <analysis_dir> <sample_analysis_dir>...
  scrna-10x-scooter.R integrate <analysis_dir> <reduction> <num_dim> <batch_analysis_dir>...
  scrna-10x-scooter.R de <analysis_dir> <resolution> <group_var>
  scrna-10x-scooter.R hto <analysis_dir>
  scrna-10x-scooter.R plot umap <analysis_dir> <genes_csv>
  scrna-10x-scooter.R --help

Options:
  --min_genes=<n>   cutoff for minimum number of genes per cell (2nd percentile if not specified)
  --max_genes=<n>   cutoff for maximum number of genes per cell (98th percentile if not specified)
  --mt=<n>          cutoff for mitochondrial genes percentage per cell [default: 10]
  -h, --help        show this screen
" -> doc


# ========== functions ==========


# load dependencies
load_libraries = function() {

  message("\n\n ========== load libraries ========== \n\n")

  suppressPackageStartupMessages({
    library(magrittr)
    library(glue)
    library(Seurat)
    library(future)
    library(Matrix)
    library(tidyverse)
    library(data.table)
    library(cowplot)
    library(scales)
    library(pheatmap)
    library(RColorBrewer)
    library(ggsci)
    library(eulerr)
    library(UpSetR)
    library(ggforce)
    # remotes::install_github("igordot/scooter")
    library(scooter)
  })

  theme_set(theme_cowplot())

}

# convert a sparse matrix of counts to a Seurat object and generate some QC plots
create_seurat_obj_qc = function(seurat_obj) {

  s_obj = seurat_obj

  message("\n\n ========== nFeature_RNA/nCount_RNA/percent_mito plots ========== \n\n")

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = s_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  suppressMessages({
    dist_unfilt_nft_plot =
      VlnPlot(
        s_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_unfilt_nct_plot =
      VlnPlot(
        s_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_unfilt_pmt_plot =
      VlnPlot(
        s_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_unfilt_plot = plot_grid(dist_unfilt_nft_plot, dist_unfilt_nct_plot, dist_unfilt_pmt_plot, ncol = 3)
    ggsave("qc.distribution.unfiltered.png", plot = dist_unfilt_plot, width = 10, height = 6, units = "in")
  })
  Sys.sleep(1)

  cor_ncr_nfr_plot =
    FeatureScatter(
      s_obj, feature1 = "num_UMIs", feature2 = "num_genes", group.by = "orig.ident", cols = colors_samples
    ) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  cor_ncr_pmt_plot =
    FeatureScatter(
      s_obj, feature1 = "num_UMIs", feature2 = "pct_mito", group.by = "orig.ident", cols = colors_samples
    ) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  cor_nfr_pmt_plot =
    FeatureScatter(
      s_obj, feature1 = "num_genes", feature2 = "pct_mito", group.by = "orig.ident", cols = colors_samples
    ) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  cor_unfilt_plot = plot_grid(cor_ncr_nfr_plot, cor_ncr_pmt_plot, cor_nfr_pmt_plot, ncol = 3)
  ggsave("qc.correlations.unfiltered.png", plot = cor_unfilt_plot, width = 18, height = 5, units = "in")
  Sys.sleep(1)

  # check distribution of gene counts and mitochondrial percentage
  low_quantiles = c(0.05, 0.02, 0.01, 0.001)
  high_quantiles = c(0.95, 0.98, 0.99, 0.999)
  message("num genes low percentiles:")
  s_obj$num_genes %>% quantile(low_quantiles) %>% round(1) %>% print()
  message(" ")
  message("num genes high percentiles:")
  s_obj$num_genes %>% quantile(high_quantiles) %>% round(1) %>% print()
  message(" ")
  message("pct mito high percentiles:")
  s_obj$pct_mito %>% quantile(high_quantiles) %>% round(1) %>% print()
  message(" ")

  # save unfiltered cell metadata
  s_obj@meta.data %>%
    rownames_to_column("cell") %>% as_tibble() %>%
    mutate(sample_name = orig.ident) %>%
    write_csv("metadata.unfiltered.csv")

  return(s_obj)

}

# filter data by number of genes and mitochondrial percentage
filter_data = function(seurat_obj, min_umis = 1000, min_genes = NULL, max_genes = NULL, max_mt = 10) {

  s_obj = seurat_obj

  message("\n\n ========== filter data matrix ========== \n\n")

  # log the unfiltered gene numbers to file
  write(glue("unfiltered min num genes: {min(s_obj$num_genes)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered max num genes: {max(s_obj$num_genes)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered mean num genes: {round(mean(s_obj$num_genes), 3)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered median num genes: {median(s_obj$num_genes)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered min num UMIs: {min(s_obj$num_UMIs)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered max num UMIs: {max(s_obj$num_UMIs)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered mean num UMIs: {round(mean(s_obj$num_UMIs), 3)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered median num UMIs: {median(s_obj$num_UMIs)}"), file = "create.log", append = TRUE)

  # convert arguments to integers (command line arguments end up as characters)
  min_genes = as.numeric(min_genes)
  max_genes = as.numeric(max_genes)
  max_mt = as.numeric(max_mt)

  # default cutoffs (gene numbers rounded to nearest 10)
  # as.numeric() converts NULLs to 0 length numerics, so can't use is.null()
  if (!length(min_genes)) min_genes = s_obj$num_genes %>% quantile(0.02, names = FALSE) %>% round(-1)
  if (!length(max_genes)) max_genes = s_obj$num_genes %>% quantile(0.98, names = FALSE) %>% round(-1)
  if (!length(max_mt)) max_mt = 10

  message(glue("min genes cutoff: {min_genes}"))
  message(glue("max genes cutoff: {max_genes}"))
  message(glue("max mitochondrial percentage cutoff: {max_mt}"))
  message(" ")

  # log the cutoffs to file
  write(glue("min genes cutoff: {min_genes}"), file = "create.log", append = TRUE)
  write(glue("max genes cutoff: {max_genes}"), file = "create.log", append = TRUE)
  write(glue("max mitochondrial percentage cutoff: {max_mt}"), file = "create.log", append = TRUE)

  message(glue("imported cells: {ncol(s_obj)}"))
  message(glue("imported genes: {nrow(s_obj)}"))

  # set a minimum UMIs cutoff
  # min_umis = 0
  # if (nrow(s_obj) > 1000) min_umis = 1000

  # filter
  cells_subset =
    seurat_obj@meta.data %>%
    rownames_to_column("cell") %>%
    filter(nCount_RNA > min_umis & nFeature_RNA > min_genes & nFeature_RNA < max_genes & pct_mito < max_mt) %>%
    pull(cell)
  s_obj = subset(s_obj, cells = cells_subset)

  message("filtered cells: ", ncol(s_obj))
  message("filtered genes: ", nrow(s_obj))

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = s_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  suppressMessages({
    dist_filt_nft_plot =
      VlnPlot(
        s_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_filt_nct_plot =
      VlnPlot(
        s_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_filt_pmt_plot =
      VlnPlot(
        s_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_filt_plot = plot_grid(dist_filt_nft_plot, dist_filt_nct_plot, dist_filt_pmt_plot, ncol = 3)
    ggsave("qc.distribution.filtered.png", plot = dist_filt_plot, width = 10, height = 6, units = "in")
  })
  Sys.sleep(1)

  # after removing unwanted cells from the dataset, normalize the data
  # LogNormalize:
  # - normalizes the gene expression measurements for each cell by the total expression
  # - multiplies this by a scale factor (10,000 by default)
  # - log-transforms the result
  s_obj = NormalizeData(s_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

  # save counts matrix as a basic gzipped text file
  # object@data stores normalized and log-transformed single cell expression
  # used for visualizations, such as violin and feature plots, most diff exp tests, finding high-variance genes
  counts_norm = GetAssayData(s_obj) %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  fwrite(counts_norm, file = "counts.normalized.csv.gz", sep = ",", nThread = 4)

  # log to file
  write(glue("filtered cells: {ncol(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("filtered genes: {nrow(s_obj)}"), file = "create.log", append = TRUE)
  write(glue("filtered mean num genes: {round(mean(s_obj$num_genes), 3)}"), file = "create.log", append = TRUE)
  write(glue("filtered median num genes: {median(s_obj$num_genes)}"), file = "create.log", append = TRUE)
  write(glue("filtered mean num UMIs: {round(mean(s_obj$num_UMIs), 3)}"), file = "create.log", append = TRUE)
  write(glue("filtered median num UMIs: {median(s_obj$num_UMIs)}"), file = "create.log", append = TRUE)

  return(s_obj)

}

# add antibody-derived tags (ADT) data to a Seurat object
# add_adt_assay = function(seurat_obj, sample_name, sample_dir) {
add_adt_assay_qc = function(seurat_obj, sample_name) {

#  message("\n\n ========== import antibody-derived tags (ADT) ========== \n\n")
#
#  message("loading ADT matrix for sample: ", sample_name)
#
#  # check if sample dir is valid
#  if (!dir.exists(sample_dir)) { stop(glue("dir {sample_dir} does not exist")) }
#  if (!file.exists(glue("{sample_dir}/matrix.mtx.gz"))) { stop(glue("dir does not contain matrix.mtx.gz")) }
#
#  adt_mat = Read10X(data.dir = sample_dir, gene.column = 1)
#  adt_mat = as.matrix(adt_mat)
#
#  # removed "unmapped" ADT
#  if (rownames(adt_mat)[length(rownames(adt_mat))] == "unmapped") {
#    adt_mat = adt_mat[1:length(rownames(adt_mat))-1, ]
#  }

  adt_mat = GetAssayData(seurat_obj, assay = "ADT", slot = "counts")

  message(glue("ADT library {sample_name} cells: {ncol(adt_mat)}"))
  message(glue("ADT library {sample_name} ADTs: {nrow(adt_mat)}"))

  # log to file
  write(glue("ADT library {sample_name} cells: {ncol(adt_mat)}"), file = "create.log", append = TRUE)
  write(glue("ADT library {sample_name} ADTs: {nrow(adt_mat)}"), file = "create.log", append = TRUE)

  # clean up hashtag matrix to match the RNA data
  # rownames(adt_mat) = str_sub(rownames(adt_mat), 1, -17)
  # adt_mat = adt_mat[sort(rownames(adt_mat)), ]
  # colnames(adt_mat) = str_c(sample_name, ":", colnames(adt_mat))

  # Seurat replaces "_" or "|" in feature names with "-"
  # rownames(adt_mat) = str_replace(rownames(adt_mat), "_", "-")

  # clean up counts matrix to match the RNA data
  common_cells = intersect(colnames(adt_mat), colnames(seurat_obj))
  if (length(common_cells) < 10) { stop("cell names do not match expression matrix") }
  adt_mat = adt_mat[, common_cells]

  message(glue("RNA library {sample_name} cells: {ncol(seurat_obj)}"))
  message(glue("RNA and ADT library {sample_name} common cells: {ncol(adt_mat)}"))
  write(glue("RNA and ADT library {sample_name} common cells: {ncol(adt_mat)}"), file = "create.log", append = TRUE)

  # create a matrix that includes all cells from the original seurat object (fill 0 for missing cells)
  if (ncol(adt_mat) < ncol(seurat_obj)) {
    adt_filled_mat = matrix(data = 0, nrow = nrow(adt_mat), ncol = ncol(seurat_obj))
    colnames(adt_filled_mat) = colnames(seurat_obj)
    rownames(adt_filled_mat) = rownames(adt_mat)
    adt_filled_mat[, common_cells] = adt_mat[, common_cells]
  }

  # add ADT data as a new assay independent from RNA
  # seurat_obj[["ADT"]] = CreateAssayObject(counts = adt_filled_mat)

  # rename nCount_RNA and nFeature_RNA slots to make them more clear
  seurat_obj$num_ADT_UMIs = seurat_obj$nCount_ADT
  seurat_obj$num_ADT_genes = seurat_obj$nFeature_ADT

  message("\n\n ========== normalize ADT data ========== \n\n")

  # normalize ADT data using centered log-ratio (CLR) transformation
  seurat_obj = NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
  seurat_obj = ScaleData(seurat_obj, assay = "ADT")

  # save raw ADT counts matrix as a text file
  counts_raw = GetAssayData(seurat_obj, assay = "ADT", slot = "counts") %>% as.matrix()
  counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  fwrite(counts_raw, file = "counts.adt.raw.csv.gz", sep = ",", nThread = 4)

  # save normalized ADT counts matrix as a text file
  counts_norm = GetAssayData(seurat_obj, assay = "ADT") %>% as.matrix() %>% round(3)
  counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  fwrite(counts_norm, file = "counts.adt.normalized.csv.gz", sep = ",", nThread = 4)

  # plot ADT metrics
  plot_adt_qc(seurat_obj = seurat_obj)

  seurat_obj = scooter::set_identity(seurat_obj = seurat_obj, identity_column = "orig.ident")

  return(seurat_obj)

}

# plot ADT metrics
plot_adt_qc = function(seurat_obj) {

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = seurat_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  suppressMessages({
    dist_adt_g_plot =
      VlnPlot(
        seurat_obj, features = "num_ADT_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_adt_u_plot =
      VlnPlot(
        seurat_obj, features = "num_ADT_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_plot = plot_grid(dist_adt_g_plot, dist_adt_u_plot, ncol = 2)
    ggsave("qc.adt.distribution.png", plot = dist_plot, width = 14, height = 6, units = "in")
  })
  Sys.sleep(1)

  cor_umis_plot =
    FeatureScatter(
      seurat_obj, feature1 = "num_genes", feature2 = "num_ADT_genes",
      group.by = "orig.ident", cols = colors_samples_named
    ) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  cor_genes_plot =
    FeatureScatter(
      seurat_obj, feature1 = "num_UMIs", feature2 = "num_ADT_UMIs",
      group.by = "orig.ident", cols = colors_samples_named
    ) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  cor_mito_plot =
    FeatureScatter(
      seurat_obj, feature1 = "pct_mito", feature2 = "num_ADT_UMIs",
      group.by = "orig.ident", cols = colors_samples_named
    ) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  cor_adt_plot = plot_grid(cor_umis_plot, cor_genes_plot, cor_mito_plot, ncol = 3)
  ggsave("qc.adt.correlations.png", plot = cor_adt_plot, width = 18, height = 5, units = "in")
  Sys.sleep(1)

  # generate ADT counts summary
  counts_raw = GetAssayData(seurat_obj, assay = "ADT", slot = "counts") %>% as.matrix()
  adt_counts_summary = rowSums(counts_raw) %>% enframe(name = "ADT", value = "total_counts")
  adt_counts_summary$mean_counts = rowMeans(counts_raw) %>% round(1)
  adt_counts_summary$median_counts = matrixStats::rowMedians(counts_raw)
  adt_counts_summary$q05_counts = matrixStats::rowQuantiles(counts_raw, probs = 0.05)
  adt_counts_summary$q95_counts = matrixStats::rowQuantiles(counts_raw, probs = 0.95)
  write_csv(adt_counts_summary, "qc.adt.counts.csv")

}

# add hashtag oligo (HTO) data to a Seurat object
# add_hto_assay = function(seurat_obj, sample_name, sample_dir) {
split_adt_hto_assay = function(seurat_obj) {

  # check that the input object already has UMAP computed
  if (is.null(seurat_obj@reductions$umap)) { stop("UMAP not computed yet") }

  message("\n\n ========== import hashtag oligos (HTOs) ========== \n\n")

  # message("loading HTO matrix for sample: ", sample_name)
  sample_name = seurat_obj@meta.data$orig.ident %>% as.character() %>% unique()

  # check if sample dir is valid
  # if (!dir.exists(sample_dir)) { stop(glue("dir {sample_dir} does not exist")) }
  # if (!file.exists(glue("{sample_dir}/matrix.mtx.gz"))) { stop(glue("dir does not contain matrix.mtx.gz")) }

  # attempt to split ADTs matrix into ADTs and HTOs
  adt_hto_features = rownames(seurat_obj[["ADT"]])
  adt_features = adt_hto_features %>% str_subset("^CD|IgG|ADT") %>% unique() %>% sort()
  # if there are no features that fit ADT criteria, assume they are all HTOs
  if (length(adt_features) > 0) {
    # hto_features = setdiff(adt_hto_features, adt_features)
    hto_features = adt_hto_features %>% str_subset("HTO") %>% unique() %>% sort()
  } else {
    hto_features = adt_hto_features
  }
  

  message("detected ADTs: ", str_c(adt_features, collapse = ", "))
  message("detected HTOs: ", str_c(hto_features, collapse = ", "))

  # check for potential problems
  if (length(c(adt_features, hto_features)) < length(adt_hto_features)) { stop("missing ADTs/HTOs detected") }
  if (length(intersect(adt_features, hto_features)) > 0) { stop("conflicting ADTs/HTOs detected") }
  if (length(hto_features) < 2) { stop("no HTOs detected") }

  # hto_mat = Read10X(data.dir = sample_dir, gene.column = 1)
  hto_mat = GetAssayData(seurat_obj, assay = "ADT", slot = "counts")
  hto_mat = as.matrix(hto_mat)
  hto_mat = hto_mat[hto_features, ]
  rownames(hto_mat) = str_remove(rownames(hto_mat), "^HTO-")

  # removed "unmapped" HTO
  # if (rownames(hto_mat)[length(rownames(hto_mat))] == "unmapped") {
  #   hto_mat = hto_mat[1:length(rownames(hto_mat))-1, ]
  # }

  message(glue("HTO library {sample_name} unfiltered cells: {ncol(hto_mat)}"))
  message(glue("HTO library {sample_name} unfiltered HTOs: {nrow(hto_mat)}"))

  # log to file
  write(glue("HTO library {sample_name} unfiltered cells: {ncol(hto_mat)}"), file = "create.log", append = TRUE)
  write(glue("HTO library {sample_name} unfiltered HTOs: {nrow(hto_mat)}"), file = "create.log", append = TRUE)

  # clean up hashtag matrix to match the RNA data
  # rownames(hto_mat) = str_sub(rownames(hto_mat), 1, -17)
  # hto_mat = hto_mat[sort(rownames(hto_mat)), ]
  # colnames(hto_mat) = str_c(sample_name, ":", colnames(hto_mat))

  # Seurat replaces "_" or "|" in feature names with "-"
  # rownames(hto_mat) = str_replace(rownames(hto_mat), "_", "-")

  # clean up counts matrix to match the RNA data
  common_cells = intersect(colnames(hto_mat), colnames(seurat_obj))
  if (length(common_cells) < 10) { stop("cell names do not match expression matrix") }
  hto_mat = hto_mat[, common_cells]

  message(glue("RNA library {sample_name} cells: {ncol(seurat_obj)}"))
  message(glue("RNA and HTO library {sample_name} common cells: {ncol(hto_mat)}"))
  write(glue("RNA and HTO library {sample_name} common cells: {ncol(hto_mat)}"), file = "create.log", append = TRUE)

  # log the detailed unfiltered HTO stats to file
  hto_col_sums = colSums(hto_mat)
  write(glue("unfiltered min HTO reads: {min(hto_col_sums)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered max HTO reads: {max(hto_col_sums)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered mean HTO reads: {round(mean(hto_col_sums), 3)}"), file = "create.log", append = TRUE)
  write(glue("unfiltered median HTO reads: {median(hto_col_sums)}"), file = "create.log", append = TRUE)

  # compute the deviation from the median of each cell (constant is 1.48 for normal distribution)
  hto_col_sums = log1p(hto_col_sums)
  hto_col_mads = (hto_col_sums - median(hto_col_sums)) / mad(hto_col_sums, constant = 1)

  # remove outlier cells based on total HTO counts
  hto_outliers = hto_col_mads > 3
  hto_outliers = names(hto_outliers[hto_outliers == TRUE])
  seurat_obj = subset(seurat_obj, cells = hto_outliers, invert = TRUE)
  hto_mat = hto_mat[, colnames(seurat_obj)]

  message(glue("HTO library {sample_name} filtered cells: {ncol(hto_mat)}"))
  message(glue("HTO library {sample_name} filtered HTOs: {nrow(hto_mat)}"))

  write(glue("HTO library {sample_name} filtered cells: {ncol(hto_mat)}"), file = "create.log", append = TRUE)
  write(glue("HTO library {sample_name} filtered HTOs: {nrow(hto_mat)}"), file = "create.log", append = TRUE)

  # log the detailed filtered HTO stats to file
  hto_col_sums = colSums(hto_mat)
  write(glue("filtered min HTO reads: {min(hto_col_sums)}"), file = "create.log", append = TRUE)
  write(glue("filtered max HTO reads: {max(hto_col_sums)}"), file = "create.log", append = TRUE)
  write(glue("filtered mean HTO reads: {round(mean(hto_col_sums), 3)}"), file = "create.log", append = TRUE)
  write(glue("filtered median HTO reads: {median(hto_col_sums)}"), file = "create.log", append = TRUE)

  # log the new RNA stats to file
  write(glue("filtered mean num genes: {round(mean(seurat_obj$num_genes), 3)}"), file = "create.log", append = TRUE)
  write(glue("filtered median num genes: {median(seurat_obj$num_genes)}"), file = "create.log", append = TRUE)

  # add HTO data as a new assay independent from RNA
  seurat_obj[["HTO"]] = CreateAssayObject(counts = hto_mat)

  # normalize HTO data using centered log-ratio (CLR) transformation
  seurat_obj = NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")

  # assign single cells back to their sample origins
  # seurat_obj = HTODemux(seurat_obj, assay = "HTO", kfunc = "kmeans", positive.quantile = 0.999)
  seurat_obj = HTODemux(seurat_obj, assay = "HTO", kfunc = "kmeans")

  # plot HTO metrics
  plot_hto_qc(seurat_obj = seurat_obj)

  # save HTO stats
  seurat_obj = scooter::set_identity(seurat_obj = seurat_obj, identity_column = "hash.ID")
  calculate_cluster_stats(seurat_obj = seurat_obj, label = "hto")

  # update metadata, setting the hashtag as the sample name
  seurat_obj@meta.data$library = factor(seurat_obj@meta.data$orig.ident)
  seurat_obj@meta.data$orig.ident = factor(seurat_obj@meta.data$hash.ID)

  seurat_obj = scooter::set_identity(seurat_obj = seurat_obj, identity_column = "orig.ident")

  return(seurat_obj)

}

# plot HTO metrics
plot_hto_qc = function(seurat_obj) {

  # normalized HTO signal combined with metadata table
  hto_tbl = GetAssayData(seurat_obj, assay = "HTO") %>% t()
  hto_tbl = hto_tbl[, sort(colnames(hto_tbl))] %>% as_tibble(rownames = "cell")
  id_tbl = seurat_obj@meta.data %>% as_tibble(rownames = "cell") %>% select(cell, hash.ID, HTO_classification.global)
  hto_tbl = full_join(hto_tbl, id_tbl, by = "cell")
  hto_tbl

  # HTO color scheme
  colors_hto_names = c(levels(hto_tbl$HTO_classification.global), levels(hto_tbl$hash.ID)) %>% unique()
  colors_hto = colors_clusters[1:length(colors_hto_names)]
  names(colors_hto) = colors_hto_names

  # visualize pairs of HTO signals
  hto_facet_plot =
    ggplot(sample_frac(hto_tbl), aes(x = .panel_x, y = .panel_y, fill = hash.ID, color = hash.ID)) +
    geom_point(shape = 16, size = 0.2) +
    geom_autodensity(color = NA, fill = "gray20") +
    geom_density2d(color = "black", alpha = 0.5) +
    scale_color_manual(values = colors_hto) +
    scale_fill_manual(values = colors_hto) +
    facet_matrix(vars(-cell, -hash.ID, -HTO_classification.global), layer.diag = 2, layer.upper = 3) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(plot.background = element_rect(fill = "white"), aspect.ratio = 1, legend.title = element_blank(), strip.background = element_blank())
  save_plot(filename = "qc.hto.correlation.png", plot = hto_facet_plot, base_height = 8, base_width = 10)
  Sys.sleep(1)
  save_plot(filename = "qc.hto.correlation.pdf", plot = hto_facet_plot, base_height = 8, base_width = 10)
  Sys.sleep(1)

  # number of UMIs for singlets, doublets and negative cells
  hto_umi_plot =
    VlnPlot(seurat_obj, features = "num_UMIs", group.by = "HTO_classification.global", pt.size = 0.1) +
    theme(plot.background = element_rect(fill = "white"), axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = colors_hto) +
    scale_y_continuous(labels = comma)
  save_plot("qc.hto.umis.png", plot = hto_umi_plot, base_height = 6, base_width = 6)
  Sys.sleep(1)

  # number of genes for singlets, doublets and negative cells
  hto_gene_plot =
    VlnPlot(seurat_obj, features = "num_genes", group.by = "HTO_classification.global", pt.size = 0.1) +
    theme(plot.background = element_rect(fill = "white"), axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = colors_hto) +
    scale_y_continuous(labels = comma)
  save_plot("qc.hto.genes.png", plot = hto_gene_plot, base_height = 6, base_width = 6)
  Sys.sleep(1)

  # number of genes for singlets, doublets and negative cells
  hto_mito_plot =
    VlnPlot(seurat_obj, features = "pct_mito", group.by = "HTO_classification.global", pt.size = 0.1) +
    theme(plot.background = element_rect(fill = "white"), axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = colors_hto) +
    scale_y_continuous(labels = comma)
  save_plot("qc.hto.mito.png", plot = hto_mito_plot, base_height = 6, base_width = 6)
  Sys.sleep(1)

  group_var = "HTO_classification.global"
  seurat_obj = scooter::set_identity(seurat_obj = seurat_obj, identity_column = group_var)
  plot_umap =
    DimPlot(
      seurat_obj, reduction = "umap",
      pt.size = get_dr_point_size(seurat_obj), cols = colors_hto, shuffle = TRUE, raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  save_plot(glue("dr.umap.{group_var}.png"), plot = plot_umap, base_height = 6, base_width = 8)
  Sys.sleep(1)
  save_plot(glue("dr.umap.{group_var}.pdf"), plot = plot_umap, base_height = 6, base_width = 8)
  Sys.sleep(1)

  group_var = "hash.ID"
  seurat_obj = scooter::set_identity(seurat_obj = seurat_obj, identity_column = group_var)
  plot_umap =
    DimPlot(
      seurat_obj, reduction = "umap",
      pt.size = get_dr_point_size(seurat_obj), cols = colors_hto, shuffle = TRUE, raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  save_plot(glue("dr.umap.{group_var}.png"), plot = plot_umap, base_height = 6, base_width = 8)
  Sys.sleep(1)
  save_plot(glue("dr.umap.{group_var}.pdf"), plot = plot_umap, base_height = 6, base_width = 8)
  Sys.sleep(1)

  if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

}

# split Seurat object (by sample by default)
split_seurat_obj = function(seurat_obj, original_wd, split_var = "orig.ident") {

  # set identity to the column used for splitting
  s_obj = seurat_obj
  DefaultAssay(s_obj) = "RNA"
  s_obj = scooter::set_identity(seurat_obj = s_obj, identity_column = split_var)

  # clean up metadata
  s_obj@assays$RNA@var.features = vector()
  s_obj@assays$RNA@scale.data = matrix()
  s_obj@reductions = list()
  s_obj@meta.data = s_obj@meta.data %>% select(-starts_with("snn_res."))
  s_obj@meta.data = s_obj@meta.data %>% select(-starts_with("res."))

  setwd(original_wd)

  # process each split group
  split_names = Idents(s_obj) %>% as.character() %>% unique() %>% sort()
  for (s in split_names) {

    message(glue("split: {s}"))
    split_obj = subset(s_obj, idents = s)
    message(glue("subset {s} cells: {ncol(split_obj)}"))
    if (ncol(split_obj) > 10) {
      split_dir = glue("split-{s}")
      if (dir.exists(split_dir)) {
        stop(glue("output dir {split_dir} already exists"))
      } else {
        dir.create(split_dir)
        setwd(split_dir)
        write(glue("subset {s} cells: {ncol(split_obj)}"), file = "create.log", append = TRUE)
        split_obj = calculate_variance(seurat_obj = split_obj, jackstraw_max_cells = 100)
        split_obj = scooter::set_identity(seurat_obj = split_obj, identity_column = "orig.ident")
        saveRDS(split_obj, file = "seurat_obj.rds")
        calculate_cluster_stats(split_obj, label = "sample")
      }
    }

    # clean up and return to the main dir before processing the next split
    if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
    setwd(original_wd)

  }

}

# merge multiple Seurat objects
combine_seurat_obj = function(original_wd, sample_analysis_dirs) {

  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")

  message("\n\n ========== combine samples ========== \n\n")

  seurat_obj_list = list()
  for (i in 1:length(sample_analysis_dirs)) {

    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
    sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")

    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))

    # load seurat object
    seurat_obj_list[[i]] = readRDS(sample_seurat_rds)

    # clean up object
    seurat_obj_list[[i]]@assays$RNA@var.features = vector()
    seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
    seurat_obj_list[[i]]@reductions = list()
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res."))
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("res."))

    # print single sample sample stats
    sample_name = seurat_obj_list[[i]]$orig.ident %>% as.character() %>% sort()
    sample_name = sample_name[1]
    message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
    write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(" ")

  }

  # merge
  merged_obj = merge(seurat_obj_list[[1]], seurat_obj_list[2:length(seurat_obj_list)])
  rm(seurat_obj_list)

  # print combined sample stats
  message(glue("combined unfiltered cells: {ncol(merged_obj)}"))
  write(glue("combined unfiltered cells: {ncol(merged_obj)}"), file = "create.log", append = TRUE)
  message(glue("combined unfiltered genes: {nrow(merged_obj)}"))
  write(glue("combined unfiltered genes: {nrow(merged_obj)}"), file = "create.log", append = TRUE)

  # filter poorly expressed genes (detected in less than 10 cells)
  filtered_genes = Matrix::rowSums(GetAssayData(merged_obj, assay = "RNA", slot = "counts") > 0)
  min_cells = 10
  if (ncol(merged_obj) > 50000) { min_cells = ncol(merged_obj) * 0.001 }
  filtered_genes = filtered_genes[filtered_genes >= min_cells] %>% names() %>% sort()
  # keep all HTO and ADT features if present
  if ("HTO" %in% names(merged_obj@assays)) { filtered_genes = c(filtered_genes, rownames(merged_obj@assays$HTO)) }
  if ("ADT" %in% names(merged_obj@assays)) { filtered_genes = c(filtered_genes, rownames(merged_obj@assays$ADT)) }
  merged_obj = subset(merged_obj, features = filtered_genes)

  # encode sample name as factor (also sets alphabetical sample order)
  merged_obj@meta.data$orig.ident = factor(merged_obj@meta.data$orig.ident)
  merged_obj = scooter::set_identity(seurat_obj = merged_obj, identity_column = "orig.ident")

  # print combined sample stats
  message(glue("combined cells: {ncol(merged_obj)}"))
  write(glue("combined cells: {ncol(merged_obj)}"), file = "create.log", append = TRUE)
  message(glue("combined genes: {nrow(merged_obj)}"))
  write(glue("combined genes: {nrow(merged_obj)}"), file = "create.log", append = TRUE)

  # print gene/cell minimum cutoffs
  min_cells = Matrix::rowSums(GetAssayData(merged_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  min_genes = Matrix::colSums(GetAssayData(merged_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  message(glue("min cells per gene: {min_cells}"))
  write(glue("min cells per gene: {min_cells}"), file = "create.log", append = TRUE)
  message(glue("min genes per cell: {min_genes}"))
  write(glue("min genes per cell: {min_genes}"), file = "create.log", append = TRUE)

  # check that the full counts table is small enough to fit into an R matrix (max around 100k x 21k)
  num_matrix_elements = GetAssayData(merged_obj, assay = "RNA", slot = "counts") %>% length()
  if (num_matrix_elements < 2^31) {

    # save raw counts matrix as a text file
    counts_raw = GetAssayData(merged_obj, assay = "RNA", slot = "counts") %>% as.matrix()
    counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    fwrite(counts_raw, file = "counts.raw.csv.gz", sep = ",", nThread = 4)

    # save normalized counts matrix as a text file
    counts_norm = GetAssayData(merged_obj, assay = "RNA") %>% as.matrix() %>% round(3)
    counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    fwrite(counts_norm, file = "counts.normalized.csv.gz", sep = ",", nThread = 4)

  }

  # create a named color scheme to ensure names and colors are in the proper order
  sample_names = merged_obj$orig.ident %>% as.character() %>% sort() %>% unique()
  colors_samples_named = colors_samples[1:length(sample_names)]
  names(colors_samples_named) = sample_names

  vln_theme =
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  suppressMessages({
    dist_nft_plot =
      VlnPlot(
        merged_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_nct_plot =
      VlnPlot(
        merged_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_pmt_plot =
      VlnPlot(
        merged_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_plot = plot_grid(dist_nft_plot, dist_nct_plot, dist_pmt_plot, ncol = 3)
    ggsave("qc.distribution.png", plot = dist_plot, width = 20, height = 6, units = "in")
  })
  Sys.sleep(1)

  return(merged_obj)

}

# integrate multiple Seurat objects
integrate_seurat_obj = function(original_wd, sample_analysis_dirs, num_dim, int_reduction = "cca") {

  # check if the inputs seems reasonable
  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")
  num_dim = as.integer(num_dim)
  if (num_dim < 5) { stop("too few dims: ", num_dim) }
  if (num_dim > 50) { stop("too many dims: ", num_dim) }

  message("\n\n ========== integrate samples ========== \n\n")

  seurat_obj_list = list()
  var_genes_list = list()
  exp_genes = c()
  for (i in 1:length(sample_analysis_dirs)) {

    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
    sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")

    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))

    # load seurat object
    seurat_obj_list[[i]] = readRDS(sample_seurat_rds)
    sample_name = seurat_obj_list[[i]]$orig.ident %>% as.character() %>% sort()
    sample_name = sample_name[1]

    # clean up object
    seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
    seurat_obj_list[[i]]@reductions = list()
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res."))
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("res."))

    # save expressed genes keeping only genes present in all the datasets (for genes to integrate in IntegrateData)
    if (length(exp_genes) > 0) {
      exp_genes = intersect(exp_genes, rownames(seurat_obj_list[[i]])) %>% sort()
    } else {
      exp_genes = rownames(seurat_obj_list[[i]])
    }

    # save variable genes
    var_genes_list[[sample_name]] = VariableFeatures(seurat_obj_list[[i]])

    # print single sample sample stats
    message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
    write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(" ")

  }

  # euler plot of variable gene overlaps (becomes unreadable and can take days for many overlaps)
  if (length(var_genes_list) < 8) {
    colors_euler = colors_samples[1:length(var_genes_list)]
    euler_fit = euler(var_genes_list, shape = "ellipse")
    euler_plot = plot(euler_fit,
                      fills = list(fill = colors_euler, alpha = 0.7),
                      edges = list(col = colors_euler))
    png("variance.vargenes.euler.png", res = 200, width = 5, height = 5, units = "in")
      print(euler_plot)
    dev.off()
  }

  # upset plot of variable gene overlaps
  upset_plot = upset(fromList(var_genes_list), nsets = 50, nintersects = 15, order.by = "freq", mb.ratio = c(0.5, 0.5))
  png("variance.vargenes.upset.png", res = 200, width = 8, height = 5, units = "in")
    print(upset_plot)
  dev.off()

  # RPCA workflow requires users to run PCA on each dataset prior to integration
  if (int_reduction == "rpca") {
    message("\n\n ========== Seurat::SelectIntegrationFeatures() ========== \n\n")
    # select features that are repeatedly variable across datasets for integration run PCA on each dataset
    int_features = SelectIntegrationFeatures(object.list = seurat_obj_list)
    message(glue("integration features: {length(int_features)}"))
    seurat_obj_list =
      lapply(seurat_obj_list, FUN = function(x) {
        x = ScaleData(x, features = int_features, vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)
        x = RunPCA(x, features = int_features, verbose = FALSE)
      })
  }

  message("\n\n ========== Seurat::FindIntegrationAnchors() ========== \n\n")

  # find the integration anchors
  if (int_reduction == "cca") {
    anchors = FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = 2000,
      reduction = int_reduction, dims = 1:num_dim)
  } else if (int_reduction == "rpca") {
    # k.anchor: how many neighbors (k) to use when picking anchors (default: 5)
    # "You can increase the strength of alignment by increasing the k.anchor parameter."
    # "Increasing this parameter to 20 will assist in aligning these populations."
    anchors = FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = 2000,
      reduction = int_reduction, dims = 1:num_dim, k.anchor = 10)
  } else {
    stop("unknown reduction")
  }
  rm(seurat_obj_list)

  message("\n\n ========== Seurat::IntegrateData() ========== \n\n")

  # integrating all genes may cause issues and may not add any relevant information
  integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim)
  rm(anchors)

  # after running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix
  # switch to integrated assay
  DefaultAssay(integrated_obj) = "integrated"

  # print integrated sample stats
  message(glue("integrated unfiltered cells: {ncol(integrated_obj)}"))
  write(glue("integrated unfiltered cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
  message(glue("integrated unfiltered genes: {nrow(integrated_obj)}"))
  write(glue("integrated unfiltered genes: {nrow(integrated_obj)}"), file = "create.log", append = TRUE)

  # filter poorly expressed genes (detected in less than 10 cells)
  filtered_genes = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0)
  min_cells = 10
  if (ncol(integrated_obj) > 50000) { min_cells = ncol(integrated_obj) * 0.001 }
  filtered_genes = filtered_genes[filtered_genes >= min_cells] %>% names() %>% sort()
  # keep all HTO and ADT features if present
  if ("HTO" %in% names(integrated_obj@assays)) { filtered_genes = c(filtered_genes, rownames(integrated_obj@assays$HTO)) }
  if ("ADT" %in% names(integrated_obj@assays)) { filtered_genes = c(filtered_genes, rownames(integrated_obj@assays$ADT)) }
  integrated_obj = subset(integrated_obj, features = filtered_genes)

  # encode sample name as factor (also sets alphabetical sample order)
  integrated_obj@meta.data$orig.ident = factor(integrated_obj@meta.data$orig.ident)

  # print integrated sample stats
  message(glue("integrated cells: {ncol(integrated_obj)}"))
  write(glue("integrated cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
  message(glue("integrated genes: {nrow(GetAssayData(integrated_obj, assay = 'RNA'))}"))
  write(glue("integrated genes: {nrow(GetAssayData(integrated_obj, assay = 'RNA'))}"), file = "create.log", append = TRUE)

  # print gene/cell minumum cutoffs
  min_cells = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  min_genes = Matrix::colSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  message(glue("min cells per gene: {min_cells}"))
  write(glue("min cells per gene: {min_cells}"), file = "create.log", append = TRUE)
  message(glue("min genes per cell: {min_genes}"))
  write(glue("min genes per cell: {min_genes}"), file = "create.log", append = TRUE)

  # check that the full counts table is small enough to fit into an R matrix (max around 100k x 21k)
  num_matrix_elements = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% length()
  if (num_matrix_elements < 2^31) {

    # save raw counts matrix as a text file
    counts_raw = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% as.matrix()
    counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    fwrite(counts_raw, file = "counts.raw.csv.gz", sep = ",", nThread = 4)

    # save normalized counts matrix as a text file
    counts_norm = GetAssayData(integrated_obj, assay = "RNA") %>% as.matrix() %>% round(3)
    counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    fwrite(counts_norm, file = "counts.normalized.csv.gz", sep = ",", nThread = 4)

  }

  vln_theme =
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  suppressMessages({
    dist_nft_plot =
      VlnPlot(
        integrated_obj, features = "num_genes", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_nct_plot =
      VlnPlot(
        integrated_obj, features = "num_UMIs", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_pmt_plot =
      VlnPlot(
        integrated_obj, features = "pct_mito", group.by = "orig.ident",
        pt.size = 0.1, sort = TRUE, combine = TRUE, cols = colors_samples
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_plot = plot_grid(dist_nft_plot, dist_nct_plot, dist_pmt_plot, ncol = 3)
    ggsave("qc.distribution.png", plot = dist_plot, width = 20, height = 6, units = "in")
  })
  Sys.sleep(1)

  return(integrated_obj)

}

# calculate various variance metrics and perform basic analysis
# PC selection approaches:
# - PCHeatmap - more supervised, exploring PCs to determine relevant sources of heterogeneity
# - PCElbowPlot - heuristic that is commonly used and can be calculated instantly
# - JackStrawPlot - implements a statistical test based on a random null model, but is time-consuming
# jackStraw procedure is very slow, so skip for large projects (>10,000 cells)
calculate_variance = function(seurat_obj, jackstraw_max_cells = 10000) {

  s_obj = seurat_obj

  message("\n\n ========== Seurat::FindVariableGenes() ========== \n\n")

  # identify features that are outliers on a 'mean variability plot'
  # Seurat v3 implements an improved method based on a variance stabilizing transformation ("vst")
  s_obj = FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  # export highly variable feature information (mean, variance, variance standardized)
  hvf_tbl = HVFInfo(s_obj) %>% round(3) %>% rownames_to_column("gene") %>% arrange(-variance.standardized)
  write_csv(hvf_tbl, "variance.csv")

  # plot variance
  var_plot = VariableFeaturePlot(s_obj, pt.size = 0.5)
  var_plot = LabelPoints(var_plot, points = head(hvf_tbl$gene, 30), repel = TRUE, xnudge = 0, ynudge = 0)
  ggsave("variance.features.png", plot = var_plot, width = 12, height = 5, units = "in")

  message("\n\n ========== Seurat::ScaleData() ========== \n\n")

  # regress out unwanted sources of variation
  # regressing uninteresting sources of variation can improve dimensionality reduction and clustering
  # could include technical noise, batch effects, biological sources of variation (cell cycle stage)
  # scaled z-scored residuals of these models are stored in scale.data slot
  # used for dimensionality reduction and clustering
  s_obj = ScaleData(s_obj, vars.to.regress = c("num_UMIs", "pct_mito"), verbose = FALSE)

  message("\n\n ========== Seurat::PCA() ========== \n\n")

  # use fewer PCs for small datasets
  num_pcs = 50
  if (ncol(s_obj) < 100) num_pcs = 20
  if (ncol(s_obj) < 25) num_pcs = 5

  # PCA on the scaled data
  # PCA calculation stored in object[["pca"]]
  s_obj = RunPCA(s_obj, assay = "RNA", features = VariableFeatures(s_obj), npcs = num_pcs, verbose = FALSE)

  # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
  pca_plot =
    DimPlot(
      s_obj, group.by = "orig.ident", reduction = "pca",
      pt.size = 0.5, cols = colors_samples, shuffle = TRUE, raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  ggsave("variance.pca.png", plot = pca_plot, width = 8, height = 6, units = "in")

  message("\n\n ========== Seurat::DimHeatmap() ========== \n\n")

  # PCHeatmap (former) allows for easy exploration of the primary sources of heterogeneity in a dataset
  if (num_pcs > 15) {
    png("variance.pca.heatmap.png", res = 300, width = 10, height = 16, units = "in")
      DimHeatmap(s_obj, reduction = "pca", dims = 1:15, nfeatures = 20, cells = 250, fast = TRUE)
    dev.off()
  }

  message("\n\n ========== Seurat::PCElbowPlot() ========== \n\n")

  # a more ad hoc method for determining PCs to use, draw cutoff where there is a clear elbow in the graph
  elbow_plot = ElbowPlot(s_obj, reduction = "pca", ndims = num_pcs)
  ggsave("variance.pca.elbow.png", plot = elbow_plot, width = 8, height = 5, units = "in")

  # resampling test inspired by the jackStraw procedure - very slow, so skip for large projects (>10,000 cells)
  if (ncol(s_obj) < jackstraw_max_cells) {

    message("\n\n ========== Seurat::JackStraw() ========== \n\n")

    # determine statistical significance of PCA scores
    s_obj = JackStraw(s_obj, assay = "RNA", reduction = "pca", dims = num_pcs, verbose = FALSE)

    # compute Jackstraw scores significance
    s_obj = ScoreJackStraw(s_obj, reduction = "pca", dims = 1:num_pcs, do.plot = FALSE)

    # plot the results of the JackStraw analysis for PCA significance
    # significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line)
    jackstraw_plot =
      JackStrawPlot(s_obj, reduction = "pca", dims = 1:num_pcs) +
      guides(col = guide_legend(ncol = 2))
    ggsave("variance.pca.jackstraw.png", plot = jackstraw_plot, width = 12, height = 6, units = "in")

  }

  return(s_obj)

}

# calculate various variance metrics and perform basic analysis (integrated analysis workflow)
# specify neighbors for UMAP (default is 30 in Seurat 2 and 3 pre-release)
calculate_variance_integrated = function(seurat_obj, num_dim, num_neighbors = 30) {

  s_obj = seurat_obj

  num_dim = as.integer(num_dim)
  if (num_dim < 5) { stop("too few dims: ", num_dim) }
  if (num_dim > 50) { stop("too many dims: ", num_dim) }

  message("\n\n ========== Seurat::ScaleData() ========== \n\n")

  # s_obj = ScaleData(s_obj, features = rownames(s_obj), verbose = FALSE)
  s_obj = ScaleData(s_obj, verbose = FALSE)

  message("\n\n ========== Seurat::PCA() ========== \n\n")

  # PCA on the scaled data
  s_obj = RunPCA(s_obj, npcs = num_dim, verbose = FALSE)

  # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
  pca_plot =
    DimPlot(
      s_obj, reduction = "pca", group.by = "orig.ident",
      pt.size = 0.5, cols = colors_samples, shuffle = TRUE, raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  ggsave("variance.pca.png", plot = pca_plot, width = 10, height = 6, units = "in")

  message("\n\n ========== Seurat::RunTSNE() ========== \n\n")

  # use tSNE as a tool to visualize, not for clustering directly on tSNE components
  # cells within the graph-based clusters determined above should co-localize on the tSNE plot
  s_obj = RunTSNE(s_obj, reduction = "pca", dims = 1:num_dim, dim.embed = 2)

  # reduce point size for larger datasets
  dr_pt_size = get_dr_point_size(s_obj)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = scooter::set_identity(seurat_obj = s_obj, identity_column = "orig.ident")
  plot_tsne =
    DimPlot(s_obj, reduction = "tsne", pt.size = dr_pt_size, cols = colors_samples, shuffle = TRUE, raster = FALSE) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  ggsave(glue("dr.tsne.{num_dim}.sample.png"), plot = plot_tsne, width = 10, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.tsne.{num_dim}.sample.pdf"), plot = plot_tsne, width = 10, height = 6, units = "in")
  Sys.sleep(1)

  message("\n\n ========== Seurat::RunUMAP() ========== \n\n")

  # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  s_obj = RunUMAP(s_obj, reduction = "pca", dims = 1:num_dim, n.neighbors = num_neighbors, verbose = FALSE)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = scooter::set_identity(seurat_obj = s_obj, identity_column = "orig.ident")
  plot_umap =
    DimPlot(s_obj, reduction = "umap", pt.size = dr_pt_size, cols = colors_samples, shuffle = TRUE, raster = FALSE) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  ggsave(glue("dr.umap.{num_dim}.sample.png"), plot = plot_umap, width = 10, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.umap.{num_dim}.sample.pdf"), plot = plot_umap, width = 10, height = 6, units = "in")
  Sys.sleep(1)

  save_metadata(seurat_obj = s_obj)

  return(s_obj)

}

# perform graph-based clustering and tSNE
# specify neighbors for UMAP and FindNeighbors (default is 30 in Seurat 2 and 3 pre-release)
calculate_clusters = function(seurat_obj, num_dim, num_neighbors = 30) {

  # check if number of dimensions seems reasonable
  if (num_dim < 5) { stop("too few dims: ", num_dim) }
  if (num_dim > 50) { stop("too many dims: ", num_dim) }

  s_obj = seurat_obj

  message("\n\n ========== Seurat::RunTSNE() ========== \n\n")

  # use tSNE as a tool to visualize, not for clustering directly on tSNE components
  # cells within the graph-based clusters determined above should co-localize on the tSNE plot
  s_obj = RunTSNE(s_obj, reduction = "pca", dims = 1:num_dim, dim.embed = 2)

  # reduce point size for larger datasets
  dr_pt_size = get_dr_point_size(s_obj)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = scooter::set_identity(seurat_obj = s_obj, identity_column = "orig.ident")
  plot_tsne =
    DimPlot(s_obj, reduction = "tsne", pt.size = dr_pt_size, cols = colors_samples, shuffle = TRUE, raster = FALSE) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  ggsave(glue("dr.tsne.{num_dim}.sample.png"), plot = plot_tsne, width = 10, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.tsne.{num_dim}.sample.pdf"), plot = plot_tsne, width = 10, height = 6, units = "in")
  Sys.sleep(1)

  message("\n\n ========== Seurat::RunUMAP() ========== \n\n")

  # runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  s_obj = RunUMAP(s_obj, reduction = "pca", dims = 1:num_dim, n.neighbors = num_neighbors, verbose = FALSE)

  # tSNE using original sample names (shuffle cells so any one group does not appear overrepresented due to ordering)
  s_obj = scooter::set_identity(seurat_obj = s_obj, identity_column = "orig.ident")
  plot_umap =
    DimPlot(s_obj, reduction = "umap", pt.size = dr_pt_size, cols = colors_samples, shuffle = TRUE, raster = FALSE) +
    theme(
      plot.background = element_rect(fill = "white"), aspect.ratio = 1,
      axis.ticks = element_blank(), axis.text = element_blank()
    )
  ggsave(glue("dr.umap.{num_dim}.sample.png"), plot = plot_umap, width = 10, height = 6, units = "in")
  Sys.sleep(1)
  ggsave(glue("dr.umap.{num_dim}.sample.pdf"), plot = plot_umap, width = 10, height = 6, units = "in")
  Sys.sleep(1)

  message("\n\n ========== Seurat::FindNeighbors() ========== \n\n")

  message("assay: ", DefaultAssay(s_obj))
  message("num dims: ", num_dim)

  # construct a Shared Nearest Neighbor (SNN) Graph for a given dataset
  s_obj =
    FindNeighbors(
      s_obj, dims = 1:num_dim, k.param = num_neighbors, compute.SNN = TRUE, force.recalc = TRUE
    )

  message("\n\n ========== Seurat::FindClusters() ========== \n\n")

  message("initial metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))

  # resolutions for graph-based clustering
  # increased resolution values lead to more clusters (recommendation: 0.6-1.2 for 3K cells, 2-4 for 33K cells)
  res_range = seq(0.1, 2.0, 0.1)
  if (ncol(s_obj) > 10000) res_range = c(res_range, 3, 5, 7, 9)

  # algorithm: 1 = original Louvain; 2 = Louvain with multilevel refinement; 3 = SLM
  # identify clusters of cells by SNN modularity optimization based clustering algorithm
  s_obj = FindClusters(s_obj, algorithm = 3, resolution = res_range, verbose = FALSE)

  # simplify clustering column names to match previous style (just "res")
  colnames(s_obj@meta.data) = str_replace(colnames(s_obj@meta.data), ".*nn_res\\.", "res.")

  # remove "seurat_clusters" column that is added automatically (added in v3 late dev version)
  s_obj@meta.data = s_obj@meta.data %>% select(-seurat_clusters)

  message("new metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))

  # create a separate sub-directory for cluster resolution plots
  clusters_dir = "clusters-resolutions"
  if (!dir.exists(clusters_dir)) { dir.create(clusters_dir) }

  # for calculated cluster resolutions: remove redundant (same number of clusters), rename, and plot
  res_cols = str_subset(colnames(s_obj@meta.data), "^res\\.")
  res_cols = sort(res_cols)
  res_num_clusters_prev = 1
  for (res in res_cols) {

    # proceed if current resolution has more clusters than previous and less than the color scheme length
    res_vector = s_obj@meta.data[, res] %>% as.character()
    res_num_clusters_cur = res_vector %>% n_distinct()
    if (res_num_clusters_cur > res_num_clusters_prev && res_num_clusters_cur < length(colors_clusters)) {

      # check if the resolution still has original labels (characters starting with 0)
      if (min(res_vector) == "0") {

        # convert to character vector
        s_obj@meta.data[, res] = as.character(s_obj@meta.data[, res])
        # relabel identities so they start with 1 and not 0
        s_obj@meta.data[, res] = as.numeric(s_obj@meta.data[, res]) + 1
        # pad with 0s to avoid sorting issues
        s_obj@meta.data[, res] = str_pad(s_obj@meta.data[, res], width = 2, side = "left", pad = "0")
        # pad with "C" to avoid downstream numeric conversions
        s_obj@meta.data[, res] = str_c("C", s_obj@meta.data[, res])
        # encode as a factor
        s_obj@meta.data[, res] = factor(s_obj@meta.data[, res])

      }

      # resolution value based on resolution column name
      res_val = str_remove(res, "^res\\.")

      # plot file name
      res_str = format(as.numeric(res_val), nsmall = 1)
      res_str = str_remove(res_str, "\\.")
      res_str = str_pad(res_str, width = 3, side = "left", pad = "0")
      res_str = str_c("res", res_str)
      dr_filename = glue("{clusters_dir}/dr.{DefaultAssay(s_obj)}.{num_dim}.{res_str}.clust{res_num_clusters_cur}")

      s_obj = plot_clusters(seurat_obj = s_obj, resolution = res_val, filename_base = dr_filename)

      # add blank line to make output easier to read
      message(" ")

    } else {

      # remove resolution if the number of clusters is same as previous
      s_obj@meta.data = s_obj@meta.data %>% select(-one_of(res))

    }

    # update resolution cluster count for next iteration
    res_num_clusters_prev = res_num_clusters_cur

  }

  message("updated metadata fields: ", str_c(colnames(s_obj@meta.data), collapse = ", "))

  save_metadata(seurat_obj = s_obj)

  return(s_obj)

}

# compile all cell metadata into a single table
save_metadata = function(seurat_obj) {

  s_obj = seurat_obj
  metadata_tbl = s_obj@meta.data %>% as_tibble(rownames = "cell") %>% mutate(sample_name = orig.ident)
  tsne_tbl = s_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  umap_tbl = s_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
  metadata_tbl = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
  metadata_tbl = metadata_tbl %>% arrange(cell)
  write_csv(metadata_tbl, "metadata.csv")

}

# plot tSNE with color-coded clusters at specified resolution
plot_clusters = function(seurat_obj, resolution, filename_base) {

  s_obj = seurat_obj

  # set identities based on specified resolution
  s_obj = scooter::set_identity(seurat_obj = s_obj, identity_column = resolution)

  # print stats
  num_clusters = Idents(s_obj) %>% as.character() %>% n_distinct()
  message("resolution: ", resolution)
  message("num clusters: ", num_clusters)

  # generate plot if there is a reasonable number of clusters
  if (num_clusters > 1 && num_clusters < length(colors_clusters)) {

    # shuffle cells so they appear randomly and one group does not show up on top
    plot_tsne =
      DimPlot(
        s_obj, reduction = "tsne",
        pt.size = get_dr_point_size(s_obj), cols = colors_clusters, shuffle = TRUE, raster = FALSE
      ) +
      theme(
        plot.background = element_rect(fill = "white"), aspect.ratio = 1,
        axis.ticks = element_blank(), axis.text = element_blank()
      )
    ggsave(glue("{filename_base}.tsne.png"), plot = plot_tsne, width = 9, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("{filename_base}.tsne.pdf"), plot = plot_tsne, width = 9, height = 6, units = "in")
    Sys.sleep(1)

    plot_umap =
      DimPlot(
        s_obj, reduction = "umap",
        pt.size = get_dr_point_size(s_obj), cols = colors_clusters, shuffle = TRUE, raster = FALSE
      ) +
      theme(
        plot.background = element_rect(fill = "white"), aspect.ratio = 1,
        axis.ticks = element_blank(), axis.text = element_blank()
      )
    ggsave(glue("{filename_base}.umap.png"), plot = plot_umap, width = 9, height = 6, units = "in")
    Sys.sleep(1)
    ggsave(glue("{filename_base}.umap.pdf"), plot = plot_umap, width = 9, height = 6, units = "in")
    Sys.sleep(1)

    if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

  }

  return(s_obj)

}

# plot marker genes in multiple ways
plot_cluster_markers_top = function(seurat_obj, genes, filename_base) {

  # color gradient for FeaturePlot-based plots
  # gradient_colors = c("gray85", "red2")
  gradient_colors = colorRampPalette(c("#d9cfcb", "#d49070", "#ca5528", "#b72600", "#981000", "#730000"))(100)

  # switch to "RNA" assay from potentially "integrated"
  DefaultAssay(seurat_obj) = "RNA"

  # UMAP plots color-coded by expression level (should be square to match the original tSNE plots)
  feat_plot =
    FeaturePlot(
      seurat_obj, features = genes, reduction = "umap", cells = sample(colnames(seurat_obj)),
      pt.size = 0.5, cols = gradient_colors, ncol = 4
    )
  ggsave(glue("{filename_base}.umap.png"), plot = feat_plot, width = 16, height = 10, units = "in", dpi = 150)
  # ggsave(glue("{filename_base}.umap.pdf"), plot = feat_plot, width = 16, height = 10, units = "in")

  # dot plot visualization
  dot_plot = DotPlot(seurat_obj, features = genes, dot.scale = 12, cols = gradient_colors)
  ggsave(glue("{filename_base}.dotplot.png"), plot = dot_plot, width = 20, height = 8, units = "in")
  ggsave(glue("{filename_base}.dotplot.pdf"), plot = dot_plot, width = 20, height = 8, units = "in")

  # gene violin plots (size.use below 0.2 doesn't seem to make a difference)
  # skip PDF since every cell has to be plotted and they become too big
  vln_plot = VlnPlot(seurat_obj, features = genes, pt.size = 0.1, combine = TRUE, cols = colors_clusters, ncol = 4)
  ggsave(glue("{filename_base}.violin.png"), plot = vln_plot, width = 16, height = 10, units = "in", dpi = 150)

  # expression levels per cluster for bar plots (averaging and output are in non-log space)
  cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", features = genes, verbose = FALSE)[["RNA"]]
  cluster_avg_exp_long = cluster_avg_exp %>% as.data.frame() %>% rownames_to_column("gene") %>% gather(cluster, avg_exp, -gene)

  # bar plots
  # create a named color scheme to ensure names and colors are in the proper order
  clust_names = levels(seurat_obj)
  color_scheme_named = colors_clusters[1:length(clust_names)]
  names(color_scheme_named) = clust_names
  barplot_plot = ggplot(cluster_avg_exp_long, aes(x = cluster, y = avg_exp, fill = cluster)) +
    geom_col(color = "black") +
    theme(plot.background = element_rect(fill = "white"), legend.position = "none") +
    scale_fill_manual(values = color_scheme_named) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot() +
    facet_wrap(~ gene, ncol = 4, scales = "free")
  ggsave(glue("{filename_base}.barplot.png"), plot = barplot_plot, width = 20, height = 10, units = "in", dpi = 150)
  # ggsave(glue("{filename_base}.barplot.pdf"), plot = barplot_plot, width = 20, height = 10, units = "in")

}

# plot gene expression overlaid on a UMAP based on a table of genes and groups
plot_dr_umap_genes = function(seurat_obj, genes_csv) {

  # check that the input object already has UMAP computed
  if (is.null(seurat_obj@reductions$umap)) { stop("UMAP not computed yet") }

  # import genes table and check that the "gene" and "group" columns exist
  if (!file.exists(genes_csv)) { stop(glue("genes table {genes_csv} does not exist")) }
  genes_tbl = read_csv(genes_csv, col_types = cols())
  if (!is.element("gene", colnames(genes_tbl))) { stop("gene table column names must include 'gene'") }
  if (!is.element("group", colnames(genes_tbl))) { stop("gene table column names must include 'group'") }

  # switch to RNA assay
  DefaultAssay(seurat_obj) = "RNA"

  # check for detected genes
  genes_tbl = distinct(genes_tbl)
  missing_genes = setdiff(genes_tbl$gene, rownames(seurat_obj))
  genes_tbl = genes_tbl %>% filter(gene %in% rownames(seurat_obj))
  message("\n gene groups: ", str_c(unique(genes_tbl$group), collapse = ", "))
  message("\n detectable genes: ", str_c(genes_tbl$gene, collapse = ", "))
  message("\n missing genes: ", str_c(missing_genes, collapse = ", "))
  if (nrow(genes_tbl) == 0) { stop("no detectable genes") }
  Sys.sleep(1)

  # plot settings
  DefaultAssay(seurat_obj) = "RNA"
  dr_point_size = get_dr_point_size(seurat_obj)
  gradient_colors = colorRampPalette(c("#d9cfcb", "#d49070", "#ca5528", "#b72600", "#981000", "#730000"))(100)

  # plot each gene in a separate directory based on group
  for (gene_group in unique(genes_tbl$group)) {
    genes_group_tbl = genes_tbl %>% filter(group == gene_group)
    plot_dir = glue("dr-umap-genes-{gene_group}")
    dir.create(plot_dir)
    for (g in sort(genes_group_tbl$gene)) {
      gene_umap =
        FeaturePlot(
          seurat_obj, features = g, reduction = "umap",
          cells = sample(colnames(seurat_obj)), pt.size = dr_point_size, cols = gradient_colors, raster = FALSE
        ) +
        labs(title = g) +
        theme(
          plot.background = element_rect(fill = "white"),
          aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
          axis.ticks = element_blank(), axis.text = element_blank()
        )
      save_plot(glue("{plot_dir}/dr.umap.gene.{g}.png"), plot = gene_umap, base_height = 6.5, base_width = 7)
      Sys.sleep(1)
      # save_plot(glue("{plot_dir}/dr.umap.gene.{g}.pdf"), plot = gene_umap, base_height = 6.5, base_width = 7)
      # Sys.sleep(1)
    }
  }

}

# gather metadata and calculate cluster stats (number of cells)
calculate_cluster_stats = function(seurat_obj, label) {

  message("\n\n ========== calculate cluster stats ========== \n\n")

  message("cluster names: ", str_c(levels(seurat_obj), collapse = ", "))

  # compile relevant cell metadata into a single table
  seurat_obj$cluster = Idents(seurat_obj)
  metadata_tbl = as_tibble(seurat_obj@meta.data, rownames = "cell")
  metadata_tbl = select(metadata_tbl, cell, num_UMIs, num_genes, pct_mito, sample_name = orig.ident, cluster)
  if (!is.null(seurat_obj@reductions$tsne)) {
    tsne_tbl = seurat_obj[["tsne"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
    umap_tbl = seurat_obj[["umap"]]@cell.embeddings %>% round(3) %>% as.data.frame() %>% rownames_to_column("cell")
    metadata_tbl = metadata_tbl %>% full_join(tsne_tbl, by = "cell") %>% full_join(umap_tbl, by = "cell")
    metadata_tbl = metadata_tbl %>% arrange(cell)
    write_csv(metadata_tbl, glue("metadata.{label}.csv"))
  }
  Sys.sleep(1)

  # get number of cells split by cluster and by sample (orig.ident)
  summary_cluster_sample =
    metadata_tbl %>%
    select(cluster, sample_name) %>%
    mutate(num_cells_total = n()) %>%
    group_by(sample_name) %>%
    mutate(num_cells_sample = n()) %>%
    group_by(cluster) %>%
    mutate(num_cells_cluster = n()) %>%
    group_by(cluster, sample_name) %>%
    mutate(num_cells_cluster_sample = n()) %>%
    ungroup() %>%
    distinct() %>%
    mutate(
      pct_cells_cluster = num_cells_cluster / num_cells_total,
      pct_cells_cluster_sample = num_cells_cluster_sample / num_cells_sample
    ) %>%
    mutate(
      pct_cells_cluster = round(pct_cells_cluster * 100, 1),
      pct_cells_cluster_sample = round(pct_cells_cluster_sample * 100, 1)
    ) %>%
    arrange(cluster, sample_name)

  # get number of cells split by cluster (ignore samples)
  summary_cluster = summary_cluster_sample %>% select(-contains("sample")) %>% distinct()
  write_csv(summary_cluster, glue("summary.{label}.csv"))
  Sys.sleep(1)

  # export results split by sample if multiple samples per cluster are present
  num_samples = metadata_tbl %>% pull(sample_name) %>% n_distinct()
  if (nrow(summary_cluster_sample) > nrow(summary_cluster)) {
    write_csv(summary_cluster_sample, glue("summary.{label}.per-sample.csv"))
    Sys.sleep(1)
  }

}

# calculate cluster average expression (non-log space)
calculate_cluster_expression = function(seurat_obj, label) {

  message("\n\n ========== calculate cluster average expression ========== \n\n")

  message("cluster names: ", str_c(levels(seurat_obj), collapse = ", "))

  # gene expression for an "average" cell in each identity class (averaging and output are in non-log space)
  cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", verbose = FALSE)[["RNA"]]
  cluster_avg_exp = cluster_avg_exp %>% round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(cluster_avg_exp, glue("expression.mean.{label}.csv"))
  Sys.sleep(1)

  # cluster averages split by sample (orig.ident)
  num_samples = n_distinct(seurat_obj@meta.data$orig.ident)
  if (num_samples > 1) {
    seurat_obj@meta.data$persample = paste0(Idents(seurat_obj), "|", seurat_obj@meta.data$orig.ident)
    sample_avg_exp = AverageExpression(seurat_obj, assay = "RNA", group.by = "persample", verbose = FALSE)[["RNA"]]
    sample_avg_exp = sample_avg_exp %>% round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    write_csv(sample_avg_exp, glue("expression.mean.{label}.per-sample.csv"))
    Sys.sleep(1)
  }

}

# calculate cluster markers (compared to all other cells) and plot top ones
# tests:
# - roc: ROC test returns the classification power (ranging from 0 - random, to 1 - perfect)
# - wilcox: Wilcoxon rank sum test (default in Seurat 2)
# - bimod: Likelihood-ratio test for single cell gene expression (McDavid, Bioinformatics, 2013) (default in Seurat 1)
# - tobit: Tobit-test for differential gene expression (Trapnell, Nature Biotech, 2014)
# - MAST: GLM-framework that treates cellular detection rate as a covariate (Finak, Genome Biology, 2015)
# pairwise option compares each cluster to each of the other clusters to yield markers that are both local and global
calculate_cluster_markers = function(seurat_obj, label, test, pairwise = FALSE) {

  message("\n\n ========== calculate cluster markers ========== \n\n")

  message("cluster set: ", label)
  message("marker test: ", test)

  # get cluster names
  clusters = Idents(seurat_obj) %>% as.character() %>% unique() %>% sort()

  # use only clusters with more than 10 cells
  clusters = clusters[table(Idents(seurat_obj)) > 10]

  if (!pairwise) {

    # standard cluster markers calculation

    markers_dir = "markers-global"

    # capture output to avoid excessive warnings
    # default logfc.threshold = 0.25 = log2(1.189), log2(1.1) = 0.137
    markers_log =
      capture.output({
        all_markers =
          FindAllMarkers(
            seurat_obj, assay = "RNA", test.use = test,
            min.pct = 0.1, logfc.threshold = log2(1.1), base = 2, fc.name = "log2FC",
            only.pos = TRUE, min.diff.pct = -Inf, verbose = FALSE
          )
      }, type = "message")

    # do some light filtering and clean up (ROC test returns slightly different output)
    if (test == "roc") {

      all_markers =
        all_markers %>%
        select(cluster, gene, log2FC = avg_diff, myAUC, power) %>%
        filter(power > 0.3) %>%
        mutate(log2FC = round(log2FC, 5), myAUC = round(myAUC, 5), power = round(power, 5)) %>%
        arrange(cluster, -power)
      top_markers = all_markers %>% filter(log2FC > 0)
      top_markers = top_markers %>% group_by(cluster) %>% top_n(50, power) %>% ungroup()

    } else {

      all_markers =
        all_markers %>%
        select(cluster, gene, log2FC, p_val, p_val_adj) %>%
        filter(p_val_adj < 0.001) %>%
        mutate(log2FC = round(log2FC, 5)) %>%
        arrange(cluster, p_val_adj, p_val)
      top_markers = all_markers %>% filter(log2FC > 0)
      top_markers = top_markers %>% group_by(cluster) %>% top_n(50, log2FC) %>% ungroup()

    }

  } else {

    # pairwise (each cluster versus each other cluster) cluster markers calculation

    markers_dir = "markers-pairwise"

    # initialize empty results tibble
    unfiltered_markers = tibble(
      cluster = character(),
      cluster2 = character(),
      gene = character(),
      log2FC = numeric(),
      p_val = numeric(),
      p_val_adj = numeric()
    )

    # check each cluster combination
    for (cluster1 in clusters) {
      for (cluster2 in setdiff(clusters, cluster1)) {

        # find differentially expressed genes between two specific clusters
        # low fold change cutoff to maximize chance of appearing in all comparisons
        # capture output to avoid excessive warnings
        markers_log =
          capture.output({
            cur_markers =
              FindMarkers(
                seurat_obj, assay = "RNA", ident.1 = cluster1, ident.2 = cluster2, test.use = test,
                min.pct = 0.1, logfc.threshold = 0, base = 2, fc.name = "log2FC",
                only.pos = TRUE, min.diff.pct = -Inf, verbose = FALSE
              )
          }, type = "message")

        # clean up markers table (would need to be modified for "roc" test)
        cur_markers =
          cur_markers %>%
          rownames_to_column("gene") %>%
          mutate(cluster = cluster1) %>%
          mutate(cluster2 = cluster2) %>%
          filter(p_val_adj < 0.01) %>%
          mutate(log2FC = round(log2FC, 5)) %>%
          select(one_of(colnames(unfiltered_markers)))

        # add current cluster combination genes to the table of all markers
        unfiltered_markers = bind_rows(unfiltered_markers, cur_markers)

      }
    }

    # adjust test name for output
    test = glue("pairwise.{test}")

    # sort the markers to make the table more readable
    unfiltered_markers =
      unfiltered_markers %>%
      distinct() %>%
      add_count(cluster, gene) %>%
      rename(cluster_gene_n = n) %>%
      arrange(cluster, gene, cluster2)

    # filter for genes that are significant compared to all other clusters
    all_markers =
      unfiltered_markers %>%
      filter(cluster_gene_n == (length(clusters) - 1)) %>%
      select(-cluster_gene_n)

    # extract the lowest and highest fold changes and p-values
    all_markers =
      all_markers %>%
      group_by(cluster, gene) %>%
      summarize_at(
        c("log2FC", "p_val", "p_val_adj"),
        list(min = min, max = max)
      ) %>%
      ungroup() %>%
      arrange(cluster, -log2FC_min)

    top_markers = all_markers %>% group_by(cluster) %>% top_n(50, log2FC_min) %>% ungroup()

  }

  # create a separate sub-directory for all markers
  if (!dir.exists(markers_dir)) { dir.create(markers_dir) }

  # filename prefix
  filename_base = glue("{markers_dir}/markers.{label}.{test}")

  # save unfiltered markers for pairwise comparisons
  if (pairwise) {
    unfiltered_markers_csv = glue("{filename_base}.unfiltered.csv")
    message("unfiltered markers: ", unfiltered_markers_csv)
    write_csv(unfiltered_markers, unfiltered_markers_csv)
    Sys.sleep(1)
  }

  all_markers_csv = glue("{filename_base}.all.csv")
  message("all markers: ", all_markers_csv)
  write_csv(all_markers, all_markers_csv)
  Sys.sleep(1)

  top_markers_csv = glue("{filename_base}.top.csv")
  message("top markers: ", top_markers_csv)
  write_csv(top_markers, top_markers_csv)
  Sys.sleep(1)

  # plot cluster markers heatmap
  plot_cluster_markers_heatmap(seurat_obj, markers_tbl = all_markers, num_genes = c(5, 10, 20), filename_base = filename_base)

  # plot top cluster markers for each cluster
  for (cluster_name in clusters) {
    filename_cluster_base = glue("{markers_dir}/markers.{label}-{cluster_name}.{test}")
    cluster_markers = top_markers %>% filter(cluster == cluster_name)
    if (nrow(cluster_markers) > 9) {
      Sys.sleep(1)
      top_cluster_markers = cluster_markers %>% head(12) %>% pull(gene)
      plot_cluster_markers_top(seurat_obj, genes = top_cluster_markers, filename_base = filename_cluster_base)
    }

  }

}

# generate cluster markers heatmap
plot_cluster_markers_heatmap = function(seurat_obj, markers_tbl, num_genes, filename_base) {

  # adjust pairwise clusters to match the standard format
  if ("log2FC_min" %in% colnames(markers_tbl)) {
    markers_tbl = markers_tbl %>% mutate(log2FC = log2FC_min)
  }

  # keep only the top cluster for each gene so each gene appears once
  markers_tbl = markers_tbl %>% filter(log2FC > 0)
  markers_tbl = markers_tbl %>% group_by(gene) %>% top_n(1, log2FC) %>% slice(1) %>% ungroup()

  num_clusters = Idents(seurat_obj) %>% as.character() %>% n_distinct()
  marker_genes = markers_tbl %>% pull(gene) %>% unique() %>% sort()
  cluster_avg_exp = AverageExpression(seurat_obj, assay = "RNA", features = marker_genes, verbose = FALSE)[["RNA"]]
  cluster_avg_exp = cluster_avg_exp %>% as.matrix() %>% log1p()
  cluster_avg_exp = cluster_avg_exp[rowSums(cluster_avg_exp) > 0, ]

  # heatmap settings
  hm_colors = colorRampPalette(c("#053061", "#FFFFFF", "#E41A1C"))(51)
  hm_width = ( num_clusters / 2 ) + 2

  for (ng in num_genes) {

    hm_base = glue("{filename_base}.heatmap.top{ng}")

    markers_top_tbl = markers_tbl %>% group_by(cluster) %>% top_n(ng, log2FC) %>% ungroup()
    markers_top_tbl = markers_top_tbl %>% arrange(cluster, -log2FC)

    # generate the scaled expression matrix and save the text version
    hm_mat = cluster_avg_exp[markers_top_tbl$gene, ]
    hm_mat = hm_mat %>% t() %>% scale() %>% t()
    hm_mat %>% round(3) %>% as_tibble(rownames = "gene") %>% write_csv(glue("{hm_base}.csv"))
    Sys.sleep(1)

    # set outliers to 95th percentile to yield a more balanced color scale
    scale_cutoff = as.numeric(quantile(abs(hm_mat), 0.95))
    hm_mat[hm_mat > scale_cutoff] = scale_cutoff
    hm_mat[hm_mat < -scale_cutoff] = -scale_cutoff

    # generate the heatmap
    ph_obj = pheatmap(
      hm_mat, scale = "none", color = hm_colors, border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE,
      fontsize = 10, fontsize_row = 8, fontsize_col = 12, show_colnames = TRUE,
      main = glue("Cluster Markers: Top {ng}")
    )

    png(glue("{hm_base}.png"), width = hm_width, height = 10, units = "in", res = 300)
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


# ========== main ==========


# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# default type for the bitmap devices such as png (should default to "cairo")
options(bitmapType = "cairo")

# retrieve the command-line arguments
suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)

# show docopt options
# print(opts)

# dependencies
load_libraries()

# set number of cores for parallel package (will use all available cores by default)
options(mc.cores = 4)
# evaluate Seurat R expressions asynchronously when possible (such as ScaleData) using future package
plan("multiprocess", workers = 4)
# increase the limit of the data to be shuttled between the processes from default 500MB to 50GB
options(future.globals.maxSize = 50e9)
# disable future seed warning (introduced in future 1.19.0, should be fixed in Seurat 4)
options(future.rng.onMisuse = "ignore")

# global settings
colors_samples = scooter::get_color_scheme("samples")
colors_clusters = scooter::get_color_scheme("clusters")

# analysis info
analysis_step = "unknown"
out_dir = opts$analysis_dir

# original working dir (before moving to analysis dir)
original_wd = getwd()

# check input files convert to canonical form in case they are relative
if (!is.null(opts$sample_dir)) {
  if (dir.exists(opts$sample_dir)) {
    opts$sample_dir = normalizePath(opts$sample_dir)
  } else {
    stop(glue("sample dir {opts$sample_dir} does not exist"))
  }
}
if (!is.null(opts$genes_csv)) {
  if (file.exists(opts$genes_csv)) {
    opts$genes_csv = normalizePath(opts$genes_csv)
  } else {
    stop(glue("genes table {opts$genes_csv} does not exist"))
  }
}

# create analysis directory if starting new analysis or exit if analysis already exists
if (opts$create || opts$combine || opts$integrate) {

  if (opts$create) analysis_step = "create"
  if (opts$combine) analysis_step = "combine"
  if (opts$integrate) analysis_step = "integrate"
  message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

  if (dir.exists(out_dir)) {
    stop(glue("output analysis dir {out_dir} already exists"))
  } else {
    dir.create(out_dir)
  }

}

# set analysis directory as working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  stop(glue("output analysis dir {out_dir} does not exist"))
}

# check which command was used
if (opts$create) {

  # log to file
  write(glue("analysis: {out_dir}"), file = "create.log", append = TRUE)
  write(glue("seurat version: {packageVersion('Seurat')}"), file = "create.log", append = TRUE)
  write(glue("scooter version: {packageVersion('scooter')}"), file = "create.log", append = TRUE)

  # create new seurat object based on input sample names and sample directories
  counts_mat = scooter::load_sample_counts_matrix(opts$sample_name, opts$sample_dir)
  write(glue("counts matrix cells: {ncol(counts_mat[['Gene Expression']])}"), file = "create.log", append = TRUE)
  write(glue("counts matrix genes: {nrow(counts_mat[['Gene Expression']])}"), file = "create.log", append = TRUE)
  # create seurat object using gene expression data
  seurat_obj = scooter::create_seurat_obj(counts_matrix = counts_mat[["Gene Expression"]], assay = "RNA")
  # rename nCount_RNA and nFeature_RNA slots to make them more clear
  seurat_obj$num_UMIs = seurat_obj$nCount_RNA
  seurat_obj$num_genes = seurat_obj$nFeature_RNA
  # QC plots
  create_seurat_obj_qc(seurat_obj)
  # add ADT data to seurat object
  if ("Antibody Capture" %in% names(counts_mat)) {
    seurat_obj = scooter::add_seurat_assay(seurat_obj, assay = "ADT", counts_matrix = counts_mat[["Antibody Capture"]])
    add_adt_assay_qc(seurat_obj, sample_name = opts$sample_name)
  }

  # filter by number of genes and mitochondrial genes percentage (optional parameters)
  seurat_obj = filter_data(seurat_obj, min_genes = opts$min_genes, max_genes = opts$max_genes, max_mt = opts$mt)

  # calculate various variance metrics
  seurat_obj = calculate_variance(seurat_obj, jackstraw_max_cells = 1000)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

} else if (opts$combine) {

  write(glue("seurat version: {packageVersion('Seurat')}"), file = "create.log", append = TRUE)

  # merge multiple samples/libraries based on previous analysis directories
  seurat_obj = combine_seurat_obj(original_wd = original_wd, sample_analysis_dirs = opts$sample_analysis_dir)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

  # calculate various variance metrics
  seurat_obj = calculate_variance(seurat_obj)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

  # save sample stats
  calculate_cluster_stats(seurat_obj, label = "sample")

} else if (opts$integrate) {

  write(glue("seurat version: {packageVersion('Seurat')}"), file = "create.log", append = TRUE)
  write(glue("integration reduction: {opts$reduction}"), file = "create.log", append = TRUE)

  # run integration
  if (opts$reduction %in% c("cca", "rpca")) {
    seurat_obj = integrate_seurat_obj(original_wd, sample_analysis_dirs = opts$batch_analysis_dir, num_dim = opts$num_dim, int_reduction = opts$reduction)
  } else {
     stop(glue("integration reduction type {opts$reduction} is not valid"))
  }

  seurat_obj = calculate_variance_integrated(seurat_obj, num_dim = opts$num_dim)

  saveRDS(seurat_obj, file = "seurat_obj.rds")

  # save sample stats
  calculate_cluster_stats(seurat_obj, label = "sample")

} else {

  # all commands besides "create", "combine", and "integrate" start with an existing seurat object
  if (file.exists("seurat_obj.rds")) {

    message("loading seurat_obj")
    seurat_obj = readRDS("seurat_obj.rds")

  } else {

    stop("seurat obj does not already exist (run 'create' step first)")

  }

  if (opts$cluster) {

    analysis_step = "cluster"
    message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

    # determine clusters
    seurat_obj = calculate_clusters(seurat_obj, num_dim = as.integer(opts$num_dim))
    saveRDS(seurat_obj, file = "seurat_obj.rds")

  }

  if (opts$hto) {

    analysis_step = "hto"
    message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

    # add HTO data
    seurat_obj = split_adt_hto_assay(seurat_obj)
    seurat_obj@meta.data$library_hash = str_c(seurat_obj@meta.data$library, "-", seurat_obj@meta.data$hash.ID)
    saveRDS(seurat_obj, file = "seurat_obj.rds")

    # split singlets by HTO sample
    Idents(seurat_obj) = "HTO_classification.global"
    seurat_obj = subset(seurat_obj, idents = "Singlet")
    split_seurat_obj(seurat_obj, original_wd = original_wd, split_var = "library_hash")

  }

  if (opts$plot) {

    if (opts$umap) {

      analysis_step = "plot umap"
      message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

      plot_dr_umap_genes(seurat_obj = seurat_obj, genes_csv = opts$genes_csv)

    } else {

      stop("unknown plot type (only umap is currently supported)")

    }

  }

  if (opts$identify || opts$de) {

    # set resolution in the seurat object
    seurat_obj = scooter::set_identity(seurat_obj = seurat_obj, identity_column = opts$resolution)

    # use a grouping-specific sub-directory for all output
    grouping_label = scooter::check_identity_column(seurat_obj, opts$resolution)
    grouping_label = gsub("\\.", "", grouping_label)
    num_clusters = Idents(seurat_obj) %>% as.character() %>% n_distinct()
    clust_label = glue("clust{num_clusters}")
    res_dir = glue("clusters-{grouping_label}-{clust_label}")
    if (!dir.exists(res_dir)) dir.create(res_dir)
    setwd(res_dir)

    if (opts$identify) {

      analysis_step = "identify"
      message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

      # create tSNE plot (should already exist in the main directory)
      dr_filename = glue("dr.{grouping_label}.{clust_label}")
      seurat_obj = plot_clusters(seurat_obj, resolution = opts$resolution, filename_base = dr_filename)

      # cluster stat tables (number of cells and average expression)
      calculate_cluster_stats(seurat_obj, label = clust_label)
      calculate_cluster_expression(seurat_obj, label = clust_label)

      # calculate and plot standard cluster markers
      # calculate_cluster_markers(seurat_obj, label = clust_label, test = "roc")
      calculate_cluster_markers(seurat_obj, label = clust_label, test = "wilcox")
      # calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST")

      # calculate and plot pairwise cluster markers (very slow, so skip for high number of clusters)
      num_clusters = Idents(seurat_obj) %>% as.character() %>% n_distinct()
      if (num_clusters < 20) {
        calculate_cluster_markers(seurat_obj, label = clust_label, test = "wilcox", pairwise = TRUE)
        # calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST", pairwise = TRUE)
      }

    }

    # differential expression
    if (opts$de) {

        analysis_step = "diff"
        message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

        # calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "wilcox", group_var = opts$group_var)
        # calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "MAST", group_var = opts$group_var)
        de_tbl = differential_expression_per_cluster(seurat_obj, cluster_column = opts$resolution, group_column = opts$group_var, test = "wilcox")

    }

  }

}

message(glue("\n\n ========== finished analysis step {analysis_step} for {out_dir} ========== \n\n"))

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
