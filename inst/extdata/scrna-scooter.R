#!/usr/bin/env Rscript

"
Command-line pipeline for analysis of single cell RNA-seq data using Seurat (version >=5.3.1).
The input is Cell Ranger output directory with filtered feature-barcode matrices in MEX or HDF5 format.

Basic steps:
  create - import counts, filter cells by QC metrics, normalize, and compute an initial PCA/tSNE/UMAP
  merge - combine multiple samples/libraries (no batch correction)
  integrate - combine multiple samples/libraries (batch correction)
  cluster - cluster cells across a resolution sweep, at a given number of PCs
  identify - compute stats, expression, and marker genes for clusters at a chosen resolution

Advanced steps:
  de - differential expression between samples/libraries within clusters
  adt - add antibody-derived tag data [retired]
  hto - add hashtag oligo data and split by hashtag
  plot umap - generate gene expression overlaid on a UMAP based on a table of genes and groups

Usage:
  scrna-scooter.R create <analysis_dir> --sample_name=<x> --sample_dir=<x> [--num_mads=<n> --min_genes=<n> --max_genes=<n> --min_counts=<n> --max_counts=<n> --mt=<n> --num_dim=<n>]
  scrna-scooter.R cluster <analysis_dir> --num_dim=<n> --reduction=<cluster_reduction>
  scrna-scooter.R identify <analysis_dir> --cluster_var=<x>
  scrna-scooter.R merge <analysis_dir> <sample_analysis_dir>...
  scrna-scooter.R integrate <analysis_dir> --batch_var=<x> --reduction=<x> [--num_dim=<n>]
  scrna-scooter.R de <analysis_dir> --cluster_var=<x> --group_var=<x>
  scrna-scooter.R hto <analysis_dir>
  scrna-scooter.R plot umap <analysis_dir> (--genes_csv=<x> | --cluster_var=<x>)
  scrna-scooter.R --help

Create options:
  --sample_name=<x>  sample/library name that will be used for plots and tables
  --sample_dir=<x>   path to the feature-barcode matrix directory for the specific sample/library
  --num_mads=<n>     median absolute deviations for outlier removal, applied first [default: 3]
  --min_genes=<n>    cutoff for minimum number of genes per cell [default: 500]
  --max_genes=<n>    cutoff for maximum number of genes per cell
  --min_counts=<n>   cutoff for minimum number of counts per cell
  --max_counts=<n>   cutoff for maximum number of counts per cell
  --mt=<n>           cutoff for maximum percentage of mitochondrial reads per cell [default: 10]

Create/cluster/integrate options:
  --num_dim=<n>      number of dimensions to use; cluster: required; create: defaults to 30 if
                     omitted; integrate: defaults to 50 if omitted
  --reduction=<x>    for integrate: integration method (`cca`, `rpca`, or `harmony`); for
                     cluster: reduction to cluster on (`pca`, or `cca`/`rpca`/`harmony` after
                     integrate)
  --batch_var=<x>    batch variable (such as `orig.ident`), integrate only

Identify/DE/plot options:
  --cluster_var=<x>  cluster variable; used by identify, de, and plot umap
  --group_var=<x>    group variable (such as `orig.ident`) for differential expression comparisons
  --genes_csv=<x>    CSV file with `gene` and `group` columns (plots will be organized by group)

Options:
  -h, --help         show this screen
" -> doc


# ========== functions ==========

# load dependencies
load_libraries <- function() {
  message("\n\n ========== load libraries ========== \n\n")

  suppressPackageStartupMessages({
    library(glue)
    library(Seurat)
    library(future)
    library(Matrix)
    library(tidyverse)
    library(scales)
    library(data.table)
    library(cowplot)
    library(ggforce)
    # remotes::install_github("igordot/scooter")
    library(scooter)
  })

  theme_set(theme_cowplot())
}

# Normalize the ADT assay (CLR) and emit its QC and counts tables.
# scooter::create_seurat_object() already attached and cell-matched the
# assay, through scooter::add_seurat_assay(). So this function adds the one
# thing nothing else in the pipeline does: ADT normalization.
add_adt_assay_qc <- function(x, sample_name) {
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

  adt_mat <- GetAssayData(x, assay = "ADT", layer = "counts")

  scooter::write_message(
    glue("ADT library {sample_name} cells: {ncol(adt_mat)}"),
    log_file = "create.log"
  )
  scooter::write_message(
    glue("ADT library {sample_name} ADTs: {nrow(adt_mat)}"),
    log_file = "create.log"
  )

  # clean up hashtag matrix to match the RNA data
  # rownames(adt_mat) = str_sub(rownames(adt_mat), 1, -17)
  # adt_mat = adt_mat[sort(rownames(adt_mat)), ]
  # colnames(adt_mat) = str_c(sample_name, ":", colnames(adt_mat))

  # Seurat replaces "_" or "|" in feature names with "-"
  # rownames(adt_mat) = str_replace(rownames(adt_mat), "_", "-")

  # clean up counts matrix to match the RNA data
  common_cells <- intersect(colnames(adt_mat), colnames(x))
  if (length(common_cells) < 10) {
    stop("cell names do not match expression matrix")
  }
  adt_mat <- adt_mat[, common_cells]

  message(glue("RNA library {sample_name} cells: {ncol(x)}"))
  scooter::write_message(
    glue("RNA and ADT library {sample_name} common cells: {ncol(adt_mat)}"),
    log_file = "create.log"
  )

  # create a matrix that includes all cells from the original seurat object (fill 0 for missing cells)
  if (ncol(adt_mat) < ncol(x)) {
    adt_filled_mat <- matrix(data = 0, nrow = nrow(adt_mat), ncol = ncol(x))
    colnames(adt_filled_mat) <- colnames(x)
    rownames(adt_filled_mat) <- rownames(adt_mat)
    adt_filled_mat[, common_cells] <- adt_mat[, common_cells]
  }

  # add ADT data as a new assay independent from RNA
  # x[["ADT"]] = CreateAssayObject(counts = adt_filled_mat)

  # rename nCount_RNA and nFeature_RNA slots to make them more clear
  x$num_ADT_counts <- x$nCount_ADT
  x$num_ADT_genes <- x$nFeature_ADT

  message("\n\n ========== normalize ADT data ========== \n\n")

  # normalize ADT data using centered log-ratio (CLR) transformation
  x <- NormalizeData(x, assay = "ADT", normalization.method = "CLR")
  x <- ScaleData(x, assay = "ADT")

  # save raw ADT counts matrix as a text file
  counts_raw <- GetAssayData(x, assay = "ADT", layer = "counts") |> as.matrix()
  counts_raw <- counts_raw |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    arrange(gene)
  fwrite(counts_raw, file = "counts-adt-raw.csv.gz", sep = ",", nThread = 4)

  # save normalized ADT counts matrix as a text file
  counts_norm <- GetAssayData(x, assay = "ADT") |> as.matrix() |> round(3)
  counts_norm <- counts_norm |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    arrange(gene)
  fwrite(
    counts_norm,
    file = "counts-adt-normalized.csv.gz",
    sep = ",",
    nThread = 4
  )

  # plot ADT metrics
  plot_adt_qc(x = x)

  x <- scooter::set_identity(x = x, identity_column = "orig.ident")

  return(x)
}

# plot ADT metrics
plot_adt_qc <- function(x) {
  # create a named color scheme to ensure names and colors are in the proper order
  sample_names <- x$orig.ident |> as.character() |> sort() |> unique()
  colors_samples_named <- colors_samples[1:length(sample_names)]
  names(colors_samples_named) <- sample_names

  vln_theme <-
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  suppressMessages({
    dist_adt_g_plot <-
      VlnPlot(
        x,
        features = "num_ADT_genes",
        group.by = "orig.ident",
        layer = "counts",
        pt.size = 0.1,
        sort = TRUE,
        combine = TRUE,
        cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_adt_u_plot <-
      VlnPlot(
        x,
        features = "num_ADT_counts",
        group.by = "orig.ident",
        layer = "counts",
        pt.size = 0.1,
        sort = TRUE,
        combine = TRUE,
        cols = colors_samples_named
      ) +
      scale_y_continuous(labels = comma) +
      vln_theme
    dist_plot <- plot_grid(dist_adt_g_plot, dist_adt_u_plot, ncol = 2)
    ggsave(
      "qc-adt-distribution.png",
      plot = dist_plot,
      width = 14,
      height = 6,
      units = "in"
    )
  })
  Sys.sleep(1)

  cor_counts_plot <-
    FeatureScatter(
      x,
      feature1 = "detected_genes",
      feature2 = "num_ADT_genes",
      group.by = "orig.ident",
      cols = colors_samples_named
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5)
    )
  cor_genes_plot <-
    FeatureScatter(
      x,
      feature1 = "total_counts",
      feature2 = "num_ADT_counts",
      group.by = "orig.ident",
      cols = colors_samples_named
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5)
    )
  cor_mito_plot <-
    FeatureScatter(
      x,
      feature1 = "pct_mito",
      feature2 = "num_ADT_counts",
      group.by = "orig.ident",
      cols = colors_samples_named
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5)
    )
  cor_adt_plot <- plot_grid(
    cor_counts_plot,
    cor_genes_plot,
    cor_mito_plot,
    ncol = 3
  )
  ggsave(
    "qc-adt-correlations.png",
    plot = cor_adt_plot,
    width = 18,
    height = 5,
    units = "in"
  )
  Sys.sleep(1)

  # generate ADT counts summary
  counts_raw <- GetAssayData(x, assay = "ADT", layer = "counts") |> as.matrix()
  adt_counts_summary <- rowSums(counts_raw) |>
    enframe(name = "ADT", value = "total_counts")
  adt_counts_summary$mean_counts <- rowMeans(counts_raw) |> round(1)
  adt_counts_summary$median_counts <- matrixStats::rowMedians(counts_raw)
  adt_counts_summary$q01_counts <- matrixStats::rowQuantiles(
    counts_raw,
    probs = 0.01
  ) |>
    round(1)
  adt_counts_summary$q05_counts <- matrixStats::rowQuantiles(
    counts_raw,
    probs = 0.05
  ) |>
    round(1)
  adt_counts_summary$q95_counts <- matrixStats::rowQuantiles(
    counts_raw,
    probs = 0.95
  ) |>
    round(1)
  adt_counts_summary$q99_counts <- matrixStats::rowQuantiles(
    counts_raw,
    probs = 0.99
  ) |>
    round(1)
  write_csv(adt_counts_summary, "qc-adt-counts.csv")
}

# add hashtag oligo (HTO) data to a Seurat object
# add_hto_assay = function(seurat_obj, sample_name, sample_dir) {
split_adt_hto_assay <- function(x, method = c("HTODemux", "MULTIseqDemux")) {
  method <- match.arg(method)

  # check that the input object already has UMAP computed
  if (is.null(x@reductions$umap)) {
    stop("UMAP not computed yet")
  }

  # plot ADT metrics (repeated here in case the names were manually modified for demultiplexing)
  plot_adt_qc(x = x)

  # reset orig.ident if demultiplexing has already been performed
  if (n_distinct(x@meta.data$orig.ident) > 1) {
    if ("hash.ID" %in% colnames(x@meta.data)) {
      if ("library" %in% colnames(x@meta.data)) {
        x@meta.data$orig.ident <- x@meta.data$library
      } else {
        stop("orig.ident has multiple values but no library column found")
      }
    } else {
      stop("orig.ident has multiple values but no hash.ID column found")
    }
  }

  message("\n\n ========== import hashtag oligos (HTOs) ========== \n\n")

  # message("loading HTO matrix for sample: ", sample_name)
  sample_name <- x@meta.data$orig.ident |> as.character() |> unique()

  # get a list of all ADTs and/or HTOs
  adt_hto_features <- rownames(x[["ADT"]])
  message(
    "all antibody capture features: ",
    str_c(adt_hto_features, collapse = ", ")
  )

  # check if ADTs and HTOs are explicitly labeled
  adt_features <- adt_hto_features |> str_subset("ADT") |> unique() |> sort()
  hto_features <- adt_hto_features |> str_subset("HTO") |> unique() |> sort()

  # if all ADTs and HTOs are not explicitly labeled
  if (length(c(adt_features, hto_features)) < length(adt_hto_features)) {
    # if some HTOs are detected, assume the rest are ADTs
    if (length(hto_features) > 1) {
      adt_features <- setdiff(adt_hto_features, hto_features)
    }
    # if some ADTs are detected, assume the rest are HTOs
    if (length(adt_features) > 10) {
      hto_features <- setdiff(adt_hto_features, adt_features)
    }
    # if still no ADTs or HTOs detected, assume all are HTOs
    if (length(c(adt_features, hto_features)) == 0) {
      hto_features <- adt_hto_features
    }
  }

  message("detected ADTs: ", str_c(adt_features, collapse = ", "))
  message("detected HTOs: ", str_c(hto_features, collapse = ", "))

  # check for potential problems
  if (length(c(adt_features, hto_features)) < length(adt_hto_features)) {
    stop("missing ADTs/HTOs detected")
  }
  if (length(intersect(adt_features, hto_features)) > 0) {
    stop("conflicting ADTs/HTOs detected")
  }
  if (length(hto_features) < 2) {
    stop("no HTOs detected")
  }

  # hto_mat = Read10X(data.dir = sample_dir, gene.column = 1)
  hto_mat <- GetAssayData(x, assay = "ADT", layer = "counts")
  hto_mat <- as.matrix(hto_mat)
  hto_mat <- hto_mat[hto_features, ]
  rownames(hto_mat) <- str_remove(rownames(hto_mat), "^HTO-")

  # removed "unmapped" HTO
  # if (rownames(hto_mat)[length(rownames(hto_mat))] == "unmapped") {
  #   hto_mat = hto_mat[1:length(rownames(hto_mat))-1, ]
  # }

  scooter::write_message(
    glue("HTO library {sample_name} cells: {ncol(hto_mat)}"),
    log_file = "hto.log"
  )
  scooter::write_message(
    glue("HTO library {sample_name} HTOs: {nrow(hto_mat)}"),
    log_file = "hto.log"
  )

  # clean up hashtag matrix to match the RNA data
  # rownames(hto_mat) = str_sub(rownames(hto_mat), 1, -17)
  # hto_mat = hto_mat[sort(rownames(hto_mat)), ]
  # colnames(hto_mat) = str_c(sample_name, ":", colnames(hto_mat))

  # Seurat replaces "_" or "|" in feature names with "-"
  # rownames(hto_mat) = str_replace(rownames(hto_mat), "_", "-")

  # clean up counts matrix to match the RNA data
  common_cells <- intersect(colnames(hto_mat), colnames(x))
  if (length(common_cells) < 10) {
    stop("cell names do not match expression matrix")
  }
  hto_mat <- hto_mat[, common_cells]

  message(glue("RNA library {sample_name} cells: {ncol(x)}"))
  scooter::write_message(
    glue("RNA and HTO library {sample_name} common cells: {ncol(hto_mat)}"),
    log_file = "hto.log"
  )

  # log the detailed unfiltered HTO stats to file
  hto_col_sums <- colSums(hto_mat)
  # write(glue("unfiltered min HTO reads: {min(hto_col_sums)}"), file = "hto.log", append = TRUE)
  # write(glue("unfiltered max HTO reads: {max(hto_col_sums)}"), file = "hto.log", append = TRUE)
  # write(glue("unfiltered mean HTO reads: {round(mean(hto_col_sums), 3)}"), file = "hto.log", append = TRUE)
  # write(glue("unfiltered median HTO reads: {median(hto_col_sums)}"), file = "hto.log", append = TRUE)
  scooter::write_message(
    glue("min HTO reads: {min(hto_col_sums)}"),
    log_file = "hto.log"
  )
  scooter::write_message(
    glue("max HTO reads: {max(hto_col_sums)}"),
    log_file = "hto.log"
  )
  scooter::write_message(
    glue("mean HTO reads: {round(mean(hto_col_sums), 3)}"),
    log_file = "hto.log"
  )
  scooter::write_message(
    glue("median HTO reads: {median(hto_col_sums)}"),
    log_file = "hto.log"
  )

  # compute the deviation from the median of each cell (constant is 1.48 for normal distribution)
  # hto_col_sums = log1p(hto_col_sums)
  # hto_col_mads = (hto_col_sums - median(hto_col_sums)) / mad(hto_col_sums, constant = 1)

  # remove outlier cells based on total HTO counts
  # hto_outliers = hto_col_mads > 3
  # hto_outliers = names(hto_outliers[hto_outliers == TRUE])
  # x = subset(x, cells = hto_outliers, invert = TRUE)
  # hto_mat = hto_mat[, colnames(x)]

  # message(glue("HTO library {sample_name} filtered cells: {ncol(hto_mat)}"))
  # message(glue("HTO library {sample_name} filtered HTOs: {nrow(hto_mat)}"))

  # write(glue("HTO library {sample_name} filtered cells: {ncol(hto_mat)}"), file = "hto.log", append = TRUE)
  # write(glue("HTO library {sample_name} filtered HTOs: {nrow(hto_mat)}"), file = "hto.log", append = TRUE)

  # log the detailed filtered HTO stats to file
  # hto_col_sums = colSums(hto_mat)
  # write(glue("filtered min HTO reads: {min(hto_col_sums)}"), file = "hto.log", append = TRUE)
  # write(glue("filtered max HTO reads: {max(hto_col_sums)}"), file = "hto.log", append = TRUE)
  # write(glue("filtered mean HTO reads: {round(mean(hto_col_sums), 3)}"), file = "hto.log", append = TRUE)
  # write(glue("filtered median HTO reads: {median(hto_col_sums)}"), file = "hto.log", append = TRUE)

  # log the new RNA stats to file
  # write(glue("filtered mean num genes: {round(mean(x$detected_genes), 3)}"), file = "hto.log", append = TRUE)
  # write(glue("filtered median num genes: {median(x$detected_genes)}"), file = "hto.log", append = TRUE)

  # normalize ADT data again since some outlier cells were removed
  # x = NormalizeData(x, assay = "ADT", normalization.method = "CLR")
  # x = ScaleData(x, assay = "ADT")

  # add HTO data as a new assay independent from RNA or ADT
  x[["HTO"]] <- CreateAssayObject(counts = hto_mat)

  # normalize HTO data
  x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR")

  # perform hashtag demultiplexing
  if (method == "HTODemux") {
    # x = HTODemux(x, assay = "HTO", kfunc = "kmeans", positive.quantile = 0.999)
    x <- HTODemux(x, assay = "HTO", kfunc = "kmeans")

    # set hash.ID and reorder factors alphabetically
    x@meta.data$hash.ID <- fct_relevel(x@meta.data$hash.ID, sort)
    x@meta.data$HTODemux.ID <- x@meta.data$hash.ID

    # set hash.class and reorder factors alphabetically
    x@meta.data$hash.class <- fct_relevel(
      x@meta.data$HTO_classification.global,
      sort
    )
  } else if (method == "MULTIseqDemux") {
    # MULTIseqDemux generates only "MULTI_ID" and "MULTI_classification"
    x <- MULTIseqDemux(x, assay = "HTO", autoThresh = TRUE)

    # set hash.ID and reorder factors alphabetically
    x@meta.data$MULTI_ID <- fct_relevel(x@meta.data$MULTI_ID, sort)
    x@meta.data$hash.ID <- x@meta.data$MULTI_ID
    x@meta.data$MULTIseqDemux.ID <- x@meta.data$MULTI_ID

    # set hash.class and reorder factors alphabetically (MULTI_ID and MULTI_classification do not have a "Singlet" class)
    x@meta.data$hash.class <- as.character(x@meta.data$MULTI_ID)
    x@meta.data[
      !x@meta.data$hash.class %in% c("Doublet", "Negative"),
      "hash.class"
    ] <- "Singlet"
    x@meta.data$hash.class <- fct_relevel(x@meta.data$hash.class, sort)
  } else {
    stop("method must be HTODemux or MULTIseqDemux")
  }

  # plot HTO metrics
  plot_hto_qc(x = x, method = method)

  # save HTO stats
  x <- scooter::set_identity(x = x, identity_column = "hash.ID")
  scooter::calculate_cluster_stats(x = x, label = glue("hto.{method}"))

  # update metadata, setting the hashtag as the sample name
  x@meta.data$library <- factor(x@meta.data$orig.ident)
  x@meta.data$orig.ident <- factor(x@meta.data$hash.ID)

  x <- scooter::set_identity(x = x, identity_column = "orig.ident")

  return(x)
}

# plot HTO metrics
plot_hto_qc <- function(x, method = c("HTODemux", "MULTIseqDemux")) {
  method <- match.arg(method)

  # normalized HTO signal combined with metadata table
  hto_tbl <- GetAssayData(x, assay = "HTO") |> t()
  outlier_cells <- rownames(hto_tbl)[
    matrixStats::rowMaxs(hto_tbl) > quantile(hto_tbl, 0.99)
  ]
  hto_tbl <- hto_tbl[, sort(colnames(hto_tbl))] |> as_tibble(rownames = "cell")
  # id_tbl = x@meta.data |> as_tibble(rownames = "cell") |> select(cell, hash.ID, HTO_classification.global)
  id_tbl <- x@meta.data |>
    as_tibble(rownames = "cell") |>
    select(cell, hash.ID, hash.class)
  hto_tbl <- full_join(hto_tbl, id_tbl, by = "cell")
  hto_tbl

  # HTO color scheme
  colors_hto_names <- c(levels(hto_tbl$hash.class), levels(hto_tbl$hash.ID)) |>
    unique()
  colors_hto <- colors_clusters[1:length(colors_hto_names)]
  names(colors_hto) <- colors_hto_names

  # visualize pairs of HTO signals
  hto_facet_plot <-
    hto_tbl |>
    filter(!cell %in% outlier_cells) |>
    sample_frac() |>
    ggplot(aes(x = .panel_x, y = .panel_y, fill = hash.ID, color = hash.ID)) +
    geom_point(shape = 16, size = 0.2) +
    geom_autodensity(color = NA, fill = "gray20") +
    geom_density2d(color = "black", alpha = 0.5) +
    scale_color_manual(values = colors_hto) +
    scale_fill_manual(values = colors_hto) +
    facet_matrix(
      vars(-cell, -hash.ID, -hash.class),
      layer.diag = 2,
      layer.upper = 3
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      legend.title = element_blank(),
      strip.background = element_blank()
    )
  save_plot(
    filename = glue("qc-hto-correlation-{method}.png"),
    plot = hto_facet_plot,
    base_height = 8,
    base_width = 10
  )
  Sys.sleep(1)
  save_plot(
    filename = glue("qc-hto-correlation-{method}.pdf"),
    plot = hto_facet_plot,
    base_height = 8,
    base_width = 10
  )
  Sys.sleep(1)

  # number of counts for singlets, doublets and negative cells
  hto_counts_plot <-
    VlnPlot(
      x,
      features = "total_counts",
      group.by = "hash.class",
      pt.size = 0.1
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = colors_hto) +
    scale_y_continuous(labels = comma)
  save_plot(
    glue("qc-hto-counts-{method}.png"),
    plot = hto_counts_plot,
    base_height = 6,
    base_width = 6
  )
  Sys.sleep(1)

  # number of genes for singlets, doublets and negative cells
  hto_gene_plot <-
    VlnPlot(
      x,
      features = "detected_genes",
      group.by = "hash.class",
      pt.size = 0.1
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = colors_hto) +
    scale_y_continuous(labels = comma)
  save_plot(
    glue("qc-hto-genes-{method}.png"),
    plot = hto_gene_plot,
    base_height = 6,
    base_width = 6
  )
  Sys.sleep(1)

  # number of genes for singlets, doublets and negative cells
  hto_mito_plot <-
    VlnPlot(x, features = "pct_mito", group.by = "hash.class", pt.size = 0.1) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = colors_hto) +
    scale_y_continuous(labels = comma)
  save_plot(
    glue("qc-hto-mito-{method}.png"),
    plot = hto_mito_plot,
    base_height = 6,
    base_width = 6
  )
  Sys.sleep(1)

  group_var <- "hash.class"
  x <- scooter::set_identity(x = x, identity_column = group_var)
  plot_umap <-
    DimPlot(
      x,
      reduction = "umap",
      pt.size = get_dr_point_size(x),
      cols = colors_hto,
      na.value = "grey90",
      shuffle = TRUE,
      raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  save_plot(
    glue("dr-umap-{group_var}-{method}.png"),
    plot = plot_umap,
    base_height = 6,
    base_width = 8
  )
  Sys.sleep(1)
  # save_plot(glue("dr-umap-{group_var}-{method}.pdf"), plot = plot_umap, base_height = 6, base_width = 8)
  # Sys.sleep(1)

  group_var <- "hash.ID"
  x <- scooter::set_identity(x = x, identity_column = group_var)
  plot_umap <-
    DimPlot(
      x,
      reduction = "umap",
      pt.size = get_dr_point_size(x),
      cols = colors_hto,
      na.value = "grey90",
      shuffle = TRUE,
      raster = FALSE
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  save_plot(
    glue("dr-umap-{group_var}-{method}.png"),
    plot = plot_umap,
    base_height = 6,
    base_width = 8
  )
  Sys.sleep(1)
  # save_plot(glue("dr-umap-{group_var}-{method}.pdf"), plot = plot_umap, base_height = 6, base_width = 8)
  # Sys.sleep(1)

  if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
}

# Split a Seurat object, by sample by default. To also cluster each split,
# specify the number of dimensions for clustering.
split_seurat_object <- function(
  x,
  original_wd,
  split_var = "orig.ident",
  dr_num_dim = 0,
  dr_num_neighbors = 30
) {
  # set identity to the column used for splitting
  DefaultAssay(x) <- "RNA"
  x <- scooter::set_identity(x = x, identity_column = split_var)

  # Clean up the object: drop scale.data, reductions, and graphs.
  # Variable features are not cleared. Assay5 keeps them in the assay
  # meta.data, and scooter::normalize_counts() recomputes them downstream.
  x <- DietSeurat(
    x,
    layers = c("counts", "data"),
    dimreducs = NULL,
    graphs = NULL
  )
  x@meta.data <- x@meta.data |> select(-starts_with("snn_res."))
  x@meta.data <- x@meta.data |> select(-starts_with("res."))

  setwd(original_wd)

  # process each split group
  split_names <- Idents(x) |> as.character() |> unique() |> sort()
  for (s in split_names) {
    message(glue("split: {s}"))
    split_obj <- subset(x, idents = s)
    if (ncol(split_obj) > 10) {
      split_dir <- glue("split-{s}")
      if (dir.exists(split_dir)) {
        stop(glue("output dir {split_dir} already exists"))
      } else {
        dir.create(split_dir)
        setwd(split_dir)
        scooter::write_message(
          glue("subset {s} cells: {ncol(split_obj)}"),
          log_file = "hto.log"
        )
        # DietSeurat() above dropped scale.data/reductions, so each split needs its own variable
        # features and PCA recomputed — see scooter::normalize_counts()/run_pca()
        split_obj <- scooter::normalize_counts(
          split_obj,
          method = "log",
          assay = "RNA"
        )
        num_pcs <- min(50, ncol(split_obj) - 1, nrow(split_obj) - 1)
        split_obj <- scooter::run_pca(
          split_obj,
          num_pcs = num_pcs,
          assay = "RNA",
          var_features = TRUE
        )
        # rough first-pass tSNE/UMAP for a quick look; cluster --num_dim=<n> computes the real ones
        num_dim <- min(30, num_pcs)
        split_obj <- scooter::run_tsne(
          split_obj,
          reduction = "pca",
          num_dim = num_dim,
          assay = "RNA",
          file_format = "png"
        )
        split_obj <- scooter::run_umap(
          split_obj,
          reduction = "pca",
          num_dim = num_dim,
          assay = "RNA",
          file_format = "png"
        )
        split_obj <- scooter::set_identity(
          x = split_obj,
          identity_column = "orig.ident"
        )
        scooter::calculate_cluster_stats(split_obj, label = "sample")
        # determine clusters
        if (dr_num_dim > 0) {
          split_obj <- scooter::cluster_seurat_object(
            split_obj,
            num_dim = dr_num_dim,
            reduction = "pca",
            num_neighbors = dr_num_neighbors,
            assay = "RNA",
            metadata_file = "metadata.csv.gz",
            log_file = "hto.log"
          )
        }
        qs2::qs_save(split_obj, file = "seurat_obj.qs2", nthreads = 4)
      }
    }

    # clean up and return to the main dir before processing the next split
    if (file.exists("Rplots.pdf")) {
      file.remove("Rplots.pdf")
    }
    setwd(original_wd)
  }
}


# plot gene expression overlaid on a UMAP based on a table of genes and groups
plot_dr_umap_genes <- function(x, genes_csv) {
  # check that the input object already has UMAP computed
  if (is.null(x@reductions$umap)) {
    stop("UMAP not computed yet")
  }

  # import genes table and check that the "gene" and "group" columns exist
  if (!file.exists(genes_csv)) {
    stop(glue("genes table {genes_csv} does not exist"))
  }
  genes_tbl <- read_csv(genes_csv, col_types = cols())
  genes_tbl <- distinct(genes_tbl)
  if (!is.element("gene", colnames(genes_tbl))) {
    stop("gene table column names must include 'gene'")
  }
  if (!is.element("group", colnames(genes_tbl))) {
    stop("gene table column names must include 'group'")
  }

  # switch to RNA assay
  DefaultAssay(x) <- "RNA"
  all_genes <- rownames(x[["RNA"]])

  # use ADT assay if it has a better overlap with the given genes
  if ("ADT" %in% Assays(x)) {
    rna_matches <- intersect(unique(genes_tbl$gene), rownames(x[["RNA"]]))
    adt_matches <- intersect(unique(genes_tbl$gene), rownames(x[["ADT"]]))
    all_genes <- c(all_genes, rownames(x[["ADT"]]))
    if (length(adt_matches) > length(rna_matches)) {
      if (length(rna_matches) > 0) {
        message(
          "\n RNA assay detectable genes: ",
          str_c(sort(rna_matches), collapse = ", ")
        )
        message(
          "\n ADT assay detectable genes: ",
          str_c(sort(adt_matches), collapse = ", ")
        )
      }
      DefaultAssay(x) <- "ADT"
    }
  }

  # check for detected genes
  missing_genes <- setdiff(genes_tbl$gene, all_genes)
  genes_tbl <- genes_tbl |> filter(gene %in% all_genes)
  message("\n gene groups: ", str_c(unique(genes_tbl$group), collapse = ", "))
  message("\n detectable genes: ", str_c(genes_tbl$gene, collapse = ", "))
  message("\n missing genes: ", str_c(missing_genes, collapse = ", "))
  if (nrow(genes_tbl) == 0) {
    stop("no detectable genes")
  }
  Sys.sleep(1)

  # plot settings
  dr_point_size <- get_dr_point_size(x)
  set.seed(99)
  shuffled_cells <- sample(colnames(x))

  # plot each gene in a separate directory based on group
  for (gene_group in unique(genes_tbl$group)) {
    genes_group_tbl <- genes_tbl |> filter(group == gene_group)
    plot_dir <- glue("dr-umap-genes-{gene_group}")
    dir.create(plot_dir)
    for (g in sort(genes_group_tbl$gene)) {
      gene_umap <- scooter::plot_dr_feature(
        x,
        feature = g,
        reduction = "umap",
        cells = shuffled_cells,
        pt_size = dr_point_size
      )
      save_plot(
        glue("{plot_dir}/dr-umap-gene-{g}.png"),
        plot = gene_umap,
        base_height = 6.5,
        base_width = 7
      )
      Sys.sleep(1)
    }
  }
}


# ========== main ==========

# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# disable scientific notation
options(scipen = 999)
# default type for the bitmap devices such as png (should default to "cairo")
options(bitmapType = "cairo")

# retrieve the command-line arguments
library(docopt)
# opts = docopt(doc)
opts <- NULL
tryCatch(
  {
    opts <- docopt::docopt(doc)
  },
  error = function(err) {
    msg <- conditionMessage(err)
    if (!grepl("Usage:", msg, fixed = TRUE)) {
      message(msg)
    }
    message(doc)
    quit(save = "no", status = 1)
  }
)

# show docopt options
# print(opts)

# dependencies
load_libraries()

# set number of cores for parallel package (will use all available cores by default)
options(mc.cores = 4)
# increase the maximum amount of data that is passed to a future process from 500MB to 75GB (1GB = 1024^3)
options(future.globals.maxSize = 75 * 1024^3)
# disable future seed warning (introduced in future 1.19.0, should be fixed in Seurat 4)
options(future.rng.onMisuse = "ignore")
options(future.seed = TRUE)

# global settings
colors_samples <- scooter::get_color_scheme("samples")
colors_clusters <- scooter::get_color_scheme("clusters")

# analysis info
analysis_step <- "unknown"
out_dir <- opts$analysis_dir

# original working dir (before moving to analysis dir)
original_wd <- getwd()

# check input files convert to canonical form in case they are relative
if (!is.null(opts$sample_dir)) {
  if (dir.exists(opts$sample_dir)) {
    opts$sample_dir <- normalizePath(opts$sample_dir)
  } else {
    stop(glue("sample dir {opts$sample_dir} does not exist"))
  }
}
if (!is.null(opts$genes_csv)) {
  if (file.exists(opts$genes_csv)) {
    opts$genes_csv <- normalizePath(opts$genes_csv)
  } else {
    stop(glue("genes table {opts$genes_csv} does not exist"))
  }
}

# create a new directory for a new analysis or exit if one already exists
if (opts$create || opts$merge || opts$integrate) {
  if (opts$create) {
    analysis_step <- "create"
  }
  if (opts$merge) {
    analysis_step <- "merge"
  }
  if (opts$integrate) {
    analysis_step <- "integrate"
    seurat_obj_qs2 <- glue("{out_dir}/seurat_obj.qs2")
    seurat_obj_qs2 <- normalizePath(seurat_obj_qs2)
    integrate_num_dim <- if (is.null(opts$num_dim)) 50 else opts$num_dim
    out_dir <- glue(
      "{out_dir}-integrated-{opts$reduction}-dim{integrate_num_dim}"
    )
  }
  message(glue(
    "\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"
  ))

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
  scooter::write_message(glue("analysis: {out_dir}"), log_file = "create.log")
  scooter::write_message(
    glue("seurat version: {packageVersion('Seurat')}"),
    log_file = "create.log"
  )
  scooter::write_message(
    glue("scooter version: {packageVersion('scooter')}"),
    log_file = "create.log"
  )

  # create_seurat_object() runs the whole first-sample pipeline in one call:
  #   - read and initialize the object
  #   - attach ADT, if present
  #   - filter, normalize, and select variable features
  #   - compute PCA and a first-pass tSNE/UMAP
  #   - write every diagnostic plot for those steps
  create_args <- list(
    sample_name = opts$sample_name,
    path = opts$sample_dir,
    num_mads = opts$num_mads,
    min_genes = opts$min_genes,
    max_genes = opts$max_genes,
    min_counts = opts$min_counts,
    max_counts = opts$max_counts,
    max_mt = opts$mt,
    normalization_method = "log",
    log_file = "create.log"
  )
  # omitted rather than passed as an explicit NULL, so create_seurat_object()'s
  # own num_dim default applies instead of hand-copying it here
  if (!is.null(opts$num_dim)) {
    create_args$num_dim <- as.integer(opts$num_dim)
  }
  seurat_obj <- do.call(scooter::create_seurat_object, create_args)

  # ADT normalization has no package-level home yet: scooter::create_seurat_object() only attaches
  # the assay (via scooter::add_seurat_assay()), it does not normalize it
  if ("ADT" %in% names(seurat_obj@assays)) {
    seurat_obj <- add_adt_assay_qc(seurat_obj, sample_name = opts$sample_name)
  }

  # The counts table is rarely used, so this does not write it by default.
  # Run scooter::save_counts() by hand against the saved object if you need
  # it. save_counts() stops if the table is too big for an R matrix — max
  # around 100k x 21k.
  # scooter::save_counts(seurat_obj, file = "counts.normalized.csv.gz")

  qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)
} else if (opts$merge) {
  scooter::write_message(
    glue("seurat version: {packageVersion('Seurat')}"),
    log_file = "merge.log"
  )

  # Merge multiple samples/libraries from their previous analysis
  # directories. scooter::merge_seurat_objects() resolves each directory to
  # its seurat_obj.qs2, loads it, strips its own per-sample reductions, and
  # merges the result.
  sample_analysis_dirs <- file.path(original_wd, opts$sample_analysis_dir)
  seurat_obj <- scooter::merge_seurat_objects(
    sample_analysis_dirs,
    log_file = "merge.log"
  )

  qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)

  # the merged object has no scale.data/reductions yet, so normalize/select variable features and
  # run PCA before the diagnostic tSNE/UMAP preview
  seurat_obj <- scooter::normalize_counts(
    seurat_obj,
    method = "log",
    assay = "RNA"
  )
  num_pcs <- min(50, ncol(seurat_obj) - 1, nrow(seurat_obj) - 1)
  seurat_obj <- scooter::run_pca(
    seurat_obj,
    num_pcs = num_pcs,
    assay = "RNA",
    var_features = TRUE
  )
  num_dim <- min(50, num_pcs)
  seurat_obj <- scooter::run_tsne(
    seurat_obj,
    reduction = "pca",
    num_dim = num_dim,
    assay = "RNA",
    file_format = "png"
  )
  seurat_obj <- scooter::run_umap(
    seurat_obj,
    reduction = "pca",
    num_dim = num_dim,
    assay = "RNA",
    file_format = "png"
  )

  scooter::save_metadata(seurat_obj)

  qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)

  # save sample stats
  seurat_obj <- scooter::set_identity(
    x = seurat_obj,
    identity_column = "orig.ident"
  )
  scooter::calculate_cluster_stats(seurat_obj, label = "sample")
} else if (opts$integrate) {
  scooter::write_message(
    glue("analysis: {out_dir}"),
    log_file = "integrate.log"
  )
  scooter::write_message(
    glue("seurat version: {packageVersion('Seurat')}"),
    log_file = "integrate.log"
  )
  scooter::write_message(
    glue("scooter version: {packageVersion('scooter')}"),
    log_file = "integrate.log"
  )
  scooter::write_message(
    glue("batch variable: {opts$batch_var}"),
    log_file = "integrate.log"
  )
  scooter::write_message(
    glue("integration reduction: {opts$reduction}"),
    log_file = "integrate.log"
  )
  num_dim <- if (is.null(opts$num_dim)) 50 else as.integer(opts$num_dim)
  scooter::write_message(
    glue("dimensions: {num_dim}"),
    log_file = "integrate.log"
  )

  message("loading seurat_obj")
  seurat_obj <- scooter::resolve_seurat_object(seurat_obj_qs2)

  # run integration
  seurat_obj <- scooter::integrate_seurat_object(
    seurat_obj,
    num_dim = num_dim,
    int_reduction = opts$reduction,
    batch_var = opts$batch_var,
    log_file = "integrate.log"
  )

  qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)

  # save sample stats
  scooter::calculate_cluster_stats(seurat_obj, label = "sample")
} else {
  # all commands besides "create", "merge", and "integrate" start with an existing object
  message("loading seurat_obj")
  seurat_obj <- scooter::resolve_seurat_object("seurat_obj.qs2")

  if (opts$cluster) {
    analysis_step <- "cluster"
    message(glue(
      "\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"
    ))

    # determine clusters
    seurat_obj <- scooter::cluster_seurat_object(
      seurat_obj,
      num_dim = as.integer(opts$num_dim),
      reduction = opts$reduction,
      assay = "RNA",
      metadata_file = "metadata.csv.gz",
      log_file = "cluster.log"
    )
    qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)
  }

  if (opts$hto) {
    analysis_step <- "hto"
    message(glue(
      "\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"
    ))

    # add HTO data
    seurat_obj <- split_adt_hto_assay(seurat_obj, method = "MULTIseqDemux")
    seurat_obj@meta.data$library_hash <- str_c(
      seurat_obj@meta.data$library,
      "-",
      seurat_obj@meta.data$hash.ID
    )
    qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)

    # split singlets by HTO sample
    Idents(seurat_obj) <- "hash.class"
    seurat_obj <- subset(seurat_obj, idents = "Singlet")
    split_seurat_object(
      seurat_obj,
      original_wd = original_wd,
      split_var = "library_hash"
    )
  }

  if (opts$plot) {
    if (opts$umap) {
      analysis_step <- "plot umap"
      message(glue(
        "\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"
      ))

      if (!is.null(opts$genes_csv)) {
        plot_dr_umap_genes(x = seurat_obj, genes_csv = opts$genes_csv)
      }
      if (!is.null(opts$cluster_var)) {
        scooter::plot_dr_umap_clusters(
          x = seurat_obj,
          cluster_var = opts$cluster_var
        )
      }
    } else {
      stop("unknown plot type (only umap is currently supported)")
    }
  }

  if (opts$identify || opts$de) {
    # set resolution in the seurat object
    seurat_obj <- scooter::set_identity(
      x = seurat_obj,
      identity_column = opts$cluster_var
    )

    # use a grouping-specific sub-directory for all output
    grouping_label <- scooter::check_identity_column(
      seurat_obj,
      opts$cluster_var
    )
    grouping_label <- gsub("\\.", "", grouping_label)
    num_clusters <- Idents(seurat_obj) |> as.character() |> n_distinct()
    clust_label <- glue("clust{num_clusters}")
    res_dir <- glue("clusters-{grouping_label}-{clust_label}")
    if (!dir.exists(res_dir)) {
      dir.create(res_dir)
    }
    setwd(res_dir)

    if (opts$identify) {
      analysis_step <- "identify"
      message(glue(
        "\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"
      ))

      # create UMAP/tSNE plot (should already exist in the main directory)
      dr_filename <- glue("dr-{grouping_label}-{clust_label}")
      qs2::qs_save(seurat_obj, file = "seurat_obj.qs2", nthreads = 4)
      seurat_obj <- scooter::plot_clusters(
        seurat_obj,
        resolution = opts$cluster_var,
        filename_base = dr_filename
      )

      # cluster stat tables (number of cells and average expression)
      scooter::calculate_cluster_stats(seurat_obj, label = clust_label)
      scooter::calculate_cluster_expression(seurat_obj, label = clust_label)

      # ignore NA clusters
      seurat_obj <- subset(
        seurat_obj,
        cells = which(!is.na(Idents(seurat_obj)))
      )

      # calculate and plot standard cluster markers
      message("\n\n ========== cluster markers ========== \n\n")
      # scooter::calculate_cluster_markers(seurat_obj, label = clust_label, test = "roc")
      scooter::calculate_cluster_markers(
        seurat_obj,
        label = clust_label,
        test = "wilcox"
      )
      # scooter::calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST")

      # calculate and plot pairwise cluster markers (very slow, so skip for high number of clusters)
      num_clusters <- Idents(seurat_obj) |> as.character() |> n_distinct()
      if (num_clusters < 20) {
        message("\n\n ========== cluster markers (pairwise) ========== \n\n")
        scooter::calculate_cluster_markers(
          seurat_obj,
          label = clust_label,
          test = "wilcox",
          pairwise = TRUE
        )
        # scooter::calculate_cluster_markers(seurat_obj, label = clust_label, test = "MAST", pairwise = TRUE)
      }
    }

    # differential expression
    if (opts$de) {
      analysis_step <- "diff"
      message(glue(
        "\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"
      ))

      # calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "wilcox", group_var = opts$group_var)
      # calculate_cluster_de_genes(seurat_obj, label = clust_label, test = "MAST", group_var = opts$group_var)
      de_tbl <- differential_expression_per_cluster(
        seurat_obj,
        cluster_column = opts$cluster_var,
        group_column = opts$group_var,
        test = "wilcox"
      )
    }
  }
}

message(glue(
  "\n\n ========== finished analysis step {analysis_step} for {out_dir} ========== \n\n"
))

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

# end
