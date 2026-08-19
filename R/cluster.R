#' Reduce dimensions, cluster, and plot the result.
#'
#' The clustering step of the pipeline in one call: run tSNE and UMAP on an
#' existing "pca" reduction ([run_tsne()], [run_umap()]), build the SNN graph
#' and cluster over a range of resolutions ([calculate_clusters()]), and plot
#' each reduction by cluster.
#'
#' This function computes both reductions before clustering, since the
#' cluster plots draw on them. [run_tsne()] and [run_umap()] already write
#' their own sample-colored scatter to the working directory, so this
#' function skips plotting "by sample" again.
#'
#' The per-resolution cluster plots also write unconditionally, into a
#' `clusters-resolutions/` subdirectory. This matches every other verb's
#' standard-output convention.
#'
#' @param x Seurat object with `reduction` computed, such as the output of
#'   [create_seurat_object()] (which only ever produces "pca") or
#'   [integrate_seurat_object()] (which also produces one named after its
#'   integration method: "cca", "rpca", or "harmony").
#' @param assay Assay to attribute the tSNE/UMAP reductions to.
#'   [create_seurat_object()] creates its own assay, so it needs no such
#'   parameter. This function creates none: it only calls
#'   [run_tsne()]/[run_umap()], which read from an existing assay. So
#'   `assay` is explicit here, rather than read back off `DefaultAssay(x)`.
#' @param reduction Reduction to compute tSNE/UMAP from and to cluster on -
#'   "pca" after [create_seurat_object()], or the integration method's own name
#'   ("cca"/"rpca"/"harmony") after [integrate_seurat_object()] to cluster on
#'   the batch-corrected embedding instead of the pre-correction one. No
#'   default: which reduction is the intended one for clustering differs by
#'   which pipeline branch got you here, so every caller states it explicitly
#'   rather than risk silently picking the wrong one.
#' @param num_dim Principal components to use, bounded to 5-50.
#' @param num_neighbors Neighbors for UMAP and for the SNN graph.
#' @param metadata_file Path to write the cell metadata and embeddings to, via
#'   [save_metadata()]. `NULL` writes nothing.
#' @param log_file Filename for the log file.
#'
#' @return A Seurat object with "tsne" and "umap" reductions and one `res.<x>`
#'   column per retained clustering resolution.
#'
#' @importFrom stringr str_c str_pad str_remove str_subset
#' @export
cluster_seurat_object <- function(
  x,
  assay = "RNA",
  reduction,
  num_dim,
  num_neighbors = 30,
  metadata_file = NULL,
  log_file = NULL
) {
  num_dim <- as.integer(num_dim)
  if (num_dim < 5) {
    stop("too few dims: ", num_dim)
  }
  if (num_dim > 50) {
    stop("too many dims: ", num_dim)
  }

  cluster_colors <- get_color_scheme("clusters")

  # tSNE and UMAP are for visualization: clustering runs on the same reduction
  # below, not on assay expression data - assay is still set explicitly rather
  # than left NULL (see run_dr()'s roxygen note)
  message("\n\n ===== run_tsne() ===== \n\n")
  x <- run_tsne(
    x,
    reduction = reduction,
    num_dim = num_dim,
    assay = assay,
    file_format = "png"
  )

  message("\n\n ===== run_umap() ===== \n\n")
  x <- run_umap(
    x,
    reduction = reduction,
    num_dim = num_dim,
    num_neighbors = num_neighbors,
    assay = assay,
    file_format = "png"
  )

  message("\n\n ===== calculate_clusters() ===== \n\n")

  # multisession can be faster here on mid-size objects, but crashed (worker
  # OOM-killed) on this session's memory allocation at 50k cells - not a fixed
  # threshold, a job with more memory could push it further out. Kept
  # sequential anyway since the failure mode is unpredictable at the scale
  # merged/integrated objects tend to reach
  future::plan(future::sequential)

  x <- calculate_clusters(
    x,
    assay = assay,
    reduction = reduction,
    num_dim = num_dim,
    num_neighbors = num_neighbors,
    max_clusters = length(cluster_colors),
    log_file = log_file
  )

  # one plot per retained resolution, on both reductions
  resolution_columns <- str_subset(colnames(x@meta.data), "^res\\.")
  resolution_columns <- resolution_columns[
    order(as.numeric(str_remove(resolution_columns, "^res\\.")))
  ]

  for (resolution_column in resolution_columns) {
    resolution <- str_remove(resolution_column, "^res\\.")
    num_clusters <- x@meta.data[[resolution_column]] |>
      as.character() |>
      n_distinct()

    write_message(
      glue("resolution: {resolution} | clusters: {num_clusters}"),
      log_file
    )

    # a single cluster says nothing, and past the palette the colors repeat
    if (num_clusters < 2 || num_clusters >= length(cluster_colors)) {
      next
    }

    x <- set_identity(x = x, identity_column = resolution_column)

    # res.0.4 -> res004, so the files sort in resolution order
    resolution_tag <- format(as.numeric(resolution), nsmall = 1) |>
      str_remove("\\.") |>
      str_pad(width = 3, side = "left", pad = "0") |>
      (\(tag) str_c("res", tag))()
    file_base <- file.path(
      "clusters-resolutions",
      glue("dr-{assay}-{num_dim}-{resolution_tag}-clust{num_clusters}")
    )

    for (reduction in c("tsne", "umap")) {
      plot_dr_group(
        x,
        reduction = reduction,
        color_scheme = cluster_colors,
        na_value = "grey90",
        file_prefix = glue("{file_base}-{reduction}")
      )
    }
  }

  if (!is.null(metadata_file)) {
    save_metadata(x, file = metadata_file)
  }

  return(x)
}

#' Identify clusters of cells by graph-based clustering.
#'
#' Builds a shared nearest neighbor (SNN) graph from reduced dimensions and
#' identifies clusters over a range of resolutions.
#'
#' For a Seurat object, the cluster assignments are added to the object
#' metadata as `res.<resolution>` columns, the labels are converted to the
#' `C01`, `C02`, ... convention (1-based, zero-padded, and `C`-prefixed to
#' block downstream numeric coercion), and resolutions that do not yield any
#' new clusters are dropped.
#'
#' @param x A matrix of cell embeddings (cells as rows) or a Seurat object.
#' @param num_dim Number of dimensions to use.
#' @param num_neighbors Number of neighbors (`k.param`) used to build the SNN
#'   graph.
#' @param res Clustering resolution(s). If `NULL`, a range of resolutions is
#'   used.
#' @param algorithm Clustering algorithm: 1 = Louvain, 2 = Louvain with
#'   multilevel refinement, 3 = SLM, 4 = Leiden (the default). See
#'   `Seurat::FindClusters()`.
#' @param log_file Filename for the log file.
#' @param ... Arguments passed to the individual methods.
#'
#' @return Cluster assignments as a data frame (default method) or the Seurat
#'   object with the cluster assignments added to its metadata (Seurat method).
#'
#' @export
calculate_clusters <- function(x, ...) {
  UseMethod("calculate_clusters")
}

#' @importFrom rlang check_dots_used
#' @export
calculate_clusters.default <- function(
  x,
  num_dim,
  num_neighbors = 20,
  res = NULL,
  algorithm = 4,
  log_file = NULL,
  ...
) {
  # TODO: allow UMAP graphs to be used

  if (num_dim < 5) {
    stop("too few dims: ", num_dim)
  }
  if (num_dim > 50) {
    stop("too many dims: ", num_dim)
  }
  if (ncol(x) < num_dim) {
    stop(glue("not enough dimensions available: {ncol(x)} < {num_dim}"))
  }

  message("\n\n ===== Seurat::FindNeighbors() ===== \n\n")

  # construct a Shared Nearest Neighbor (SNN) graph from the reduced dimensions
  # FindNeighbors()'s annoy search (the nn.method default) forces itself to
  # sequential internally and restores the caller's plan on exit - true since
  # Seurat 3.2.2 (2020-09-25)
  snn_graph <- FindNeighbors(
    x[, 1:num_dim],
    k.param = num_neighbors,
    compute.SNN = TRUE,
    verbose = FALSE
  )$snn

  message("\n\n ===== Seurat::FindClusters() ===== \n\n")

  # resolutions for graph-based clustering increased resolution values lead to
  # more clusters (recommendation: 0.6-1.2 for 3K cells, 2-4 for 33K cells)
  if (is.null(res)) {
    res <- seq(0.1, 2.0, 0.1)
    if (nrow(x) > 10000) res <- c(res, 3, 5, 7, 9)
  }

  # leiden_method "igraph" avoids casting large data to a dense matrix.
  # Unlike the "leidenbase" default, it only needs igraph, which Seurat
  # already imports.
  #
  # The `method` argument that used to select this was deprecated in Seurat
  # 5.2.0. The "igraph" option arrived in Seurat 5.3.1, hence the
  # DESCRIPTION minimum.
  #
  # random.seed must be > 0 for Leiden. Otherwise FindClusters() resets it
  # and warns.
  #
  # The default method names the resulting columns "res.<resolution>".
  clusters <- FindClusters(
    snn_graph,
    algorithm = algorithm,
    leiden_method = "igraph",
    resolution = res,
    random.seed = 1,
    verbose = FALSE
  )

  # this method takes no extra params at all - see the matching note in
  # run_pca.default()
  check_dots_used()

  clusters
}

#' @param assay Assay logged as the one this clustering run is over (Seurat
#'   method). Purely for the log message - clustering itself runs on
#'   `Embeddings(x, reduction = reduction)`, not on assay data, so nothing here
#'   reads it back off the object. No default: `DefaultAssay(x)` is mutable
#'   state that need not match the assay `reduction` actually came from, so the
#'   caller states it explicitly rather than have the log silently report the
#'   wrong one.
#' @param reduction Reduction to take the cell embeddings from (Seurat method).
#'   No default - always state it explicitly.
#' @param max_clusters Drop resolutions yielding this many clusters or more.
#'   Defaults to the length of the cluster color scheme, since resolutions that
#'   exceed it cannot be plotted (Seurat method).
#'
#' @rdname calculate_clusters
#' @importFrom rlang check_dots_used
#' @export
calculate_clusters.Seurat <- function(
  x,
  assay,
  reduction,
  num_dim,
  num_neighbors = 20,
  res = NULL,
  algorithm = 4,
  log_file = NULL,
  max_clusters = length(get_color_scheme("clusters")),
  ...
) {
  message_str <- glue(
    "assay: {assay}
                      reduction: {reduction}
                      num dims: {num_dim}"
  )
  write_message(message_str, log_file)

  # cluster on the reduced dimensions, not on the expression values
  clusters <- calculate_clusters(
    x = Embeddings(x, reduction = reduction),
    num_dim = num_dim,
    num_neighbors = num_neighbors,
    res = res,
    algorithm = algorithm,
    log_file = log_file
  )

  x <- AddMetaData(x, metadata = clusters)

  # for the calculated resolutions: relabel clusters and drop the redundant
  # resolutions (sort by resolution value, since sorting the column names is
  # only correct by accident)
  res_cols <- colnames(clusters)
  res_cols <- res_cols[order(as.numeric(sub("^res\\.", "", res_cols)))]
  num_clusters_prev <- 1

  for (res_col in res_cols) {
    res_vector <- as.character(x@meta.data[, res_col])
    num_clusters_cur <- length(unique(res_vector))

    if (
      num_clusters_cur > num_clusters_prev && num_clusters_cur < max_clusters
    ) {
      # Relabel identities so they start at 1, are zero-padded to avoid
      # sorting issues, and are "C"-prefixed to avoid downstream numeric
      # conversions.
      #
      # SLM returns 0-based labels; Leiden returns 1-based. Shift from the
      # observed minimum instead of assuming either one. Already relabelled
      # columns stay untouched.
      if (all(grepl("^[0-9]+$", res_vector))) {
        res_num <- as.numeric(res_vector)
        x@meta.data[, res_col] <- factor(sprintf(
          "C%02d",
          res_num - min(res_num) + 1
        ))
      }
    } else {
      # remove the resolution if it has the same number of clusters as the
      # previous one
      x@meta.data[, res_col] <- NULL
    }

    # update the cluster count for the next iteration, whether or not the
    # column was kept
    num_clusters_prev <- num_clusters_cur
  }

  kept_cols <- intersect(res_cols, colnames(x@meta.data))
  message_str <- glue(
    "cluster resolutions: {paste(kept_cols, collapse = ', ')}"
  )
  write_message(message_str, log_file)

  # this method takes no extra params at all - see the matching note in
  # run_pca.default()
  check_dots_used()

  return(x)
}

#' Per-cluster cell counts and a joined metadata+embeddings table.
#'
#' Writes `metadata-<label>.csv.gz` (cell metadata joined with the tSNE/UMAP
#' embeddings, skipped if neither reduction is present) and
#' `summary-<label>.csv` (cell counts and percentages per cluster, split
#' further by sample into `summary-<label>-per-sample.csv` when that differs).
#'
#' @param x Seurat object with an identity set, such as via [set_identity()].
#' @param label Suffix for the output filenames.
#'
#' @return `x`, invisibly. Called for its file output.
#'
#' @importFrom readr write_csv
#' @importFrom stringr str_c
#' @export
calculate_cluster_stats <- function(x, label) {
  message("cluster names: ", str_c(levels(x), collapse = ", "))

  # compile relevant cell metadata into a single table
  x$cluster <- Idents(x)
  metadata_tbl <- as_tibble(x@meta.data, rownames = "cell")
  metadata_tbl <- select(
    metadata_tbl,
    cell,
    total_counts,
    detected_genes,
    pct_mito,
    sample_name = orig.ident,
    cluster
  )
  if (!is.null(x@reductions$tsne)) {
    tsne_tbl <- x[["tsne"]]@cell.embeddings |>
      round(3) |>
      as.data.frame() |>
      rownames_to_column("cell")
    umap_tbl <- x[["umap"]]@cell.embeddings |>
      round(3) |>
      as.data.frame() |>
      rownames_to_column("cell")
    metadata_tbl <- metadata_tbl |>
      full_join(tsne_tbl, by = "cell") |>
      full_join(umap_tbl, by = "cell")
    metadata_tbl <- metadata_tbl |> arrange(cell)
    write_csv(metadata_tbl, glue("metadata-{label}.csv.gz"))
  }
  Sys.sleep(1)

  # get number of cells split by cluster and by sample (orig.ident)
  summary_cluster_sample <-
    metadata_tbl |>
    select(cluster, sample_name) |>
    mutate(num_cells_total = n()) |>
    group_by(sample_name) |>
    mutate(num_cells_sample = n()) |>
    group_by(cluster) |>
    mutate(num_cells_cluster = n()) |>
    group_by(cluster, sample_name) |>
    mutate(num_cells_cluster_sample = n()) |>
    ungroup() |>
    distinct() |>
    mutate(
      pct_cells_cluster = num_cells_cluster / num_cells_total,
      pct_cells_cluster_sample = num_cells_cluster_sample / num_cells_sample
    ) |>
    mutate(
      pct_cells_cluster = round(pct_cells_cluster * 100, 1),
      pct_cells_cluster_sample = round(pct_cells_cluster_sample * 100, 1)
    ) |>
    arrange(cluster, sample_name)

  # get number of cells split by cluster (ignore samples)
  summary_cluster <- summary_cluster_sample |>
    select(-contains("sample")) |>
    distinct()
  write_csv(summary_cluster, glue("summary-{label}.csv"))
  Sys.sleep(1)

  # export results split by sample if multiple samples per cluster are present
  if (nrow(summary_cluster_sample) > nrow(summary_cluster)) {
    write_csv(summary_cluster_sample, glue("summary-{label}-per-sample.csv"))
    Sys.sleep(1)
  }

  invisible(x)
}

#' Per-cluster average expression (non-log space).
#'
#' Writes `expression-mean-<label>.csv`, and additionally
#' `expression-mean-<label>-per-sample.csv` when `x` has more than one sample.
#'
#' @param x Seurat object with an identity set, such as via [set_identity()].
#' @param label Suffix for the output filenames.
#'
#' @return `x`, invisibly. Called for its file output.
#'
#' @importFrom readr write_csv
#' @export
calculate_cluster_expression <- function(x, label) {
  message("cluster names: ", str_c(levels(x), collapse = ", "))

  # gene expression for an "average" cell in each identity class (averaging and
  # output are in non-log space)
  cluster_avg_exp <- AverageExpression(x, assay = "RNA", verbose = FALSE)[[
    "RNA"
  ]]
  cluster_avg_exp <- cluster_avg_exp |>
    round(3) |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    arrange(gene)
  write_csv(cluster_avg_exp, glue("expression-mean-{label}.csv"))
  Sys.sleep(1)

  # cluster averages split by sample (orig.ident)
  num_samples <- n_distinct(x@meta.data$orig.ident)
  if (num_samples > 1) {
    x@meta.data$persample <- paste0(Idents(x), "|", x@meta.data$orig.ident)
    sample_avg_exp <- AverageExpression(
      x,
      assay = "RNA",
      group.by = "persample",
      verbose = FALSE
    )[["RNA"]]
    sample_avg_exp <- sample_avg_exp |>
      round(3) |>
      as.data.frame() |>
      rownames_to_column("gene") |>
      arrange(gene)
    write_csv(sample_avg_exp, glue("expression-mean-{label}-per-sample.csv"))
    Sys.sleep(1)
  }

  invisible(x)
}

#' Calculate cluster markers (versus all other clusters, or pairwise) and plot
#' them.
#'
#' `test`, forwarded to [Seurat::FindAllMarkers()]/[Seurat::FindMarkers()]:
#' `"wilcox"` (Wilcoxon rank sum, the default), `"roc"` (classification power
#' from 0 = random to 1 = perfect), `"bimod"` (McDavid et al. 2013), `"tobit"`
#' (Trapnell et al. 2014), or `"MAST"` (Finak et al. 2015).
#'
#' `pairwise = TRUE` compares each cluster against each other cluster
#' individually (rather than each cluster against all others combined), keeping
#' only genes significant in every comparison - local and global markers both,
#' at the cost of being much slower for many clusters.
#'
#' @param x Seurat object with an identity set.
#' @param label Suffix for the output filenames.
#' @param test Statistical test - see above.
#' @param pairwise Compare clusters pairwise instead of each cluster against all
#'   others - see above.
#'
#' @return `x`, invisibly. Called for its file and plot output.
#'
#' @importFrom readr write_csv
#' @export
calculate_cluster_markers <- function(x, label, test, pairwise = FALSE) {
  message("cluster set: ", label)
  message("marker test: ", test)

  # get cluster names, keeping only clusters with more than 10 cells
  clusters <- Idents(x) |> as.character() |> unique() |> sort()
  clusters <- clusters[table(Idents(x)) > 10]

  if (!pairwise) {
    # standard cluster markers calculation
    markers_dir <- "markers-global"

    # Capture output to avoid excessive warnings.
    #
    # logfc.threshold = log2(1.1) = 0.137, below Seurat's default of 0.25
    # (log2(1.189)).
    markers_log <-
      capture.output(
        {
          all_markers <-
            FindAllMarkers(
              x,
              assay = "RNA",
              test.use = test,
              min.pct = 0.1,
              logfc.threshold = log2(1.1),
              base = 2,
              fc.name = "log2FC",
              only.pos = TRUE,
              min.diff.pct = -Inf,
              verbose = FALSE
            )
        },
        type = "message"
      )

    # do some light filtering and clean up (ROC test returns slightly different
    # output)
    if (test == "roc") {
      all_markers <-
        all_markers |>
        select(cluster, gene, log2FC = avg_diff, myAUC, power) |>
        filter(power > 0.3) |>
        mutate(
          log2FC = round(log2FC, 5),
          myAUC = round(myAUC, 5),
          power = round(power, 5)
        ) |>
        arrange(cluster, -power)
      top_markers <- all_markers |> filter(log2FC > 0)
      top_markers <- top_markers |>
        group_by(cluster) |>
        top_n(50, power) |>
        ungroup()
    } else {
      all_markers <-
        all_markers |>
        select(cluster, gene, log2FC, p_val, p_val_adj) |>
        filter(p_val_adj < 0.001) |>
        mutate(log2FC = round(log2FC, 5)) |>
        arrange(cluster, p_val_adj, p_val)
      top_markers <- all_markers |> filter(log2FC > 0)
      top_markers <- top_markers |>
        group_by(cluster) |>
        top_n(50, log2FC) |>
        ungroup()
    }
  } else {
    # pairwise (each cluster versus each other cluster) cluster markers
    # calculation
    markers_dir <- "markers-pairwise"

    unfiltered_markers <- tibble(
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
        # Find differentially expressed genes between these two clusters.
        # Use a low fold-change cutoff, to maximize the chance each gene
        # appears in every pairwise comparison. Capture output to avoid
        # excessive warnings.
        markers_log <-
          capture.output(
            {
              cur_markers <-
                FindMarkers(
                  x,
                  assay = "RNA",
                  ident.1 = cluster1,
                  ident.2 = cluster2,
                  test.use = test,
                  min.pct = 0.1,
                  logfc.threshold = 0,
                  base = 2,
                  fc.name = "log2FC",
                  only.pos = TRUE,
                  min.diff.pct = -Inf,
                  verbose = FALSE
                )
            },
            type = "message"
          )

        # clean up markers table (would need to be modified for "roc" test)
        cur_markers <-
          cur_markers |>
          rownames_to_column("gene") |>
          mutate(cluster = cluster1) |>
          mutate(cluster2 = cluster2) |>
          filter(p_val_adj < 0.01) |>
          mutate(log2FC = round(log2FC, 5)) |>
          select(one_of(colnames(unfiltered_markers)))

        unfiltered_markers <- bind_rows(unfiltered_markers, cur_markers)
      }
    }

    # adjust test name for output
    test <- glue("pairwise-{test}")

    # sort the markers to make the table more readable
    unfiltered_markers <-
      unfiltered_markers |>
      distinct() |>
      add_count(cluster, gene) |>
      rename(cluster_gene_n = n) |>
      arrange(cluster, gene, cluster2)

    # filter for genes that are significant compared to all other clusters
    all_markers <-
      unfiltered_markers |>
      filter(cluster_gene_n == (length(clusters) - 1)) |>
      select(-cluster_gene_n)

    # extract the lowest and highest fold changes and p-values
    all_markers <-
      all_markers |>
      group_by(cluster, gene) |>
      summarize_at(
        c("log2FC", "p_val", "p_val_adj"),
        list(min = min, max = max)
      ) |>
      ungroup() |>
      arrange(cluster, -log2FC_min)

    top_markers <- all_markers |>
      group_by(cluster) |>
      top_n(50, log2FC_min) |>
      ungroup()
  }

  if (!dir.exists(markers_dir)) {
    dir.create(markers_dir)
  }
  filename_base <- glue("{markers_dir}/markers-{label}-{test}")

  # save unfiltered markers for pairwise comparisons
  if (pairwise) {
    unfiltered_markers_csv <- glue("{filename_base}-unfiltered.csv")
    message("unfiltered markers: ", unfiltered_markers_csv)
    write_csv(unfiltered_markers, unfiltered_markers_csv)
    Sys.sleep(1)
  }

  all_markers_csv <- glue("{filename_base}-all.csv")
  message("all markers: ", all_markers_csv)
  write_csv(all_markers, all_markers_csv)
  Sys.sleep(1)

  top_markers_csv <- glue("{filename_base}-top.csv")
  message("top markers: ", top_markers_csv)
  write_csv(top_markers, top_markers_csv)
  Sys.sleep(1)

  plot_cluster_markers_heatmap(
    x,
    markers_tbl = all_markers,
    num_genes = c(5, 10, 20),
    filename_base = filename_base
  )

  # plot top cluster markers for each cluster
  for (cluster_name in clusters) {
    filename_cluster_base <- glue(
      "{markers_dir}/markers-{label}-{cluster_name}-{test}"
    )
    cluster_markers <- top_markers |> filter(cluster == cluster_name)
    if (nrow(cluster_markers) > 9) {
      Sys.sleep(1)
      top_cluster_markers <- cluster_markers |> head(12) |> pull(gene)
      plot_cluster_markers_top(
        x,
        genes = top_cluster_markers,
        filename_base = filename_cluster_base
      )
    }
  }

  invisible(x)
}
