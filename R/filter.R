#' Filter cells based on the number of genes, counts, and mitochondrial reads.
#'
#' Filtering runs in two stages, and the log records how many cells each one
#' removed.
#'
#' 1. **Outlier removal.** Cells further than `num_mads` median absolute
#'    deviations from the median number of genes or counts are dropped. The
#'    median and MAD are computed on the log scale and the bounds
#'    back-transformed, which suits these right-skewed counts. `num_mads =
#'    NULL` skips it.
#' 2. **Fixed cutoffs.** `min_genes` / `max_genes` / `min_counts` /
#'    `max_counts` are applied to whatever survived stage 1, then `max_mt`.
#'    Each is skipped when `NULL`, so out of the box only `min_genes` (500) and
#'    `max_mt` (10) do anything.
#'
#' Applying the fixed cutoffs second makes them a quality floor the data-driven
#' bounds cannot undercut: a wide or contaminated distribution can put the MAD
#' bound below 500 genes, and those cells are still removed. All cutoffs are
#' inclusive.
#'
#' Genes that are no longer detected in any of the retained cells are dropped,
#' since filtering cells leaves them as all-zero rows.
#'
#' @param x A tibble with metadata (default method) or a Seurat object (Seurat
#'   method). Must have `detected_genes`/`total_counts` columns, as added by
#'   [initialize_seurat_object()].
#' @param num_mads Median absolute deviations for the outlier stage. `NULL`
#'   skips the stage.
#' @param min_genes Minimum number of genes per cell. `NULL` for no cutoff.
#' @param max_genes Maximum number of genes per cell. `NULL` for no cutoff.
#' @param min_counts Minimum number of counts per cell. `NULL` for no cutoff.
#' @param max_counts Maximum number of counts per cell. `NULL` for no cutoff.
#' @param max_mt Maximum percentage of mitochondrial reads per cell. `NULL` for
#'   no cutoff.
#' @param min_cells Minimum number of cells that must survive filtering. Below
#'   this, filtering `stop()`s rather than handing back an object too small to
#'   analyze meaningfully. `NULL` skips the check.
#' @param log_file Log file.
#'
#' @return Filtered data. The Seurat method also writes a pre-filtering QC
#'   violin plot (via [plot_metrics_distribution()]) and the unfiltered
#'   metadata (via [save_metadata()], as `metadata-unfiltered.csv.gz`), plus
#'   the post-filtering QC plots (via [plot_metrics_distribution()] and
#'   [plot_metrics_correlations()]) to the working directory.
#'
#' @export
filter_cells <- function(
  x,
  num_mads = 3,
  min_genes = 500,
  max_genes = NULL,
  min_counts = NULL,
  max_counts = NULL,
  max_mt = 10,
  min_cells = 50,
  log_file = NULL
) {
  UseMethod("filter_cells")
}

#' @export
filter_cells.default <- function(
  x,
  num_mads = 3,
  min_genes = 500,
  max_genes = NULL,
  min_counts = NULL,
  max_counts = NULL,
  max_mt = 10,
  min_cells = 50,
  log_file = NULL
) {
  # command line arguments arrive as characters, and as.numeric(NULL) is a
  # zero-length numeric, which is why everything below tests length() rather
  # than is.null()
  num_mads <- as.numeric(num_mads)
  min_genes <- as.numeric(min_genes)
  max_genes <- as.numeric(max_genes)
  min_counts <- as.numeric(min_counts)
  max_counts <- as.numeric(max_counts)
  max_mt <- as.numeric(max_mt)

  # glue() on a zero-length value returns character(0), swallowing the whole
  # message
  show <- function(cutoff) if (length(cutoff)) as.character(cutoff) else "none"

  metadata <- x |>
    as.data.frame() |>
    rownames_to_column("cell")

  # stage 1: data-driven outliers on both metrics
  num_cells <- nrow(metadata)
  if (length(num_mads)) {
    gene_bounds <- round(mad_outlier_bounds(
      metadata$detected_genes,
      nmads = num_mads,
      log = TRUE
    ))
    counts_bounds <- round(mad_outlier_bounds(
      metadata$total_counts,
      nmads = num_mads,
      log = TRUE
    ))
    metadata <- metadata |>
      filter(
        .data$detected_genes >= gene_bounds[["min"]],
        .data$detected_genes <= gene_bounds[["max"]],
        .data$total_counts >= counts_bounds[["min"]],
        .data$total_counts <= counts_bounds[["max"]]
      )
    write_message(glue("MAD outlier bounds ({num_mads} MADs)"), log_file)
    write_message(
      glue("genes: {gene_bounds[['min']]} to {gene_bounds[['max']]}"),
      log_file
    )
    write_message(
      glue("counts: {counts_bounds[['min']]} to {counts_bounds[['max']]}"),
      log_file
    )
    write_message(
      glue("cells removed as MAD outliers: {num_cells - nrow(metadata)}"),
      log_file
    )
  }

  # stage 2: fixed cutoffs, applied to whatever survived the outlier stage
  num_cells <- nrow(metadata)
  if (length(min_genes)) {
    metadata <- filter(metadata, .data$detected_genes >= min_genes)
  }
  if (length(max_genes)) {
    metadata <- filter(metadata, .data$detected_genes <= max_genes)
  }
  if (length(min_counts)) {
    metadata <- filter(metadata, .data$total_counts >= min_counts)
  }
  if (length(max_counts)) {
    metadata <- filter(metadata, .data$total_counts <= max_counts)
  }
  write_message("gene/count cutoffs", log_file)
  write_message(glue("genes: {show(min_genes)} to {show(max_genes)}"), log_file)
  write_message(
    glue("counts: {show(min_counts)} to {show(max_counts)}"),
    log_file
  )
  write_message(
    glue("cells removed by gene/count cutoffs: {num_cells - nrow(metadata)}"),
    log_file
  )

  num_cells <- nrow(metadata)
  if (length(max_mt)) {
    metadata <- filter(metadata, .data$pct_mito <= max_mt)
  }
  write_message(
    glue("max mitochondrial percentage cutoff: {show(max_mt)}"),
    log_file
  )
  write_message(
    glue(
      "cells removed by mitochondrial percentage: {num_cells - nrow(metadata)}"
    ),
    log_file
  )

  if (length(min_cells) && nrow(metadata) < min_cells) {
    stop(glue(
      "only {nrow(metadata)} cells passed filtering, below the minimum of ",
      "{min_cells} required to proceed"
    ))
  }

  pull(metadata, .data$cell)
}

#' @export
#' @importFrom Matrix rowSums
filter_cells.Seurat <- function(
  x,
  num_mads = 3,
  min_genes = 500,
  max_genes = NULL,
  min_counts = NULL,
  max_counts = NULL,
  max_mt = 10,
  min_cells = 50,
  log_file = NULL
) {
  write_message(glue("imported cells: {ncol(x)}"), log_file)
  write_message(glue("imported genes: {nrow(x)}"), log_file)
  write_message(glue("unfiltered min genes: {min(x$detected_genes)}"), log_file)
  write_message(glue("unfiltered max genes: {max(x$detected_genes)}"), log_file)
  write_message(
    glue("unfiltered mean num genes: {round(mean(x$detected_genes), 3)}"),
    log_file
  )
  write_message(
    glue("unfiltered median num genes: {median(x$detected_genes)}"),
    log_file
  )
  write_message(
    glue("unfiltered min num counts: {round(min(x$total_counts), 3)}"),
    log_file
  )
  write_message(
    glue("unfiltered max num counts: {round(max(x$total_counts), 3)}"),
    log_file
  )
  write_message(
    glue("unfiltered mean num counts: {round(mean(x$total_counts), 3)}"),
    log_file
  )
  write_message(
    glue("unfiltered median num counts: {round(median(x$total_counts), 3)}"),
    log_file
  )

  # a pre-filtering look at the same metric set the filtered plot uses, so the
  # two are comparable
  plot_metrics_distribution(
    x,
    metrics = c(
      "detected_genes",
      "total_counts",
      "pct_mito",
      "pct_ribo",
      "pct_hb"
    ),
    file = "metrics-distribution-unfiltered.png"
  )
  save_metadata(x, file = "metadata-unfiltered.csv.gz")

  # the cutoffs (including the min_cells stop()) are resolved and logged by the
  # default method, against the same metadata table
  cells_subset <- filter_cells(
    x@meta.data,
    log_file = log_file,
    num_mads = num_mads,
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    max_mt = max_mt,
    min_cells = min_cells
  )

  # subset based on the cells that passed the filtering
  x <- subset(x, cells = cells_subset)

  # genes are only filtered when the object is created, against the unfiltered
  # cells, so dropping cells strands genes that are no longer detected in any
  # of them
  detected <- Matrix::rowSums(GetAssayData(x, layer = "counts") > 0)
  keep_features <- names(detected)[detected > 0]

  # `detected` only describes the default assay, and subset(features=) drops
  # any assay that has no feature in the list, so ADT/HTO have to be carried
  # through by name or they vanish entirely. The cutoffs above are all RNA
  # metrics; they were never meant to filter the other assays.
  for (extra_assay in setdiff(names(x@assays), DefaultAssay(x))) {
    keep_features <- c(keep_features, rownames(x@assays[[extra_assay]]))
  }

  x <- subset(x, features = keep_features)

  write_message(glue("filtered cells: {ncol(x)}"), log_file)
  write_message(glue("filtered genes: {nrow(x)}"), log_file)
  write_message(glue("filtered min genes: {min(x$detected_genes)}"), log_file)
  write_message(glue("filtered max genes: {max(x$detected_genes)}"), log_file)
  write_message(
    glue("filtered mean num genes: {round(mean(x$detected_genes), 3)}"),
    log_file
  )
  write_message(
    glue("filtered median num genes: {median(x$detected_genes)}"),
    log_file
  )
  write_message(
    glue("filtered min num counts: {round(min(x$total_counts), 3)}"),
    log_file
  )
  write_message(
    glue("filtered max num counts: {round(max(x$total_counts), 3)}"),
    log_file
  )
  write_message(
    glue("filtered mean num counts: {round(mean(x$total_counts), 3)}"),
    log_file
  )
  write_message(
    glue("filtered median num counts: {round(median(x$total_counts), 3)}"),
    log_file
  )

  # filter_cells() only ever runs on a single sample (part of
  # create_seurat_object()'s pipeline; merge_seurat_objects() filters
  # per-sample before merging), so the plot can afford the full metric set
  plot_metrics_distribution(
    x,
    metrics = c(
      "detected_genes",
      "total_counts",
      "pct_mito",
      "pct_ribo",
      "pct_hb"
    ),
    file = "metrics-distribution-filtered.png"
  )
  plot_metrics_correlations(x, file = "metrics-correlations-filtered.png")

  return(x)
}

#' Outlier bounds from the median absolute deviation.
#'
#' Returns the lower and upper cutoffs at `nmads` median absolute deviations
#' from the median, the standard MAD-based outlier rule. With `log = TRUE` the
#' median and MAD are computed on the natural-log scale and the bounds are
#' back-transformed with `exp()`, which suits right-skewed counts such as the
#' number of genes or counts per cell; values must be positive in that case.
#'
#' @param values A numeric vector.
#' @param nmads Number of median absolute deviations from the median.
#' @param log Compute the median and MAD on the natural-log scale.
#'
#' @return A named numeric vector, `c(min, max)`.
#'
#' @importFrom stats median mad
#' @noRd
mad_outlier_bounds <- function(values, nmads = 3, log = FALSE) {
  x <- if (log) log(values) else values
  center <- median(x)
  spread <- mad(x)
  bounds <- c(min = center - nmads * spread, max = center + nmads * spread)
  if (log) {
    bounds <- exp(bounds)
  }
  bounds
}
