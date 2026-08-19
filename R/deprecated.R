# Deprecated functions, kept as thin wrappers so old code keeps working (with a
# warning) after the Seurat 5 conversion renamed them. Each just forwards to
# its replacement - no logic lives here. Do not add new functionality to this
# file; fix the replacement instead.

#' Create a new Seurat object from a matrix.
#'
#' @description
#' Deprecated: renamed to [initialize_seurat_object()] - not
#' `create_seurat_object()`, despite the similar name. That function was later
#' repurposed into the full read/initialize/filter/normalize/PCA pipeline and
#' takes `sample_name`/`path`, not `counts_matrix`, so it errors on an
#' old-style call.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param assay Seurat assay to add the data to.
#' @param min_cells Include genes/features detected in at least this many cells.
#' @param min_genes Include cells where at least this many genes/features are
#'   detected.
#' @param log_file Filename for the logfile.
#' @param project Project name for Seurat object.
#'
#' @return Seurat object.
#'
#' @keywords internal
#' @export
create_seurat_obj <- function(
  counts_matrix,
  assay = "RNA",
  min_cells = 1,
  min_genes = 1,
  log_file = NULL,
  project = "proj"
) {
  lifecycle::deprecate_warn(
    "0.0.0.9005",
    "create_seurat_obj()",
    "initialize_seurat_object()"
  )
  initialize_seurat_object(
    counts_matrix = counts_matrix,
    assay = assay,
    min_cells = min_cells,
    min_genes = min_genes,
    project = project,
    log_file = log_file
  )
}

#' Read in Gene Expression and Antibody Capture data from a 10x Genomics Cell
#' Ranger sparse matrix or from a text file.
#'
#' @description
#' Deprecated: renamed to [read_counts_file()].
#'
#' @param sample_name A character that will be used as a prefix for all cell
#'   names.
#' @param path Path to directory containing 10x matrix, or path to a text file.
#' @param log_file Filename for the log file.
#'
#' @return Named list of matrices. One matrix for each data type.
#'
#' @keywords internal
#' @export
load_sample_counts_matrix <- function(sample_name, path, log_file = NULL) {
  lifecycle::deprecate_warn(
    "0.0.0.9005",
    "load_sample_counts_matrix()",
    "read_counts_file()"
  )
  read_counts_file(sample_name = sample_name, path = path, log_file = log_file)
}

#' Save a Seurat object's counts matrix to a csv file.
#'
#' @description
#' Deprecated: renamed to [save_counts()], which also renamed `slot` to `layer`
#' (Seurat 5 terminology) and replaced `out_dir`/`proj_name`/`label` with an
#' explicit `file` path.
#'
#' @param seurat_obj A Seurat object.
#' @param proj_name Name of the project that will be the prefix of the file
#'   name.
#' @param label An optional label for the file.
#' @param out_dir Directory in which to save csv.
#' @param assay The assay within the Seurat object to retrieve data from.
#' @param slot The layer within the Seurat object to retrieve data from.
#' @param log_file Unused. Kept for signature compatibility.
#'
#' @return A csv file in `out_dir`.
#'
#' @keywords internal
#' @export
save_seurat_counts_matrix <- function(
  seurat_obj,
  proj_name = "",
  label = "",
  out_dir = ".",
  assay = "RNA",
  slot = "data",
  log_file = NULL
) {
  lifecycle::deprecate_warn(
    "0.0.0.9005",
    "save_seurat_counts_matrix()",
    "save_counts()"
  )
  save_counts(
    x = seurat_obj,
    assay = assay,
    layer = slot,
    file = glue("{out_dir}/{proj_name}.{label}.counts.csv.gz")
  )
}

#' Calculate mitochondrial percentage from Seurat object.
#'
#' @description
#' Deprecated: renamed to [add_gene_class_percent()], which also adds
#' `pct_ribo` and `pct_hb`.
#'
#' @param x A Seurat object.
#'
#' @return Seurat object.
#'
#' @keywords internal
#' @export
calculate_mito_pct <- function(x) {
  lifecycle::deprecate_warn(
    "0.0.0.9005",
    "calculate_mito_pct()",
    "add_gene_class_percent()"
  )
  add_gene_class_percent(x)
}

#' Filter cells based on the number of genes, counts, and mitochondrial reads.
#'
#' @description
#' Deprecated: renamed to [filter_cells()] to read as "filter" (verb) + "cells"
#' (object) rather than the vague "data", and to avoid leading with the bare
#' word "filter", which collides conceptually with
#' `dplyr::filter()`/`stats::filter()`.
#'
#' @param x A tibble with metadata (default method) or a Seurat object (Seurat
#'   method).
#' @param num_mads Median absolute deviations for the outlier stage. `NULL`
#'   skips the stage.
#' @param min_genes Minimum number of genes per cell. `NULL` for no cutoff.
#' @param max_genes Maximum number of genes per cell. `NULL` for no cutoff.
#' @param min_counts Minimum number of counts per cell. `NULL` for no cutoff.
#' @param max_counts Maximum number of counts per cell. `NULL` for no cutoff.
#' @param max_mt Maximum percentage of mitochondrial reads per cell. `NULL` for
#'   no cutoff.
#' @param min_cells Minimum number of cells that must survive filtering. `NULL`
#'   skips the check.
#' @param log_file Log file.
#'
#' @return Filtered data. See [filter_cells()].
#'
#' @keywords internal
#' @export
filter_data <- function(
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
  lifecycle::deprecate_warn("0.0.0.9005", "filter_data()", "filter_cells()")
  filter_cells(
    x = x,
    num_mads = num_mads,
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    max_mt = max_mt,
    min_cells = min_cells,
    log_file = log_file
  )
}

#' Normalize the counts, select variable features, and scale.
#'
#' @description
#' Deprecated: renamed to [normalize_counts()] - "data" didn't say what was
#' being normalized; "counts" matches the vocabulary already used elsewhere in
#' the package (`save_counts()`, `total_counts`, `read_counts_file()`).
#'
#' @param x A Seurat object.
#' @param method "log" for `NormalizeData()` + `FindVariableFeatures()` +
#'   `ScaleData()` (log-normalized to a scale factor of 10,000), or "sct" for
#'   `SCTransform()`.
#' @param num_variable_genes Number of variable features.
#' @param assay Assay to normalize.
#' @param log_file Log file.
#'
#' @return The processed Seurat object. See [normalize_counts()].
#'
#' @keywords internal
#' @export
normalize_data <- function(
  x,
  method = "log",
  num_variable_genes = 3000,
  assay = "RNA",
  log_file = NULL
) {
  lifecycle::deprecate_warn(
    "0.0.0.9005",
    "normalize_data()",
    "normalize_counts()"
  )
  normalize_counts(
    x = x,
    method = method,
    num_variable_genes = num_variable_genes,
    assay = assay,
    log_file = log_file
  )
}
