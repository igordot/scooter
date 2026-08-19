# Filter cells based on the number of genes, counts, and mitochondrial reads.

Filtering runs in two stages, and the log records how many cells each
one removed.

## Usage

``` r
filter_cells(
  x,
  num_mads = 3,
  min_genes = 500,
  max_genes = NULL,
  min_counts = NULL,
  max_counts = NULL,
  max_mt = 10,
  min_cells = 50,
  log_file = NULL
)
```

## Arguments

- x:

  A tibble with metadata (default method) or a Seurat object (Seurat
  method). Must have `detected_genes`/`total_counts` columns, as added
  by
  [`initialize_seurat_object()`](https://igordot.github.io/scooter/reference/initialize_seurat_object.md).

- num_mads:

  Median absolute deviations for the outlier stage. `NULL` skips the
  stage.

- min_genes:

  Minimum number of genes per cell. `NULL` for no cutoff.

- max_genes:

  Maximum number of genes per cell. `NULL` for no cutoff.

- min_counts:

  Minimum number of counts per cell. `NULL` for no cutoff.

- max_counts:

  Maximum number of counts per cell. `NULL` for no cutoff.

- max_mt:

  Maximum percentage of mitochondrial reads per cell. `NULL` for no
  cutoff.

- min_cells:

  Minimum number of cells that must survive filtering. Below this,
  filtering [`stop()`](https://rdrr.io/r/base/stop.html)s rather than
  handing back an object too small to analyze meaningfully. `NULL` skips
  the check.

- log_file:

  Log file.

## Value

Filtered data. The Seurat method also writes a pre-filtering QC violin
plot (via
[`plot_metrics_distribution()`](https://igordot.github.io/scooter/reference/plot_metrics_distribution.md))
and the unfiltered metadata (via
[`save_metadata()`](https://igordot.github.io/scooter/reference/save_metadata.md),
as `metadata-unfiltered.csv.gz`), plus the post-filtering QC plots (via
[`plot_metrics_distribution()`](https://igordot.github.io/scooter/reference/plot_metrics_distribution.md)
and
[`plot_metrics_correlations()`](https://igordot.github.io/scooter/reference/plot_metrics_correlations.md))
to the working directory.

## Details

1.  **Outlier removal.** Cells further than `num_mads` median absolute
    deviations from the median number of genes or counts are dropped.
    The median and MAD are computed on the log scale and the bounds
    back-transformed, which suits these right-skewed counts.
    `num_mads = NULL` skips it.

2.  **Fixed cutoffs.** `min_genes` / `max_genes` / `min_counts` /
    `max_counts` are applied to whatever survived stage 1, then
    `max_mt`. Each is skipped when `NULL`, so out of the box only
    `min_genes` (500) and `max_mt` (10) do anything.

Applying the fixed cutoffs second makes them a quality floor the
data-driven bounds cannot undercut: a wide or contaminated distribution
can put the MAD bound below 500 genes, and those cells are still
removed. All cutoffs are inclusive.

Genes that are no longer detected in any of the retained cells are
dropped, since filtering cells leaves them as all-zero rows.
