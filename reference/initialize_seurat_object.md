# Create a new Seurat object from a matrix.

For an RNA assay the per-cell QC metrics are added at the same time:
`pct_mito`/`pct_ribo`/ `pct_hb` (via
[`add_gene_class_percent()`](https://igordot.github.io/scooter/reference/add_gene_class_percent.md)),
plus `detected_genes` and `total_counts` as clearer aliases of
`nFeature_RNA`/`nCount_RNA`, and `sample_name` as a clearer alias of
`orig.ident`. All six are metadata columns as soon as the object exists,
so
[`plot_metrics_distribution()`](https://igordot.github.io/scooter/reference/plot_metrics_distribution.md)
(which takes its `metrics` explicitly, no default) and
[`plot_metrics_correlations()`](https://igordot.github.io/scooter/reference/plot_metrics_correlations.md)
can use any of them without further preparation.

## Usage

``` r
initialize_seurat_object(
  counts_matrix,
  assay = "RNA",
  min_cells = 1,
  min_genes = 1,
  project = "proj",
  log_file = NULL
)
```

## Arguments

- counts_matrix:

  A matrix of raw counts.

- assay:

  Seurat assay to add the data to.

- min_cells:

  Include genes/features detected in at least this many cells.

- min_genes:

  Include cells where at least this many genes/features are detected.

- project:

  Project name for Seurat object.

- log_file:

  Filename for the logfile.

## Value

Seurat object.
