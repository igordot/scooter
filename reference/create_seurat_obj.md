# Create a new Seurat object from a matrix.

Deprecated: renamed to
[`initialize_seurat_object()`](https://igordot.github.io/scooter/reference/initialize_seurat_object.md) -
not
[`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md),
despite the similar name. That function was later repurposed into the
full read/initialize/filter/normalize/PCA pipeline and takes
`sample_name`/`path`, not `counts_matrix`, so it errors on an old-style
call.

## Usage

``` r
create_seurat_obj(
  counts_matrix,
  assay = "RNA",
  min_cells = 1,
  min_genes = 1,
  log_file = NULL,
  project = "proj"
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

- log_file:

  Filename for the logfile.

- project:

  Project name for Seurat object.

## Value

Seurat object.
