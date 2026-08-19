# Save a Seurat object's counts matrix to a csv file.

Deprecated: renamed to
[`save_counts()`](https://igordot.github.io/scooter/reference/save_counts.md),
which also renamed `slot` to `layer` (Seurat 5 terminology) and replaced
`out_dir`/`proj_name`/`label` with an explicit `file` path.

## Usage

``` r
save_seurat_counts_matrix(
  seurat_obj,
  proj_name = "",
  label = "",
  out_dir = ".",
  assay = "RNA",
  slot = "data",
  log_file = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- proj_name:

  Name of the project that will be the prefix of the file name.

- label:

  An optional label for the file.

- out_dir:

  Directory in which to save csv.

- assay:

  The assay within the Seurat object to retrieve data from.

- slot:

  The layer within the Seurat object to retrieve data from.

- log_file:

  Unused. Kept for signature compatibility.

## Value

A csv file in `out_dir`.
