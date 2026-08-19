# Split a Seurat object into per-batch layers and prepare it for integration.

Drops the scaled data, reductions, graphs, and clustering columns, then
splits the assay layers by a metadata column — splitting the layers
rather than the object is what the Seurat 5 `IntegrateLayers()` workflow
expects.

## Usage

``` r
split_layers_by_batch(
  x,
  batch_var = "orig.ident",
  assay = "RNA",
  num_variable_genes = 3000,
  num_dim = 50,
  log_file = NULL
)
```

## Arguments

- x:

  Seurat object.

- batch_var:

  Metadata column to split the layers by.

- assay:

  Assay to split.

- num_variable_genes:

  Number of variable features. Matches the single-sample default.

- num_dim:

  Number of principal components to compute.

- log_file:

  Filename for the log file.

## Value

A Seurat object with split layers and a "pca" reduction.

## Details

Everything after the split runs through the same verbs as the
single-sample path
([`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md),
[`run_pca()`](https://igordot.github.io/scooter/reference/run_pca.md)),
so the two workflows cannot drift apart on parameters.
[`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md)
then consumes the resulting "pca" reduction.
