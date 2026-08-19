# Save the cell metadata and the reduction embeddings as a single table.

Joins the requested reductions onto the cell metadata on a `cell`
column. Reductions that are not present in the object are skipped, so
this works before and after the dimensionality reduction steps have run.

## Usage

``` r
save_metadata(
  x,
  reduction = c("tsne", "umap"),
  digits = 3,
  file = "metadata.csv.gz"
)
```

## Arguments

- x:

  A Seurat object.

- reduction:

  Reductions to join onto the metadata.

- digits:

  Number of digits to round the embeddings to.

- file:

  Path to write the table to as a csv. Set to `NULL` to skip writing and
  just return the tibble.

## Value

A tibble with one row per cell.
