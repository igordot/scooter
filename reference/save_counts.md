# Save the counts matrix as a single table.

Layers left split by a [`merge()`](https://rdrr.io/r/base/merge.html)
are rejoined first, since `LayerData()` would otherwise return only the
first of them, and warn rather than error while doing it.

## Usage

``` r
save_counts(x, assay = "RNA", layer = "data", digits = 3, file)
```

## Arguments

- x:

  A Seurat object.

- assay:

  Assay to retrieve the counts from.

- layer:

  Layer to retrieve. The default "data" layer holds the normalized and
  log-transformed expression used for visualization and most
  differential expression tests.

- digits:

  Number of digits to round the values to. Ignored for a layer that
  already holds whole numbers, such as raw counts.

- file:

  Path to write the table to as a csv, or `NULL` to skip writing and
  just return the tibble. Has no default: every caller states where the
  table goes.

## Value

A tibble with one row per gene.
