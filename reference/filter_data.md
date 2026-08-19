# Filter cells based on the number of genes, counts, and mitochondrial reads.

Deprecated: renamed to
[`filter_cells()`](https://igordot.github.io/scooter/reference/filter_cells.md)
to read as "filter" (verb) + "cells" (object) rather than the vague
"data", and to avoid leading with the bare word "filter", which collides
conceptually with
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)/[`stats::filter()`](https://rdrr.io/r/stats/filter.html).

## Usage

``` r
filter_data(
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
  method).

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

  Minimum number of cells that must survive filtering. `NULL` skips the
  check.

- log_file:

  Log file.

## Value

Filtered data. See
[`filter_cells()`](https://igordot.github.io/scooter/reference/filter_cells.md).
