# Normalize the counts, select variable features, and scale.

Deprecated: renamed to
[`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md) -
"data" didn't say what was being normalized; "counts" matches the
vocabulary already used elsewhere in the package
([`save_counts()`](https://igordot.github.io/scooter/reference/save_counts.md),
`total_counts`,
[`read_counts_file()`](https://igordot.github.io/scooter/reference/read_counts_file.md)).

## Usage

``` r
normalize_data(
  x,
  method = "log",
  num_variable_genes = 3000,
  assay = "RNA",
  log_file = NULL
)
```

## Arguments

- x:

  A Seurat object.

- method:

  "log" for `NormalizeData()` + `FindVariableFeatures()` + `ScaleData()`
  (log-normalized to a scale factor of 10,000), or "sct" for
  `SCTransform()`.

- num_variable_genes:

  Number of variable features.

- assay:

  Assay to normalize.

- log_file:

  Log file.

## Value

The processed Seurat object. See
[`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md).
