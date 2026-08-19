# Calculate differentially expressed genes within each subpopulation/cluster

Calculate differentially expressed genes within each
subpopulation/cluster

## Usage

``` r
differential_expression_per_cluster(
  x,
  cluster_column,
  group_column,
  test = "wilcox",
  out_path = ".",
  write = TRUE,
  log_file = NULL
)
```

## Arguments

- x:

  A Seurat object.

- cluster_column:

  Metadata column specifying the groups to split by.

- group_column:

  Metadata column specifying the groups for differential expressin
  within each split.

- test:

  Statistical method to use.

- out_path:

  Output path.

- write:

  Boolean to save results to disk.

- log_file:

  log file.

## Value

.
