# Read in Gene Expression and Antibody Capture data from a 10x Genomics Cell Ranger sparse matrix or from a text file.

Deprecated: renamed to
[`read_counts_file()`](https://igordot.github.io/scooter/reference/read_counts_file.md).

## Usage

``` r
load_sample_counts_matrix(sample_name, path, log_file = NULL)
```

## Arguments

- sample_name:

  A character that will be used as a prefix for all cell names.

- path:

  Path to directory containing 10x matrix, or path to a text file.

- log_file:

  Filename for the log file.

## Value

Named list of matrices. One matrix for each data type.
