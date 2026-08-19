# Read in Gene Expression and Antibody Capture data from a 10x Genomics Cell Ranger sparse matrix or from a text file.

Read in Gene Expression and Antibody Capture data from a 10x Genomics
Cell Ranger sparse matrix or from a text file.

## Usage

``` r
read_counts_file(path, sample_name, log_file = NULL)
```

## Arguments

- path:

  Path to directory containing 10x matrix, or path to a text file.

- sample_name:

  A character that will be used as a prefix for all cell names.

- log_file:

  Filename for the log file.

## Value

Named list of matrices. One matrix for each data type. Currently the
only two data types are 'Gene Expression' and 'Antibody Capture'
