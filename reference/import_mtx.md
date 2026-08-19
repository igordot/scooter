# Read in 10x Genomics Cell Ranger Matrix Market format data.

Read in 10x Genomics Cell Ranger Matrix Market format data.

## Usage

``` r
import_mtx(data_path, gene_column = 2, log_file = NULL)
```

## Arguments

- data_path:

  Path to directory that holds the files output from 10x.

- gene_column:

  The column with the gene names.

- log_file:

  Filename for the log file.

## Value

Named list of matrices. One matrix for each data type as specified in
the third column of the features.tsv file. As of Oct 3rd 2019, the two
options are `Gene Expression` and `Antibody Capture`
