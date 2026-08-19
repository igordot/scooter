# Add assay to Seurat object.

Add assay to Seurat object.

## Usage

``` r
add_seurat_assay(x, assay, counts_matrix, log_file = NULL)
```

## Arguments

- x:

  Seurat object.

- assay:

  Seurat assay to add the matrix to.

- counts_matrix:

  Raw counts matrix.

- log_file:

  Filename for the log file.

## Value

Seurat object of cells found in both the existing object and new data.
