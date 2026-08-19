# Run PCA

The default method runs
[`irlba::prcomp_irlba()`](https://rdrr.io/pkg/irlba/man/prcomp_irlba.html)
on a matrix; the Seurat method runs
[`Seurat::RunPCA()`](https://satijalab.org/seurat/reference/RunPCA.html)
on the object and stores a `"pca<suffix>"` reduction.

## Usage

``` r
run_pca(x, ...)

# Default S3 method
run_pca(x, num_pcs = 50, suffix = "", ...)

# S3 method for class 'Seurat'
run_pca(
  x,
  num_pcs = 50,
  suffix = "",
  assay = "RNA",
  features = NULL,
  var_features = FALSE,
  ...
)
```

## Arguments

- x:

  A matrix (default) or a Seurat object.

- ...:

  Passed to the underlying function.

- num_pcs:

  Number of principal components.

- suffix:

  Tag inserted into the reduction/column names.

- assay:

  Assay to use (Seurat method). Always set explicitly.

- features:

  Explicit features to use (Seurat method).

- var_features:

  Use the assay's variable features (Seurat method).

## Value

A named list (default) or a Seurat object with the reduction added
(Seurat method). The Seurat method also writes the diagnostic plots to
the working directory: a sample-colored scatter and an elbow plot, plus
a heatmap past 15 components.
