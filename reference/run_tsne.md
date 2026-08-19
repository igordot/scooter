# Run TSNE

The default method runs
[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) on a
matrix; the Seurat method runs
[`Seurat::RunTSNE()`](https://satijalab.org/seurat/reference/RunTSNE.html)
on the object and stores a `"tsne<suffix>"` reduction.

## Usage

``` r
run_tsne(x, ...)

# Default S3 method
run_tsne(x, seed.use = 1, dim.embed = 2, suffix = "", ...)

# S3 method for class 'Seurat'
run_tsne(
  x,
  seed.use = 1,
  dim.embed = 2,
  suffix = "",
  assay = "RNA",
  reduction = NULL,
  num_dim = NULL,
  features = NULL,
  var_features = FALSE,
  file_format = c("png", "pdf"),
  ...
)
```

## Arguments

- x:

  A matrix (default) or a Seurat object.

- ...:

  Passed to the underlying function.

- seed.use:

  Seed to use.

- dim.embed:

  Number of tSNE embeddings to return.

- suffix:

  Tag inserted into the reduction/column names.

- assay:

  Assay to use (Seurat method). Always set explicitly.

- reduction:

  Existing reduction to take dimensions from (Seurat method).

- num_dim:

  Number of dimensions to take from `reduction` (Seurat method).

- features:

  Explicit features to use (Seurat method).

- var_features:

  Use the assay's variable features (Seurat method).

- file_format:

  File formats to write the scatter plot in (Seurat method).

## Value

A tSNE embedding matrix (default) or a Seurat object with the reduction
added. The Seurat method also writes a sample-colored scatter of the
result to the working directory, as `dr.<reduction>.png`/`.pdf` (or just
`.png` if `file_format = "png"`).
