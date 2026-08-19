# Run UMAP

The default method runs
[`uwot::umap2()`](https://jlmelville.github.io/uwot/reference/umap2.html)
on a matrix; the Seurat method runs
[`Seurat::RunUMAP()`](https://satijalab.org/seurat/reference/RunUMAP.html)
with `umap.method = "uwot2"` (also
[`uwot::umap2()`](https://jlmelville.github.io/uwot/reference/umap2.html))
and stores a `"umap<suffix>"` reduction. Defaults match
[`Seurat::RunUMAP()`](https://satijalab.org/seurat/reference/RunUMAP.html)
(30 neighbors, minimum distance 0.3) rather than the lower `uwot`
defaults.

## Usage

``` r
run_umap(x, ...)

# Default S3 method
run_umap(x, num_neighbors = 30, min_dist = 0.3, suffix = "", ...)

# S3 method for class 'Seurat'
run_umap(
  x,
  num_neighbors = 30,
  min_dist = 0.3,
  suffix = "",
  assay = "RNA",
  reduction = NULL,
  num_dim = NULL,
  features = NULL,
  var_features = FALSE,
  graph = NULL,
  file_format = c("png", "pdf"),
  ...
)
```

## Arguments

- x:

  A matrix (default) or a Seurat object.

- ...:

  Passed to the underlying function.

- num_neighbors:

  Number of neighbors.

- min_dist:

  Minimum distance between points in the embedding.

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

- graph:

  Graph to use as input (Seurat method).

- file_format:

  File formats to write the scatter plot in (Seurat method).

## Value

A UMAP embedding matrix (default) or a Seurat object with the reduction
added. The Seurat method also writes a sample-colored scatter of the
result to the working directory, as `dr.<reduction>.png`/`.pdf` (or just
`.png` if `file_format = "png"`).
