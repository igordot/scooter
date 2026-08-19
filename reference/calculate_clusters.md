# Identify clusters of cells by graph-based clustering.

Builds a shared nearest neighbor (SNN) graph from reduced dimensions and
identifies clusters over a range of resolutions.

## Usage

``` r
calculate_clusters(
  x,
  num_dim,
  num_neighbors = 20,
  res = NULL,
  algorithm = 4,
  log_file = NULL,
  ...
)

# S3 method for class 'Seurat'
calculate_clusters(
  x,
  assay,
  reduction,
  num_dim,
  num_neighbors = 20,
  res = NULL,
  algorithm = 4,
  log_file = NULL,
  max_clusters = length(get_color_scheme("clusters")),
  ...
)
```

## Arguments

- x:

  A matrix of cell embeddings (cells as rows) or a Seurat object.

- num_dim:

  Number of dimensions to use.

- num_neighbors:

  Number of neighbors (`k.param`) used to build the SNN graph.

- res:

  Clustering resolution(s). If `NULL`, a range of resolutions is used.

- algorithm:

  Clustering algorithm: 1 = Louvain, 2 = Louvain with multilevel
  refinement, 3 = SLM, 4 = Leiden (the default). See
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).

- log_file:

  Filename for the log file.

- ...:

  Arguments passed to the individual methods.

- assay:

  Assay logged as the one this clustering run is over (Seurat method).
  Purely for the log message - clustering itself runs on
  `Embeddings(x, reduction = reduction)`, not on assay data, so nothing
  here reads it back off the object. No default: `DefaultAssay(x)` is
  mutable state that need not match the assay `reduction` actually came
  from, so the caller states it explicitly rather than have the log
  silently report the wrong one.

- reduction:

  Reduction to take the cell embeddings from (Seurat method). No
  default - always state it explicitly.

- max_clusters:

  Drop resolutions yielding this many clusters or more. Defaults to the
  length of the cluster color scheme, since resolutions that exceed it
  cannot be plotted (Seurat method).

## Value

Cluster assignments as a data frame (default method) or the Seurat
object with the cluster assignments added to its metadata (Seurat
method).

## Details

For a Seurat object, the cluster assignments are added to the object
metadata as `res.<resolution>` columns, the labels are converted to the
`C01`, `C02`, ... convention (1-based, zero-padded, and `C`-prefixed to
block downstream numeric coercion), and resolutions that do not yield
any new clusters are dropped.
