# Calculate cluster markers (versus all other clusters, or pairwise) and plot them.

`test`, forwarded to
[`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)/[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html):
`"wilcox"` (Wilcoxon rank sum, the default), `"roc"` (classification
power from 0 = random to 1 = perfect), `"bimod"` (McDavid et al. 2013),
`"tobit"` (Trapnell et al. 2014), or `"MAST"` (Finak et al. 2015).

## Usage

``` r
calculate_cluster_markers(x, label, test, pairwise = FALSE)
```

## Arguments

- x:

  Seurat object with an identity set.

- label:

  Suffix for the output filenames.

- test:

  Statistical test - see above.

- pairwise:

  Compare clusters pairwise instead of each cluster against all others -
  see above.

## Value

`x`, invisibly. Called for its file and plot output.

## Details

`pairwise = TRUE` compares each cluster against each other cluster
individually (rather than each cluster against all others combined),
keeping only genes significant in every comparison - local and global
markers both, at the cost of being much slower for many clusters.
