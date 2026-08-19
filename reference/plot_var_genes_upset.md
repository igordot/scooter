# UpSet plot of the variable genes shared between batches.

UpSet plot of the variable genes shared between batches.

## Usage

``` r
plot_var_genes_upset(var_genes_list, nsets = 50, nintersects = 15)
```

## Arguments

- var_genes_list:

  Named list of variable gene vectors, one element per batch.

- nsets:

  Maximum number of sets to include.

- nintersects:

  Maximum number of intersections to show.

## Value

An UpSet plot.
