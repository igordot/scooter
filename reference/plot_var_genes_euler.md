# Euler diagram of the variable genes shared between batches.

Becomes unreadable and can take a very long time to fit for many sets,
so returns `NULL` when there are too many batches.

## Usage

``` r
plot_var_genes_euler(var_genes_list, color_scheme = NULL, max_sets = 8)
```

## Arguments

- var_genes_list:

  Named list of variable gene vectors, one element per batch.

- color_scheme:

  (optional) Vector of colors.

- max_sets:

  Return `NULL` when there are more sets than this.

## Value

A plot, or `NULL` if there are too many sets to draw.
