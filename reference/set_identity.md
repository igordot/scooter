# Set identity of the Seurat object.

Wrapper for SeuratObject::Idents() with extra safety checks.

## Usage

``` r
set_identity(x, identity_column)
```

## Arguments

- x:

  A Seurat object.

- identity_column:

  The name of the identity column to pull from object metadata.

## Value

A Seurat object with an updated identity set.
