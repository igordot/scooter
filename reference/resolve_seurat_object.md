# Resolve a Seurat object, given either the object itself or a path to one on disk.

Meant for parameters that accept either an in-memory Seurat object or a
serialized one - resolving here means callers do not need to check which
they got themselves.

## Usage

``` r
resolve_seurat_object(x)
```

## Arguments

- x:

  A Seurat object, a path to one saved as `.rds` or `.qs2`, or a
  directory containing exactly one such file.

## Value

A Seurat object.
