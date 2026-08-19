# Function to extract data from Seurat object.

Function to extract data from Seurat object.

## Usage

``` r
as_data_frame_seurat(
  x,
  assay = NULL,
  slot = NULL,
  features = NULL,
  reduction = NULL,
  metadata = TRUE
)
```

## Arguments

- x:

  A Seurat object.

- assay:

  Assay such as RNA.

- slot:

  Slot such as counts. Default is scale.data.

- features:

  Features from assay.

- reduction:

  Character vector of reduction types.

- metadata:

  Boolean. To grab metadata or not

## Value

A metadata file merged on cell identifiers.
