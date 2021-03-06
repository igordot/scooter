% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{as_data_frame_seurat}
\alias{as_data_frame_seurat}
\title{Function to extract data from Seurat object.}
\usage{
as_data_frame_seurat(
  seurat_obj,
  assay = NULL,
  slot = NULL,
  features = NULL,
  reduction = NULL,
  metadata = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object.}

\item{assay}{Assay such as RNA.}

\item{slot}{Slot such as counts. Default is scale.data.}

\item{features}{Features from assay.}

\item{reduction}{Character vector of reduction types.}

\item{metadata}{Boolean. To grab metadata or not}
}
\value{
A metadata file merged on cell identifiers.
}
\description{
Function to extract data from Seurat object.
}
