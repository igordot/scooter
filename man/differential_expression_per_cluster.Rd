% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de.R
\name{differential_expression_per_cluster}
\alias{differential_expression_per_cluster}
\title{Calculate differentially expressed genes within each subpopulation/cluster}
\usage{
differential_expression_per_cluster(
  seurat_obj,
  cluster_column,
  group_column,
  test = "wilcox",
  out_path = ".",
  write = TRUE,
  log_file = NULL
)
}
\arguments{
\item{seurat_obj}{Gene expression data.}

\item{cluster_column}{Metadata column specifying the groups to split by.}

\item{group_column}{Metadata column specifying the groups for differential expressin within each split.}

\item{test}{Statistical method to use.}

\item{out_path}{Output path.}

\item{write}{Boolean to save results to disk.}

\item{log_file}{log file.}
}
\value{
.
}
\description{
Calculate differentially expressed genes within each subpopulation/cluster
}
