% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{alive_identification}
\alias{alive_identification}
\title{Alive Cell Identification}
\usage{
alive_identification(
  input_read_RNA_assay,
  empty_droplets_tbl = NULL,
  cell_type_ensembl_harmonised_tbl = NULL,
  cell_type_column = NULL,
  assay = NULL,
  feature_nomenclature
)
}
\arguments{
\item{input_read_RNA_assay}{A \code{SingleCellExperiment} or \code{Seurat} object containing RNA assay data.}

\item{empty_droplets_tbl}{A tibble identifying empty droplets.}

\item{cell_type_ensembl_harmonised_tbl}{A tibble with annotated cell type label data.}

\item{cell_type_column}{A character vector indicating the cell type column used for grouping during quality control and dead cell removal.}

\item{assay}{assay used, default = "RNA"}
}
\value{
A tibble identifying alive cells.
}
\description{
\code{alive_identification} filters out dead cells by analyzing mitochondrial and ribosomal gene expression percentages.
}
