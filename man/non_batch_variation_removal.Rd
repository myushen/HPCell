% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{non_batch_variation_removal}
\alias{non_batch_variation_removal}
\title{Non-Batch Variation Removal}
\usage{
non_batch_variation_removal(
  input_read_RNA_assay,
  empty_droplets_tbl,
  alive_identification_tbl,
  cell_cycle_score_tbl,
  assay = NULL,
  factors_to_regress = NULL
)
}
\arguments{
\item{input_read_RNA_assay}{Path to demultiplexed data.}

\item{empty_droplets_tbl}{Path to empty droplets data.}

\item{alive_identification_tbl}{A tibble from alive cell identification.}

\item{cell_cycle_score_tbl}{A tibble from cell cycle scoring.}

\item{assay}{assay used, default = "RNA"}
}
\value{
Normalized and adjusted data.
}
\description{
Regresses out variations due to mitochondrial content, ribosomal content, and
cell cycle effects.
}
