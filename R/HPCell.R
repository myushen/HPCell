#' HPCell: Massively-Parallel R Native Pipeline for Single-Cell Analysis
#'
#' @importFrom methods is
#' @importFrom stats formula quantile setNames
#' @importFrom utils data head tail
"_PACKAGE"

# Register tidyseurat / tidySingleCellExperiment S3 methods (e.g. left_join.Seurat)
# without import(tidyseurat) in NAMESPACE (avoids dplyr export conflicts on check).
.onLoad <- function(libname, pkgname) {
  loadNamespace("tidyseurat")
  loadNamespace("tidySingleCellExperiment")
}

.myDataEnv <- new.env(parent = emptyenv()) # not exported

.data_internal <- function(dataset) {
  if (!exists(dataset, envir = .myDataEnv)) {
    utils::data(list = c(dataset), envir = .myDataEnv)
  }
}