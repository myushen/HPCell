#' Flag Samples with Nonsensical Count Values After Transformation
#'
#' @description
#' Adds a sanity-check step to the HPCell pipeline that inspects the count
#' distribution of each transformed sample. A sample is flagged when its
#' maximum count is at most \code{max_count_threshold} **and** the majority of
#' non-zero counts are below 1 — a pattern that indicates the data were never
#' properly back-transformed out of log/fractional space.
#'
#' The check mirrors the light subsample already computed inside
#' \code{\link{transform_utility}} (up to \code{n_cells_sample} randomly
#' chosen cells, evaluated after mode-subtraction).
#'
#' @param input_hpc An \code{HPCell} object.
#' @param target_input Name of the targets target providing the transformed
#'   \code{SingleCellExperiment} objects. Default: \code{"sce_transformed"}.
#' @param target_output Name of the targets target to write flagged sample IDs
#'   to. Default: \code{"flagged_sample_id"}.
#' @param n_cells_sample Integer. Maximum number of cells to subsample per
#'   sample when computing the check statistics. Default: \code{5000L}.
#' @param max_count_threshold Numeric. Samples whose count maximum exceeds this
#'   value are considered sensical and are **not** flagged.  Default: \code{10}.
#' @param frac_below_one Numeric in \eqn{(0, 1)}. Minimum fraction of non-zero
#'   counts that must be below 1 (together with \code{max_count_threshold})
#'   to trigger a flag. Default: \code{0.5}.
#' @param ... Additional arguments.
#'
#' @return The updated \code{HPCell} object with the sanity-check step
#'   appended.  The new target (\code{target_output}) is a list whose elements
#'   are either the \code{sample_id} character value(s) for flagged samples or
#'   \code{NULL} for clean samples.
#' @export
flag_transform_samples <- function(input_hpc,
                                   target_input  = "sce_transformed",
                                   target_output = "flagged_sample_id",
                                   n_cells_sample    = 5000L,
                                   max_count_threshold = 10,
                                   frac_below_one    = 0.5,
                                   ...) {
  UseMethod("flag_transform_samples")
}

#' @rdname flag_transform_samples
#' @importFrom SummarizedExperiment assay colData
#' @export
flag_transform_samples.HPCell <- function(input_hpc,
                                          target_input  = "sce_transformed",
                                          target_output = "flagged_sample_id",
                                          n_cells_sample    = 5000L,
                                          max_count_threshold = 10,
                                          frac_below_one    = 0.5,
                                          ...) {
  input_hpc |>
    hpc_iterate(
      target_output = target_output,
      user_function = flag_counts_sanity_utility |> quote(),
      sce                 = target_input |> is_target(),
      n_cells_sample      = n_cells_sample,
      max_count_threshold = max_count_threshold,
      frac_below_one      = frac_below_one
    )
}

#' Check Whether a Transformed Sample Has Nonsensical Count Values
#'
#' Worker function called by the \code{\link{flag_transform_samples}} pipeline
#' step.  It subsamples up to \code{n_cells_sample} cells, then checks two
#' conditions simultaneously:
#'
#' \enumerate{
#'   \item The maximum count in the subsample is \eqn{\le} \code{max_count_threshold}.
#'   \item More than \code{frac_below_one} of non-zero counts are \eqn{< 1}.
#' }
#'
#' Both conditions must hold for the sample to be flagged.  This detects
#' datasets that were never back-transformed from log/fractional space (values
#' mostly in \eqn{[0, 1]}, max \eqn{\le 10}) and would produce nonsensical
#' downstream results in CellNexus.
#'
#' @param sce A \code{SingleCellExperiment} after transformation produced by
#'   \code{\link{transform_utility}}, or \code{NULL}.
#' @param n_cells_sample Integer. Maximum number of cells to subsample.
#'   Default: \code{5000L}.
#' @param max_count_threshold Numeric. Counts with a maximum above this value
#'   are not flagged. Default: \code{10}.
#' @param frac_below_one Numeric. Fraction threshold for non-zero values
#'   below 1. Default: \code{0.5}.
#'
#' @return A character vector of unique \code{sample_id} values from
#'   \code{colData(sce)} when the sample is flagged; otherwise \code{NULL}.
#'
#' @importFrom SummarizedExperiment assay colData
#' @export
flag_counts_sanity_utility <- function(sce,
                                       n_cells_sample      = 5000L,
                                       max_count_threshold = 10,
                                       frac_below_one      = 0.5) {
  if (is.null(sce) || ncol(sce) == 0L) return(NULL)
  
  set.seed(42)
  idx <- sample(seq_len(ncol(sce)),
                size = min(as.integer(n_cells_sample), ncol(sce)))
  counts_light <- assay(sce)[, idx, drop = FALSE]
  
  # Materialise the subsampled block – it is already small enough to be safe
  counts_vec <- as.numeric(counts_light)
  
  max_val  <- max(counts_vec, na.rm = TRUE)
  nonzero  <- counts_vec[counts_vec > 0]
  frac_lt1 <- if (length(nonzero) > 0L) mean(nonzero < 1) else 0
  
  is_nonsensical <- (max_val <= max_count_threshold) && (frac_lt1 > frac_below_one)
  
  if (is_nonsensical) {
    message(
      "[flag_transform_samples] Sample flagged: max = ", round(max_val, 4),
      ", fraction non-zero < 1 = ", round(frac_lt1, 4)
    )
    unique(colData(sce)[["sample_id"]])
  } else {
    NULL
  }
}
