#' sanity check Samples with Nonsensical Count Values After Transformation
#'
#' @description
#' Adds a per-sample sanity-check step to the HPCell pipeline, followed by a
#' cross-sample merge and an HTML quality-control report.  Three independent
#' binary sanity checks are computed for each sample:
#'
#' \describe{
#'   \item{`sanity_check_max_lt_10`}{Maximum count value is < 10, suggesting the data
#'     were never back-transformed out of log/fractional space.}
#'   \item{`sanity_check_min_lt_0`}{Minimum count value is negative, which should not
#'     occur after the mode-subtraction and floor-clipping in
#'     \code{\link{transform_utility}}.}
#'   \item{`sanity_check_rounding_error`}{Counts are not all exact integers but every
#'     value is within \eqn{10^{-4}} of an integer — tiny floating-point
#'     artefacts from numerical operations.  Exact-integer samples are excluded.}
#' }
#'
#' A sample can carry any combination of sanity checks.  The resulting target
#' (`target_output`) is a list of 1-row tibbles (or `NULL` for empty samples);
#' an aggregate tibble is written to `<target_output>_combined`, and an HTML
#' report is written to `<target_output>_report`.
#'
#' @param input_hpc An \code{HPCell} object.
#' @param target_input Name of the targets target providing the transformed
#'   \code{SingleCellExperiment} objects. Default: \code{"sce_transformed"}.
#' @param target_output Base name for the output targets.  The pipeline appends
#'   `_combined` for the merged tibble and `_report` for the HTML report.
#'   Default: \code{"sanity_checked_sample_stats"}.
#' @param n_cells_sample Integer. Maximum number of cells to subsample per
#'   sample when computing the check statistics. Default: \code{5000L}.
#' @param ... Additional arguments (unused; for method dispatch).
#'
#' @return The updated \code{HPCell} object with three new pipeline steps
#'   appended (iterate, merge, report).
#' @export
sanity_check_transform_samples <- function(input_hpc,
                                   target_input  = "sce_transformed",
                                   target_output = "sanity_checked_sample_stats",
                                   n_cells_sample = 5000L,
                                   ...) {
  UseMethod("sanity_check_transform_samples")
}

#' @rdname sanity_check_transform_samples
#' @importFrom SummarizedExperiment assay colData
#' @importFrom glue glue
#' @importFrom readr write_lines
#' @importFrom purrr set_names
#' @importFrom here here
#' @export
sanity_check_transform_samples.HPCell <- function(input_hpc,
                                          target_input  = "sce_transformed",
                                          target_output = "sanity_checked_sample_stats",
                                          n_cells_sample = 5000L,
                                          ...) {
  target_combined <- paste0(target_output, "_combined")
  target_report   <- paste0(target_output, "_report")
  rmd_path        <- system.file("rmd/sanity_checked_sample_stats_report.qmd", package = "HPCell")
  
  input_hpc |>
    
    # 1. Per-sample sanity check tibble
    hpc_iterate(
      target_output  = target_output,
      user_function  = sanity_check_counts_utility |> quote(),
      sce            = target_input |> is_target(),
      n_cells_sample = n_cells_sample
    ) |>
    
    # 2. Merge all per-sample tibbles into one
    hpc_merge(
      target_output = target_combined,
      user_function = sanity_check_transform_merge |> quote(),
      sanity_check_list     = target_output |> is_target()
    ) |>
    
    # 3. Render QC report via tarchetypes::tar_render (rmarkdown, no Quarto needed)
    append_sanity_check_report(
      target_output = target_report,
      rmd_path      = rmd_path,
      sanity_check_tbl_target = target_combined
    )
}

# Private helper: appends a tarchetypes::tar_quarto_raw target to the pipeline
# script.  It also sniffs the Quarto binary location at pipeline-construction
# time and embeds a Sys.setenv(QUARTO_PATH = ...) line so that worker
# processes (which may not inherit the interactive PATH) can find the CLI.
append_sanity_check_report <- function(input_hpc, target_output, rmd_path, sanity_check_tbl_target) {
  target_script <- glue("{input_hpc$initialisation$store}.R")
  external_dir  <- glue("{input_hpc$initialisation$store}/external") |> as.character()
  dir.create(external_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Copy the .qmd template into the writable external dir.
  # Quarto's --output-file sanity check only accepts a filename (not a path), so we
  # point path= to the copy; Quarto then writes the .html right alongside it.
  local_qmd <- file.path(external_dir, paste0(target_output, ".qmd"))
  file.copy(rmd_path, local_qmd, overwrite = TRUE)
  
  target_output |> delete_lines_with_word(target_script)
  
  # deployment = "main" forces the render on the login/RStudio node where
  # Quarto is installed, not on Slurm compute workers.
  glue(
    'target_list = c(target_list, tarchetypes::tar_quarto_raw(',
    '  "{target_output}",',
    '  path = "{local_qmd}",',
    '  execute_params = quote(list(sanity_check_tbl = {sanity_check_tbl_target})),',
    '  deployment = "main",',
    '  quiet = FALSE',
    '))',
    .sep = "\n"
  ) |>
    readr::write_lines(target_script, append = TRUE)
  
  input_hpc |>
    c(list(list(target_output = target_output, iterate = "single")) |>
        set_names(target_output)) |>
    add_class("HPCell")
}

#' Per-Sample Count Sanity Check
#'
#' Worker function called once per sample by the
#' \code{\link{sanity_check_transform_samples}} pipeline step.  It subsamples up to
#' \code{n_cells_sample} cells and evaluates three independent binary sanity checks on
#' the post-mode-subtraction counts stored in \code{sce}:
#'
#' \describe{
#'   \item{`sanity_check_max_lt_10`}{`max(counts) < 10`}
#'   \item{`sanity_check_min_lt_0`}{`min(counts) < 0`}
#'   \item{`sanity_check_rounding_error`}{Counts are not all exact integers but every
#'     value is within \eqn{10^{-4}} of an integer — tiny floating-point
#'     artefacts.  Exact-integer samples are excluded via an intermediate
#'     \code{all_int} guard so this sanity check only fires on data that has been
#'     through numerical operations introducing sub-\eqn{10^{-4}} residuals.}
#' }
#'
#' @param sce A \code{SingleCellExperiment} produced by
#'   \code{\link{transform_utility}}, or \code{NULL}.
#' @param n_cells_sample Integer. Maximum number of cells to subsample for the
#'   check. Default: \code{5000L}.
#'
#' @return A one-row \code{\link[tibble]{tibble}} with columns \code{sample_id},
#'   \code{sanity_check_max_lt_10}, \code{sanity_check_min_lt_0}, and
#'   \code{sanity_check_rounding_error} (all sanity checks are 0/1 integers), or \code{NULL}
#'   when \code{sce} is \code{NULL} or has no cells.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom tibble tibble
#' @export
sanity_check_counts_utility <- function(sce, n_cells_sample = 5000L) {
  if (is.null(sce) || ncol(sce) == 0L) return(NULL)
  
  set.seed(42)
  idx         <- sample(seq_len(ncol(sce)), size = min(as.integer(n_cells_sample), ncol(sce)))
  counts_vec  <- as.numeric(assay(sce)[, idx, drop = FALSE])
  
  max_val  <- max(counts_vec, na.rm = TRUE)
  min_val  <- min(counts_vec, na.rm = TRUE)
  tol      <- 1e-4
  
  all_int  <- all(counts_vec == floor(counts_vec), na.rm = TRUE)
  
  tibble::tibble(
    sample_id           = unique(colData(sce)[["sample_id"]]),
    sanity_check_max_lt_10      = as.integer(max_val < 10),
    sanity_check_min_lt_0       = as.integer(min_val < 0),
    sanity_check_rounding_error = as.integer(
      !all_int && all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE)
    )
  )
}

#' Merge Per-Sample sanity check Tibbles into a Single Table
#'
#' Utility function used by the \code{\link{sanity_check_transform_samples}} merge step.
#' It discards \code{NULL} entries (empty / skipped samples) and row-binds the
#' remaining one-row tibbles.
#'
#' @param sanity_check_list A list of one-row \code{\link[tibble]{tibble}}s as returned
#'   by \code{\link{sanity_check_counts_utility}}, possibly containing
#'   \code{NULL} elements.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per non-null sample and
#'   columns \code{sample_id}, \code{sanity_check_max_lt_10}, \code{sanity_check_min_lt_0},
#'   \code{sanity_check_rounding_error}.
#'
#' @importFrom purrr compact
#' @importFrom dplyr bind_rows
#' @export
sanity_check_transform_merge <- function(sanity_check_list) {
  sanity_check_list |> purrr::compact() |> dplyr::bind_rows()
}
