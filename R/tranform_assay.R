#' Apply a Custom Assay Transformation
#'
#' @description
#' Adds a custom per-cell assay transformation step to the HPCell pipeline.
#' The caller declares one of `identity`, `expm1`, or `safe_expm1` per sample.
#' `identity` does not scale; `expm1` and `safe_expm1` cap at `scale_max`.
#'
#' @param input_hpc An `HPCell` object.
#' @param fx A list of per-sample transformation methods, already chosen
#'   outside HPCell: `"safe_expm1"`, `"expm1"`, or `"identity"`.
#' @param target_input Name of the targets target providing the data object.
#' @param target_output Name of the targets target to write the transformed
#'   data to.
#' @param scale_max A list of numeric upper bounds (one per sample) used by
#'   `expm1` and `safe_expm1`. Ignored for `identity`. Default: `10` per sample.
#' @param ... Additional arguments (unused; for method dispatch).
#' @return The updated `HPCell` object with the transformation step appended.
#' @export
transform_assay <- function(input_hpc, fx = input_hpc$initialisation$input_hpc |> map(~identity), 
                            target_input = "data_object", target_output = "sce_transformed", 
                            # make every element to be 20 as a placeholder unless further indicate
                            scale_max = input_hpc$initialisation$input_hpc |> map(~10), ...) {
  UseMethod("transform_assay")
}

#' @rdname transform_assay
#' @importFrom purrr map
#' @export
transform_assay.HPCell = function(
    input_hpc,
    
    # This might be carrying the environment
    fx = input_hpc$initialisation$input_hpc |> map(~"identity"), 
    target_input = "data_object", 
    target_output = "sce_transformed", 
    scale_max = input_hpc$initialisation$input_hpc |> map(~10),
    ...
) {
  
  fx |> saveRDS("temp_fx.rds")
  scale_max |> saveRDS("count_upper_bound.rds")
  
  input_hpc |> 
    
    # Track the file
    hpc_single("transform_file", "temp_fx.rds", format = "file") |> 
    hpc_iterate(
      target_output = "transform", 
      user_function = readRDS |> quote() ,
      file = "transform_file" |> is_target() 
      # ,
      # iteration = "list", 
      # deployment = "main"
    ) |> 
    
    # Track the file
    hpc_single("count_upper_bound_file", "count_upper_bound.rds", format = "file") |> 
    hpc_iterate(
      target_output = "count_upper_bound", 
      user_function = readRDS |> quote() ,
      file = "count_upper_bound_file" |> is_target() 
    ) |> 
    
    hpc_iterate(
      target_output = target_output, 
      user_function = transform_utility |> quote() , 
      input_read_RNA_assay = "data_object" |> is_target(), 
      transform_fx = "transform" |> is_target()  ,
      external_path = glue("{input_hpc$initialisation$store}/external") |> as.character(),
      container_type = "data_container_type" |> is_target(),
      scale_max = "count_upper_bound" |> is_target()
      
    )
  
}

#' Infer the approximate counts distribution of a sample
#'
#' Decision tree used to infer Census sample distribution and to classify each sample's
#' expression assay before choosing an HPCell transform
#' (`safe_expm1`, `expm1`, or `identity`).
#'
#' Required columns: `has_negative`, `all_integer`, `counts_gap_min_mean`,
#' `positive_mode`, and either `is_gene_expression_likely_log`, `has_rounding_error`
#'
#' @param df A data frame of per-sample count metrics.
#' @param counts_gap_threshold Minimum positive-count gap (min/mean) that
#'   distinguishes double-log from single-log data. Default: `0.25`.
#' @param pos_mode_threshold Positive-count mode must exceed this value (together
#'   with the gap threshold) to call `double_log1p`. Default: `1`.
#'
#' @return `df` with an `inferred_distribution` column. Possible values:
#'   `double_log1p`, `log1p`, `raw_scaled`, `raw`, `log1p_negative_max_10`,
#'   `raw_negative_scaled`, or `NA` when no rule matches.
#' @importFrom rlang .data
#' @export
impute_x_approximate_distribution <- function(df,
                                              counts_gap_threshold,
                                              pos_mode_threshold) {
  df |>
    dplyr::mutate(
      inferred_distribution = dplyr::case_when(
        
        # 0) When counts gap between 0 and next min value >= threshold
        !has_negative & is_gene_expression_likely_log & !all_integer & !has_rounding_error &
          (counts_gap_min_mean >= counts_gap_threshold) & (positive_mode > pos_mode_threshold) ~ "double_log1p",
        
        # 1) Small counts gap
        !has_negative & is_gene_expression_likely_log & !all_integer & !has_rounding_error &
          !(
            (counts_gap_min_mean >= counts_gap_threshold) &
              (positive_mode > pos_mode_threshold)
            
          ) ~ "log1p",
        
        # 2) No negatives, has large values
        !has_negative & !is_gene_expression_likely_log & !all_integer & !has_rounding_error ~ "raw_scaled",
        
        # 3) Large values, integer counts
        !has_negative & !is_gene_expression_likely_log & all_integer & !has_rounding_error ~ "raw",
        
        # 4) Has negatives, compressed range
        has_negative & is_gene_expression_likely_log & !all_integer & !has_rounding_error ~ "log1p_negative_max_10",
        
        # 5) Has negatives and large values
        has_negative & !is_gene_expression_likely_log & !all_integer & !has_rounding_error ~ "raw_negative_scaled",
        
        # fallback
        TRUE ~ NA_character_
      )
    )
}

#' Limit Maximum Value of Counts by Scaling
#'
#' Scales a numeric vector down so that its maximum value does not exceed
#' \code{scale_max}. If the maximum of \code{counts} is already within the
#' limit, the vector is returned unchanged.
#'
#' @param counts A numeric vector of count values to be scaled.
#' @param scale_max A numeric scalar specifying the upper bound for the maximum
#'   value of \code{counts}.
#'
#' @return A numeric vector of the same length as \code{counts}, with all
#'   values scaled so that the maximum does not exceed \code{scale_max}.
#' @export
limit_max_to_scale <- function(counts, scale_max) {
  if (max(counts) > scale_max) {
    scale_factor <- scale_max / max(counts)
    counts <- counts * scale_factor
  }
  counts
}

#' Safe Exponential Minus One Transformation
#'
#' Applies a two-stage \code{\link[base]{expm1}} transformation to a numeric
#' vector with scaling applied before each stage to prevent overflow. Prior to
#' each call to \code{expm1()}, \code{\link{limit_max_to_scale}} is used to
#' ensure the maximum value does not exceed \code{scale_max}, making the
#' transformation numerically stable for large count values.
#'
#' The transformation pipeline is:
#' \deqn{counts \rightarrow \text{scale} \rightarrow expm1 \rightarrow \text{scale} \rightarrow expm1}
#'
#'#' Declared outside HPCell for case 0 (large gap / double log1p).
#'
#' @param counts A numeric vector of count values to be transformed.
#' @param scale_max A numeric scalar passed to \code{\link{limit_max_to_scale}}
#'   specifying the upper bound applied before each \code{expm1} step.
#'
#' @return A numeric vector of the same length as \code{counts} with the
#'   two-stage \code{expm1} transformation applied.
#' @seealso \code{\link{expm1}}, \code{\link{limit_max_to_scale}}
#' @export
safe_expm1 <- function(counts, scale_max) {
  counts <- counts |> limit_max_to_scale(scale_max) |> expm1() |> limit_max_to_scale(scale_max) |> expm1()
  counts
}

# Caller-declared methods from outside HPCell:
# identity  — raw / integer, no scaling
# expm1     — scale to scale_max then inverse log1p
# safe_expm1 — scale → expm1 → scale → expm1
apply_declared_transform <- function(counts, transform_fx, scale_max) {
  fx <- match.fun(transform_fx)
  if (identical(fx, identity)) {
    fx(counts)
  } else if (identical(fx, expm1)) {
    fx(limit_max_to_scale(counts, scale_max))
  } else if (identical(fx, safe_expm1)) {
    fx(counts, scale_max)
  } else {
    stop(
      "Unsupported transform_fx '", deparse(substitute(transform_fx)), "'. ",
      "Inspect further"
    )
  }
}

#' Apply a transformation to an assay and save as HDF5
#'
#' This function applies a caller-declared transformation (`identity`, `expm1`,
#' or `safe_expm1`) to the assay of a SummarizedExperiment object and saves
#' the transformed object.
#'
#' @param input_read_RNA_assay A SummarizedExperiment object to be transformed.
#' @param transform_fx A caller-declared method: `"identity"` (raw, no
#'   scaling), `"expm1"` (scale to `scale_max` then expm1), or `"safe_expm1"`
#'   (scale → expm1 → scale → expm1).
#' @param external_path A character string specifying the directory path to save the transformed object.
#' @param container_type A character vector specifying the output file type. Ideally it should match to the input file type.
#' @param scale_max Upper bound for `limit_max_to_scale()`. Used by `expm1`
#'   and `safe_expm1`; ignored for `identity`.
#' @return The function does not return an object. It saves the transformed SummarizedExperiment object to the specified path.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom dplyr select
#' @importFrom glue glue
#' @importFrom digest digest
#' @importFrom stats density
#'
#' @export
transform_utility  = function(input_read_RNA_assay, transform_fx, 
                              external_path, container_type, scale_max) {
  
  numer_of_cells_to_sample = 5e3
  
  if(ncol(input_read_RNA_assay) == 0) return(NULL)
  
  # Rename assay to "X" for consistency; warn if the original name differs
  if (length(names(assays(input_read_RNA_assay))) == 1 &&
      names(assays(input_read_RNA_assay)) != "X") {
    warning(sprintf(
      "Input assay is named '%s', not 'X'. Renaming to 'X' for downstream consistency.",
      names(assays(input_read_RNA_assay))
    ))
    names(assays(input_read_RNA_assay)) <- "X"
  }
  
  # strip metadata that we don't need
  input_read_RNA_assay =
    input_read_RNA_assay |>
    select(
      any_of(c(".cell", "observation_joinid", "observation_originalid", 
               "donor_id", "dataset_id", "sample_id", "cell_type")),
      starts_with("cell_type")
    )

  # Remove reduced dimensions
  reducedDim(input_read_RNA_assay) = NULL
  
  # Remove row data to avoid downstream binding errors
  rowData(input_read_RNA_assay) <- NULL
  
  # Clear memory  
  gc()
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  # Always use "X" as the canonical assay name regardless of the input assay name
  assay_name <- "X"
  
  # Extract the counts matrix from the first assay of the input
  counts <- assay(input_read_RNA_assay)
  
  # Convert transform_method to a function if it is a character string
  transform_function <- match.fun(transform_fx)
  
  counts <- apply_declared_transform(counts, transform_fx, scale_max)
  
  # Clear memory
  gc()
  
  # This is to avoid memory explosion
  set.seed(42)
  counts_light_for_checks = counts[,sample(seq_len(ncol(counts)), size = min(numer_of_cells_to_sample, ncol(counts))),drop=FALSE]
  
  # Compute the density estimate of the counts. This needs ~13Gb to run for 5000+ cell datasets
  density_est <- counts_light_for_checks |> as.matrix() |> density()
  
  # Clear memory
  gc()
  
  # Find the mode (peak) value of the counts
  mode_value <- density_est$x[which.max(density_est$y)]
  
  # If the mode value is negative, shift counts to make the mode zero
  if (mode_value < 0) {
    counts <- counts + abs(mode_value)
    counts_light_for_checks = counts_light_for_checks + abs(mode_value)
  }
  
  # Round counts to avoid potential subtraction errors due to floating-point precision
  counts <- round(counts, 5)
  counts_light_for_checks = round(counts_light_for_checks, 5)
  
  # Find the most frequent count value (mode) in the counts
  majority_gene_counts <- compute_mode_delayedarray(counts_light_for_checks)$mode
  
  # Subtract the mode value from counts if it is not zero
  if (majority_gene_counts != 0) {
    counts <- counts - majority_gene_counts
    counts_light_for_checks <- counts_light_for_checks - majority_gene_counts
  }
  
  # Replace negative counts with zero to avoid downstream failures
  if (min(counts_light_for_checks) < 0) {
    counts[counts < 0] <- 0
  }
  
  # Clear memory
  gc()
  
  # Assign the modified counts back to the data object
  assay(input_read_RNA_assay, assay_name) <- counts
  
  # Remove cells with zero total counts
  # !!! MAYBE WE SHOULD LKEEP THESE CELLS AND LEAVE THEM TO THE FILTERING STEP
  input_read_RNA_assay <- input_read_RNA_assay[, colSums(counts) > 0]
  
  if (ncol(input_read_RNA_assay) == 0) return(NULL)
  
  # Rebuild the SCE to stay light, and to set the assay with the right name
  input_read_RNA_assay = SingleCellExperiment(
    assays = setNames(list(input_read_RNA_assay |> assay()), assay_name), 
    colData = colData(input_read_RNA_assay)
  )
  
  # Return the modified data object
  input_read_RNA_assay |> 
    
    save_experiment_data(
      dir = glue("{external_path}/{digest(input_read_RNA_assay)}"), 
      container_type = container_type
    )
  
  # extension <- switch(container_type,
  #                     
  #                     "sce_rds" = ".rds",
  #                     "seurat_rds" = ".rds",
  #                     
  #                     "seurat_h5" = ".h5Seurat",
  #                     
  #                     "anndata" = ".h5ad",
  #                     
  #                     "sce_hdf5" = "")
  
  # file_name = paste0(file_name, extension)
  # 
  # # Return data as target instead of file_name pointer
  # 
  # input_read_RNA_assay
  # 
  # 
  # extension <- switch(container_type,
  #                     "sce_rds" = ".rds",
  #                     "seurat_rds" = ".rds",
  #                     "seurat_h5" = ".h5Seurat",
  #                     "anndata" = ".h5ad",
  #                     "sce_hdf5" = "")
  # file_name = paste0(file_name, extension)
  
  # Return data as target instead of file_name pointer
  
}
