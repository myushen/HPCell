# This function performs several steps to calculate dispersion:
# Filters out samples with only one instance in the categories of sex or ethnicity to avoid issues with statistical indeterminacy.
# Drops samples that are complete confounders based on sex and ethnicity.
# Calculates dispersion estimates differently based on the size of the dataset:
# For smaller datasets (<1000 samples), it uses all samples to calculate dispersion.
# For larger datasets, it samples a subset of the data (up to 2000 samples) and then calculates dispersion on this subset to reduce computational load.
# Dispersion is calculated using functions from RNA-seq analysis packages (like estimateDisp and estimateTrendedDisp), which handle the complexity of dispersion estimation in gene expression data. This approach is critical in differential expression analysis, as it helps to account for the variability inherent in such data.
# # Define a function 'se_add_dispersion' to add dispersion estimates to single-cell experiment data
#' @param se_df A data frame containing SingleCellExperiment objects.
#' @param my_formula A formula used to model the design matrix.
#' @param my_assay The assay to be used for dispersion estimation.
#' @return The input data frame with updated SingleCellExperiment objects containing dispersion values.
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr enframe
#' @importFrom purrr map
#' @importFrom purrr map_int
#' @importFrom stringr str_subset
#' @importFrom stringr str_c
#' @importFrom stats terms
#' @importFrom stats model.matrix
#' @importFrom magrittr %$%
#' 
#' @export
#' 
se_add_dispersion = function(se, my_formula, my_assay){


        
        # Handle cases where the SingleCellExperiment object is empty
        if(ncol(se) == 0) {
          rowData(se)$dispersion = rep(NA, nrow(se))
        }
        else if(ncol(se) < 1000) {
          # Calculate dispersion for datasets with fewer than 1000 columns
          # Create a model design matrix based on sex and ethnicity
          factors =
            my_formula |> 
            terms() |> 
            attr("term.labels") |> 
            str_subset("\\|", negate = TRUE) |> 
            enframe(value = "factor") |>
            mutate(n = map_int(
              factor, ~ se |> select(.x) |> distinct() |> nrow()
            )) |>
            filter(n > 1) |>
            pull(factor) |>
            str_c(collapse = " + ")
          
          my_design = my_formula |> remove_random_effects() |> model.matrix(data = colData(se) |> droplevels())
          rowData(se)$dispersion = se |> assay(my_assay) |> estimateDisp(design = my_design) %$% tagwise.dispersion
        }
        else {
          # For larger datasets, sample a subset of columns
          sampled_samples = sample(seq_len(ncol(se)), size = min(ncol(se), 2000))
          
          # Create a model design matrix for the sampled subset
          factors =
            my_formula |> 
            terms() |> 
            attr("term.labels") |> 
            str_subset("\\|", negate = TRUE) |> 
            enframe(value = "factor") |>
            mutate(n = map_int(
              factor, ~ se[,sampled_samples, drop=FALSE] |> select(.x) |> distinct() |> nrow()
            )) |>
            filter(n > 1) |>
            pull(factor) |>
            str_c(collapse = " + ")
          
          my_design = my_formula |> remove_random_effects() |> model.matrix(data = se[,sampled_samples, drop=FALSE] |> colData() |> droplevels())
          
          # Estimate trended dispersion based on the sampled subset
          rowData(se)$dispersion =
            assay(se[,sampled_samples, drop=FALSE], my_assay) |>
            estimateTrendedDisp(design = my_design, subset=1000, rowsum.filter=10)
        }
        
        # Return the processed SingleCellExperiment object
        se

  
}



#' Split a SummarizedExperiment object into chunks of specified size
#'
#' This function splits a SummarizedExperiment object into smaller chunks,
#' each containing a specified number of genes.
#'
#' @param se A SummarizedExperiment object to be split.
#' @param chunk_size The number of genes per chunk (default is 100).
#' @return A list of SummarizedExperiment objects, each containing a chunk of genes.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' # Assume `se` is your SummarizedExperiment object
#' # se <- your SummarizedExperiment object here
#' # Split the SummarizedExperiment object into chunks of 100 genes each
#' # se_chunks <- split_summarized_experiment(se)
#' @export
split_summarized_experiment <- function(se, chunk_size = 100) {
  # Get the total number of genes
  total_genes <- nrow(se)
  
  # Calculate the number of chunks
  num_chunks <- ceiling(total_genes / chunk_size)
  
  # Split into chunks
  chunk_list <- lapply(seq_len(num_chunks), function(i) {
    start_index <- (i - 1) * chunk_size + 1
    end_index <- min(i * chunk_size, total_genes)
    
    se_chunk <- se[start_index:end_index, ]
    return(se_chunk)
  })
  
  return(chunk_list)
}


# Define a function 'map_quantile_scale_abundance' for scaling and normalizing abundance in single-cell experiment data
map_quantile_scale_abundance = function(se_df){
  # Print a message indicating the start of the scaling process
  print("Start scale abundance")
  # Perform garbage collection to free up memory
  gc()
  
  # Process the input data frame
  se_df |>
    # Apply quantile normalization to each element of the 'data' column
    mutate(se = map(se, quantile_normalise_abundance, method = "preprocesscore_normalize_quantiles_use_target")) |>
    
    # Convert the scaled counts to a sparse matrix format
    mutate(se = map(se, ~ {
      # Convert the 'counts_scaled' matrix to a sparse matrix
      .x@assays@data$counts_scaled = as(.x@assays@data$counts_scaled, "sparseMatrix")
      
      # Clean up the environment attribute to avoid memory leaks
      attr(.x, "internals")$tt_columns$.abundance_scaled |> attr(".Environment") = NULL # new_environment()
      
      # Return the processed SingleCellExperiment object
      .x
    }))
}


#' @importFrom tidybulk test_differential_abundance
#' @importFrom tidybulk pivot_transcript
#' @importFrom purrr safely
#' 
#' @export
map_de = function(se, my_formula, assay, method, max_rows_for_matrix_multiplication = NULL, cores = 1, .scaling_factor = NULL, .group_by = NULL){
  
  .scaling_factor = enquo(.scaling_factor)
  .group_by = enquo(.group_by)
  
  group_label = se |> distinct(!!.group_by) |> pull(!!.group_by)
  
  # Catch failures
  test_differential_abundance_safe = safely(test_differential_abundance, otherwise = NULL, quiet = TRUE)

    
    # Return prematurely
    if(ncol(se) == 0) return(se)
    
    # Use fast method but does not have dispersion
    if(ncol(se) > 2000 & method |> str_detect("glmmseq")) method = "glmmseq_glmmTMB"
    
    se = se |>
      
      # Test
      test_differential_abundance_safe(
        as.formula(my_formula),
        .abundance = !!sym(assay),
        method = method,
        cores = min(nrow(se), cores),
        max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
        .dispersion = dispersion,
        .scaling_factor = !!.scaling_factor
      ) 
    
    # Throw warning if failed
    # if(se$error |> is.null() |> not())
    #   warning(se$error)
    
    se = se %$% result
    
    # Add the label for joining
    if(se |> is.null() |> not() && .group_by |> quo_is_symbolic())
      se = se |> pivot_transcript() |> mutate(!!.group_by := group_label)
    
    attr(se, "internals")$glmmseq_glmmTMB = NULL
    
    se

  
}


