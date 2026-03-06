# Step3
# Group samples by dataset_id, cell_type

# This script sets up a robust and scalable data processing pipeline for single-cell RNA sequencing (scRNA-seq) datasets using the targets package in R, which facilitates reproducible and efficient workflows. Specifically, the code orchestrates the ingestion and preprocessing of multiple SingleCellExperiment objects corresponding to different datasets (dataset_id) and targets (target_name). It leverages high-performance computing resources through the crew package, configuring multiple SLURM-based controllers (tier_1 to tier_4) to handle varying computational loads efficiently.
# 
# The pipeline performs several key steps:
#   
#   1.	Data Retrieval: It reads raw SingleCellExperiment objects for each target, ensuring that only successfully loaded data proceeds further.
# 2.	Normalization: Calculates Counts Per Million (CPM) for each cell to normalize gene expression levels across cells and samples.
# 3.	Data Aggregation: Groups the data by dataset_id and tar_group, then combines the SingleCellExperiment objects within each group into a single object, effectively consolidating the data for each dataset.
# 4.	Metadata Integration: Joins additional metadata, such as cell types, by connecting to a DuckDB database and fetching relevant information from a Parquet file. This enriches the single-cell data with essential annotations.
# 5.	Cell Type Segmentation: Splits the combined SingleCellExperiment objects into separate objects based on cell_type, facilitating downstream analyses that are specific to each cell type.
# 6.	Data Saving with Error Handling: Generates unique identifiers for each cell type within a dataset and saves both the raw counts and CPM-normalized data to specified directories. It includes special handling for cases where a cell type has only one cell, duplicating the data to prevent errors during the saving process.
# 
# By integrating targets, crew, and various data manipulation packages (dplyr, tidyverse, SingleCellExperiment), this script ensures that large-scale scRNA-seq data processing is efficient, reproducible, and capable of leveraging parallel computing resources. It is designed to handle edge cases gracefully and provides a clear framework for preprocessing scRNA-seq data, which is essential for subsequent analyses such as clustering, differential expression, and cell type identification.


library(arrow)
library(dplyr)
library(duckdb)

#


# Get Dharmesh metadata consensus
#system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/cell_annotation_new_substitute_cell_type_na_to_unknown.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/")

job::job({
  
  get_file_ids = function(cell_annotation 
                          #cell_type_consensus_parquet
  ){
    
    # cell_consensus = 
    #   tbl(
    #     dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    #     sql(glue::glue("SELECT * FROM read_parquet('{cell_type_consensus_parquet}')"))
    #   ) |>
    #   select(cell_, dataset_id, cell_type_unified_ensemble, cell_type_unified) 
    
    # This because f7c1c579-2dc0-47e2-ba19-8165c5a0e353 includes 13K samples
    # It affects only very few datasets
    sample_chunk_df = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
      ) |> 
      # Define chunks
      dplyr::count(dataset_id, sample_id, name = "cell_count") |>  # Ensure unique dataset_id and sample_id combinations
      distinct(dataset_id, sample_id, cell_count) |>  # Ensure unique dataset_id and sample_id combinations
      group_by(dataset_id) |> 
      dbplyr::window_order(dataset_id, cell_count, sample_id) |>  # Ensure order. Note: order cell_count only is not enough because it needs a secondary tie-breaker
      mutate(sample_index = row_number()) |>  # Create sequential index within each dataset
      mutate(sample_chunk = (sample_index - 1) %/% 1000 + 1) |>  # Assign chunks (up to 1000 samples per chunk)
      mutate(sample_pseudobulk_chunk = (sample_index - 1) %/% 250 + 1) |> # Max combination of dataset_id, sample_pseudobulk_chunk and file_id_pseudobulk up to 10000
      mutate(cell_chunk = cumsum(cell_count) %/% 100000 + 1) |> # max 20K cells per sample
      ungroup() 
    
    # Test whether cell_chunk and sample_chunk are unique for this sample
    run_chunk_once <- function(column_name, id) {
      sample_chunk_df |> filter(sample_id == id) |> pull(!!column_name)
    }
    
    sample_chunk_results <- replicate(20, run_chunk_once("sample_chunk", "d6e942a09a140ee8bb6f0c3da8defea4___exp7-human-150well."), simplify = FALSE)
    sample_chunk_identical <- all(sapply(sample_chunk_results[-1], function(x) identical(x, sample_chunk_results[[1]])))
    if (!sample_chunk_identical) {
      stop("Inconsistent sample chunk value was generated in multiple runs, this will lead to file id changes")
    }
    
    cell_chunk_results <- replicate(20, run_chunk_once("cell_chunk",  "d6e942a09a140ee8bb6f0c3da8defea4___exp7-human-150well."), simplify = FALSE)
    cell_chunk_identical <- all(sapply(cell_chunk_results[-1], function(x) identical(x, cell_chunk_results[[1]])))
    if (!cell_chunk_identical) {
      stop("Inconsistent cell chunk value was generated in multiple runs, this will lead to file id changes")
    }
    
    
    tbl(
      dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
      sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
    ) |> 
      # left_join(cell_consensus, copy=TRUE) |>
      
      # Cells in cell_annotation could be more than cells in cell_consensus. In order to avoid NA happens in cell_consensus cell_type column
      mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(),
                                                 "Unknown",
                                                 cell_type_unified_ensemble)) |>
      
      left_join(sample_chunk_df |> select(dataset_id, sample_chunk, sample_pseudobulk_chunk, cell_chunk, sample_id), copy=TRUE) |> 
      # # Make sure I cover cell type even if consensus of harmonisation is not present (it should be the vast minority)
      # mutate(temp_cell_type_label_for_file_id = if_else(cell_type_unified_ensemble |> is.na(), cell_type, cell_type_unified)) |> 
      # mutate(temp_cell_type_label_for_file_id = if_else(temp_cell_type_label_for_file_id |> is.na(), cell_type, temp_cell_type_label_for_file_id)) |> 
      # 
      # Define chunks
      group_by(dataset_id, sample_chunk, cell_chunk, sample_pseudobulk_chunk, cell_type_unified_ensemble, sample_id) |>
      summarise(cell_count = n(), .groups = "drop") |>
      group_by(dataset_id, sample_chunk, cell_chunk, cell_type_unified_ensemble) |>
      dbplyr::window_order(desc(cell_count)) |>
      mutate(chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
      ungroup() |> 
      as_tibble() |> 
      
      # Single cell file ID
      mutate(file_id_cellNexus_single_cell = 
               glue::glue("{dataset_id}___{sample_chunk}___{cell_chunk}___{cell_type_unified_ensemble}") |> 
               sapply(digest::digest) |> 
               paste0("___", chunk, ".h5ad") 
      ) |> 
      
      # seudobulk file id
      #mutate(file_id_cellNexus_pseudobulk = paste0(dataset_id, ".h5ad"))
      mutate(file_id_cellNexus_pseudobulk = 
               glue::glue("{dataset_id}___{sample_pseudobulk_chunk}") |> 
               sapply(digest::digest) |>
               paste0("___", chunk, ".h5ad"))
    
  }
  
  
  # FOR MENGYUAN CELL_METADATA COULD BE BIGGER THAN CELL_ANNOTATION
  
  get_file_ids(
    "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet" # MODIFY HERE: input cell annotation parquet
  )  |> 
    write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cellNexus_single_cell_2024_Jul.parquet") # MODIFY HERE: output file_id parquet
  
  gc()
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  dir.create("/vast/scratch/users/shen.m/duckdb_tmp", showWarnings = FALSE) # MODIFY HERE: duckdb temp directory
  
  DBI::dbExecute(
    con,
    "SET temp_directory='/vast/scratch/users/shen.m/duckdb_tmp';" # MODIFY HERE: duckdb temp directory (must match dir.create above)
  )
  
  # Create a view for cell_annotation in DuckDB
  # MODIFY HERE: cell_metadata parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT 
    CONCAT(cell_, '___', dataset_id) AS cell_,
    dataset_id,
    *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata.parquet')
")
  
  #   dbExecute(con, "
  #   CREATE VIEW cell_annotation AS
  #   SELECT cell_, blueprint_first_labels_fine, monaco_first_labels_fine, azimuth_predicted_celltype_l2
  #   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/annotation_tbl_light.parquet')
  # ")
  
  # MODIFY HERE: cell_annotation parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW empty_droplet_df AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation_2024_Jul.parquet')
")
  
  # MODIFY HERE: file_id parquet path inside the SQL string below (should match the write_parquet output above)
  dbExecute(con, "
  CREATE VIEW file_id_cellNexus_single_cell AS
  SELECT dataset_id, sample_chunk, cell_chunk, sample_pseudobulk_chunk, cell_type_unified_ensemble, sample_id, file_id_cellNexus_single_cell, file_id_cellNexus_pseudobulk
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cellNexus_single_cell_2024_Jul.parquet')
")
  
  #   # This DF is needed to filter out unmatched sample-cell-type combo. Otherwise, cellNexus get_pseudobulk will slice cell names out of bounds.
  #   dbExecute(con, "
  #   CREATE VIEW sample_cell_type_combo AS
  #   SELECT dataset_id, sample_id, cell_type_unified_ensemble
  #   FROM read_parquet('/vast/scratch/users/shen.m/Census_final_run/cell_type_concensus_tbl_from_hpcell.parquet')
  # ")
  
  # Perform the left join and save to Parquet
  copy_query <- "
  COPY (
     SELECT 
        cell_metadata.cell_ AS cell_id, -- Rename cell_ to cell_id
        COALESCE(empty_droplet_df.alive, FALSE) AS alive, -- Set alive column NULL to FALSE
        cell_metadata.*,              -- Include all other columns from cell_metadata
        empty_droplet_df.*,           -- Include all columns from empty_droplet_df
        file_id_cellNexus_single_cell.*, -- Include all columns from file_id_cellNexus_single_cell
        atlas_id                      -- Specify the atlas name 
      FROM cell_metadata

      LEFT JOIN empty_droplet_df
        ON empty_droplet_df.cell_ = cell_metadata.cell_
        AND empty_droplet_df.dataset_id = cell_metadata.dataset_id
    
      LEFT JOIN file_id_cellNexus_single_cell
        ON file_id_cellNexus_single_cell.sample_id = empty_droplet_df.sample_id
        AND file_id_cellNexus_single_cell.dataset_id = empty_droplet_df.dataset_id
        AND file_id_cellNexus_single_cell.cell_type_unified_ensemble = empty_droplet_df.cell_type_unified_ensemble
        
      WHERE cell_metadata.dataset_id NOT IN ('99950e99-2758-41d2-b2c9-643edcdf6d82', '9fcb0b73-c734-40a5-be9c-ace7eea401c9') -- (THESE TWO DATASETS DOESNT contain meaningful data - no observation_joinid etc), thus was excluded in the final metadata.
         
  ) TO  '/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_2_2_mengyuan.parquet' -- MODIFY HERE: output merged metadata parquet (v1_2_2)
  (FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  #system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_4.parquet box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
  
  print("Done.")
})

# We decided to make cell_id lighter without re-run everything in HPCell pipeline. Here to swap cell_id in the metadata
# cell_map is processed in a separate target script dataset_cell_map.R
job::job({
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create a view for cell_annotation in DuckDB
  # MODIFY HERE: v1_2_2 merged metadata parquet path inside the SQL string below (should match the COPY TO output above)
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_2_2_mengyuan.parquet')
")
  
  # MODIFY HERE: cell_id dictionary parquet path inside the SQL string below
  dbExecute(con, "
  CREATE VIEW cell_map AS
  SELECT *
  FROM read_parquet('/vast/projects//cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cell_id_dict_v1_0_0_Jul_2024.parquet')
")
  
  # Perform the left join and save to Parquet
  copy_query <- "
  COPY (
     SELECT *
      FROM cell_metadata
    
      LEFT JOIN cell_map
        ON cell_metadata.cell_id = cell_map.cell_id
        AND cell_metadata.file_id_cellNexus_single_cell = cell_map.file_id_cellNexus_single_cell

  ) TO  '/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_3_2_mengyuan.parquet' -- MODIFY HERE: output final metadata parquet with new cell IDs (v1_3_2)
  (FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  #system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_4.parquet box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
  
  print("Done.")
  
  
})



# MODIFY HERE: final metadata parquet path used for the targets pipeline (should match the COPY TO output above)
cell_metadata = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_3_2_mengyuan.parquet')")
  )

library(targets)
library(tidyverse)
store_file_cellNexus = "/vast/scratch/users/shen.m/targets_prepare_database_split_datasets_chunked_1_3_0_single_cell" # MODIFY HERE: targets store directory for this pipeline

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "never"),
    #debug = "dataset_id_sce", 
    
    workspace_on_error = TRUE,
    controller = crew_controller_group(
      list(
        crew_controller_slurm(
          name = "elastic",
          workers = 300,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_error = 10,
          options_cluster = crew_options_slurm(
            #memory_gigabytes_required = c(10, 20, 40, 80, 160), 
            memory_gigabytes_required = c(30, 45, 65, 80, 160), 
            cpus_per_task = c(2, 2, 5, 10, 20), 
            time_minutes = c(30, 30, 30, 60*4, 60*24),
            verbose = T
          )
        ),
        
        crew_controller_slurm(
          name = "tier_1", 
          script_lines = "#SBATCH --mem 8G",
          slurm_cpus_per_task = 1, 
          workers = 300, 
          tasks_max = 50,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        
        crew_controller_slurm(
          name = "tier_2",
          script_lines = "#SBATCH --mem 10G",
          slurm_cpus_per_task = 1,
          workers = 300,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_3",
          script_lines = "#SBATCH --mem 20G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_4",
          workers = 200,
          tasks_max = 10,
          crashes_error = 5, 
          seconds_idle = 30,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = c(60, 90, 120, 240), 
            cpus_per_task = c(2), 
            time_minutes = c(60*24),
            verbose = T
          )
        ),
        crew_controller_slurm(
          name = "tier_5",
          script_lines = "#SBATCH --mem 150G",
          slurm_cpus_per_task = 1,
          workers = 2,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        )
      )
    ), 
    trust_object_timestamps = TRUE
    #workspaces = "dataset_id_sce_52dbec3c15f98d66"
  )
  
  
  save_anndata = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    # # Parallelise
    # cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    # bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    # 
    
    
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "counts"
    
    # # Check if the 'sce' has only one cell (column)
    # if(ncol(assay(.x)) == 1) {
    #   
    #   # Duplicate the assay to prevent saving errors due to single-column matrices
    #   my_assay = cbind(assay(.x), assay(.x))
    #   # Rename the second column to distinguish it
    #   colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
    #   
    #   cd = colData(.x)
    #   cd = cd |> rbind(cd)
    #   rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
    #   
    #   
    #   
    #   .x =  SingleCellExperiment(assay = list( counts = my_assay ), colData = cd) 
    # } 
    # 
    # 
    # # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
    # .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
    # 
    # .x =  SingleCellExperiment(assay = list( counts = .x |> assay()), colData = colData(.x)) 
    
    
    # My attempt to save a integer, sparse, delayed matrix (with zellkonverter it is not possible to save integers)
    # .x |> assay() |> type() <- "integer"
    # .x |> saveHDF5SummarizedExperiment("~/temp", as.sparse = T, replace = T)
    
    # Save the experiment data to the specified counts cache directory
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    return(TRUE)  # Indicate successful saving
    
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata <- purrr::insistently(save_anndata, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  save_anndata_cpm = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    # # Parallelise
    # cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    # bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    # 
    dataset_id_sce |> 
      purrr::transpose() |> 
      lapply(
        FUN = function(x) {
          
          .x = x[[2]]
          .y = x[[1]] |> str_remove("\\.h5ad")
          
          # Check if the 'sce' has only one cell (column)
          if(ncol(assay(.x)) == 1) {
            
            # Duplicate the assay to prevent saving errors due to single-column matrices
            my_assay = cbind(assay(.x), assay(.x))
            # Rename the second column to distinguish it
            colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            
            cd = colData(.x)
            cd = cd |> rbind(cd)
            rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
            
            
            
            .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
          } 
          
          
          # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
          .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
          
          # CALCULATE CPM
          .x =  SingleCellExperiment(assay = list( cpm = calculateCPM(.x, assay.type = names(assays(.x))[1])), colData = colData(.x)) 
          
          # My attempt to save a integer, sparse, delayed matrix (with zellkonverter it is not possible to save integers)
          # .x |> assay() |> type() <- "integer"
          # .x |> saveHDF5SummarizedExperiment("~/temp", as.sparse = T, replace = T)
          
          # Save the experiment data to the specified counts cache directory
          .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
          
          return(TRUE)  # Indicate successful saving
        }
        #,
        #BPPARAM = bp  # Use the defined parallel backend
      )
    
    return("saved")
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata_cpm <- purrr::insistently(save_anndata_cpm, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  # Function to process matrix in vertical slices
  process_matrix_in_slices <- function(h5_matrix, output_filepath, output_filepath_temp, chunk_size = 1000) {
    # Load the HDF5 matrix
    n_rows <- dim(h5_matrix)[1]
    n_cols <- dim(h5_matrix)[2]
    
    if (file.exists(output_filepath)) {
      file.remove(output_filepath)
      cat("Existing output file removed.\n")
    }
    if (file.exists(output_filepath_temp)) {
      file.remove(output_filepath_temp)
      cat("Existing output file removed.\n")
    }
    
    # Create an empty list to hold the slices
    slice_list <- list()
    
    # Loop through the matrix in chunks
    for (start_col in seq(1, n_cols, by = chunk_size)) {
      end_col <- min(start_col + chunk_size - 1, n_cols)
      cat("Processing columns", start_col, "to", end_col, "\n")
      
      # Extract a slice of the matrix
      matrix_slice <- as.matrix(h5_matrix[, start_col:end_col, drop=FALSE])
      
      # Calculate ranks for the slice
      ranked_slice <- singscore::rankGenes(matrix_slice)  %>% `-` (1) 
      
      # Convert the ranked slice to sparse format
      sparse_ranked_slice <- as(ranked_slice, "CsparseMatrix")
      
      # Write the slice to the output HDF5 file
      HDF5Array::writeHDF5Array(
        sparse_ranked_slice,
        filepath = output_filepath_temp,
        name = paste0("rank_", start_col, "_to_", end_col),
        as.sparse = TRUE,
        H5type = "H5T_STD_I32LE"
      ) 
      
      # Store the slice name for later binding
      slice_list[[length(slice_list) + 1]] <- paste0("rank_", start_col, "_to_", end_col)
    }
    
    
    slice_list |> map(~HDF5Array::HDF5Array(output_filepath_temp, name =.x)) |> do.call(cbind, args=_)
    
    # # Bind all slices into a single HDF5 dataset
    # for (i in seq_along(slice_list)) {
    #   slice <- HDF5Array::HDF5Array(output_filepath_temp, name = slice_list[[i]])
    #   if (i == 1) {
    #     final_matrix <- slice
    #   } else {
    #     final_matrix <- cbind(final_matrix, slice)
    #   }
    # }
    
    # # Save the final matrix back to the HDF5 file
    # result_matrix = 
    #   HDF5Array::writeHDF5Array(
    #     final_matrix,
    #     filepath = output_filepath,
    #     name = "final_ranked_matrix",
    #     as.sparse = TRUE,
    #     H5type = "H5T_STD_I32LE"
    #   )
    # 
    # file.remove(output_filepath_temp)
    # 
    # result_matrix
  }
  
  save_rank_per_cell = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
    
    # # Parallelise
    # cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    # bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    # 
    
    .x = dataset_id_sce |> pull(sce) |> _[[1]]
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    #.y = dataset_id_sce |> pull(file_id_cellNexus_pseudobulk) |> _[[1]] |> str_remove("\\.h5ad")
    
    # Check if the 'sce' has only one cell (column)
    if(ncol(assay(.x)) == 1) {
      
      # Duplicate the assay to prevent saving errors due to single-column matrices
      my_assay = cbind(assay(.x), assay(.x))
      # Rename the second column to distinguish it
      colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
      
      cd = colData(.x)
      cd = cd |> rbind(cd)
      rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
      
      
      
      .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
    } 
    
    
    # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
    .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1], subset_up_to_number_of_cells = 100)
    
    print("start ranking")
    
    # CALCULATE rank
    rank_assay = 
      .x |>
      assay() |> 
      
      # This because some datasets are still > 1M cells
      process_matrix_in_slices(
        paste(c(cache_directory, "/", .y, "_rank_matrix.HDF5Array"), collapse = ""), 
        paste(c(cache_directory, "/", .y, "_rank_matrix_temp.HDF5Array"), collapse = ""), 
        chunk_size = 1000
      )
    
    print("creating SCE")
    
    .x =  SingleCellExperiment(assay = list( rank = rank_assay), colData = colData(.x)) 
    
    print("saving")
    
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    # Delete the temp file
    file.remove(paste(c(cache_directory, "/", .y, "_rank_matrix_temp.HDF5Array"), collapse = ""))
    
    return(TRUE)  # Indicate successful saving
    
    
    
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_rank_per_cell <- purrr::insistently(save_rank_per_cell, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  save_anndata_sct = function(dataset_id_sce, cache_directory){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    

    .x = dataset_id_sce |> pull(sct) |> _[[1]]
    
    # Some sct have 0 cells after QC
    if (ncol(.x) == 0 || is.null(.x)) return(NULL)
    
    .y = dataset_id_sce |> pull(file_id_cellNexus_single_cell) |> _[[1]] |> str_remove("\\.h5ad")
    
    .x |> assays() |> names() = "sct"
    
    # Save the experiment data to the specified counts cache directory
    .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
    
    return(TRUE)  # Indicate successful saving
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_save_anndata_sct <- purrr::insistently(save_anndata_sct, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  cbind_sce_by_dataset_id = function(target_name_grouped_by_dataset_id, file_id_db_file, cell_id_dict, my_store){
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    my_file_id = unique(target_name_grouped_by_dataset_id$file_id_cellNexus_single_cell) 
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      filter(dataset_id == my_dataset_id) |>
      select(cell_id, sample_id, dataset_id, file_id_cellNexus_single_cell) 
    # |> 
    #   
    #   # Drop extension because it is added later 
    #   mutate(file_id_cellNexus_single_cell = file_id_cellNexus_single_cell |> str_remove("\\.h5ad")) |> 
    #   as_tibble()
    
    file_id_db = 
      target_name_grouped_by_dataset_id |> 
      left_join(file_id_db, copy = TRUE)
    
    
    dataset_cell_dict = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{cell_id_dict}')"))
      )  |> 
      filter(file_id_cellNexus_single_cell == my_file_id)
    # |> 
    #   filter(dataset_id == my_dataset_id) |> 
    #   select(cell_id, dataset_id, new_cell_id) 
    
    
    file_id_db = 
      file_id_db |> 
      left_join(dataset_cell_dict, by = c("file_id_cellNexus_single_cell", "cell_id" ), copy=T  )
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sce_df = 
      file_id_db |> 
      nest(cells = c(cell_id, new_cell_id)) |> 
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = bplapply(
          sce_target_name,
          FUN = function(x) tar_read_raw(x, store = my_store) |> 
            select(.cell, donor_id, dataset_id, sample_id, cell_type) |> 
            mutate(sample_id = as.factor(sample_id)) # lighter
            ,  # Read the raw SingleCellExperiment object
          BPPARAM = bp  # Use the defined parallel backend
        )) |> 
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null)) |> 
      mutate(sce = map2(sce, cells, ~ {
        
        cell_map <- setNames(.y$new_cell_id, .y$cell_id) 
        
        .x |> filter(.cell %in% names(cell_map)) %>% 
          {
            colnames(.) <- cell_map[colnames(.)]
            .
          } |>
          
          # TEMPORARY FIX. NEED TO INVESTIGATE WHY THE SUFFIX HAPPENS
          mutate(sample_id = stringr::str_replace(sample_id, ".h5ad$",""))
        
      }, .progress = TRUE))
  
      # mutate(
      #   # Because sct is calculated post QC, NA target_name cant find sample_id and dataset_id internally,
      #   #    meaning poor samples.
      #   sct = if (!is.na(sct_target_name)) {
      #     sct_data = bplapply(
      #       sct_target_name,
      #       FUN = function(x) tar_read_raw(x, store = my_store),  # Read the raw SingleCellExperiment object
      #       BPPARAM = bp  # Use the defined parallel backend
      #     )
      #     
      #     map2(sct_data, cells, ~ if (is.null(.x)) {NULL} else {
      #       .x |> filter(.cell %in% .y$cell_id) |>
      #         
      #         # TEMPORARY FIX. NEED TO INVESTIGATE WHY THE SUFFIX HAPPENS
      #         mutate(sample_id = stringr::str_replace(sample_id, ".h5ad$",""))},
      #       
      #       .progress = TRUE
      #     )
      #   } else {list(NULL)}
      # ) 
      # 
    if(nrow(sce_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
    # plan(multisession, workers = 20)
    sce_df |> 
      
      # # Step 4: Group the data by 'dataset_id' and 'tar_group' for further summarization
      # group_by(dataset_id, tar_group, chunk) |>
      # 
      
      # FORCEFULLY drop all but counts and metadata 
      # int_colData(.x) = DataFrame(row.names = colnames(.x))
      # Creates error
      # THIS SHOULD HAVE BEEN DONE IN THE TRANFORM HPCell
      mutate(sce = map(sce, ~  SingleCellExperiment(assay = assays(.x), colData = colData(.x)) )) |> 
      
      # Step 5: Combine all 'sce' objects within each group into a single 'sce' object
      group_by(file_id_cellNexus_single_cell) |> 
      summarise( sce =  list(do.call(cbind, args = sce) ),
                 # A steo to check missing cells 
                 cells = list(do.call(rbind, args = cells))) 
    
      # mutate(
      #   sct = map(sct, ~  if (!is.null(.x)) SingleCellExperiment(assay = assays(.x), colData = colData(.x)) )
      # ) |> 
      # group_by(file_id_cellNexus_single_cell) |> 
      # summarise( sct =  list(do.call(cbind, args = sct) ),
      #            
      #            # A steo to check missing cells 
      #            cells = list(do.call(rbind, args = cells))) 
    
    # mutate(sce = map(sce,
    #                  ~ { .x = 
    #                    .x  |> 
    #                    left_join(file_id_db, by = join_by(.cell==cell_id, dataset_id==dataset_id, sample_id==sample_id)) 
    #                  .x |> 
    #                    HPCell:::splitColData(colData(.x)$file_id_cellNexus_single_cell) |>  # Split 'sce' by 'cell_type'
    #                    enframe(name = "file_id_cellNexus_single_cell", value = "sce")  # Convert to tibble with 'cell_type' and 'sce' columns
    #                  })) |> 
    # Step 8: Unnest the list of 'sce' objects to have one row per 'cell_type'
    # unnest_single_cell_experiment(sce) 
    
    
  }
  
  cbind_sct_by_dataset_id = function(target_name_grouped_by_dataset_id, file_id_db_file, cell_id_dict, my_store){
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    my_file_id = unique(target_name_grouped_by_dataset_id$file_id_cellNexus_single_cell) 
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      filter(dataset_id == my_dataset_id) |>
      select(cell_id, sample_id, dataset_id, file_id_cellNexus_single_cell) 
    # |> 
    #   
    #   # Drop extension because it is added later 
    #   mutate(file_id_cellNexus_single_cell = file_id_cellNexus_single_cell |> str_remove("\\.h5ad")) |> 
    #   as_tibble()
    
    file_id_db = 
      target_name_grouped_by_dataset_id |> 
      left_join(file_id_db, copy = TRUE)
    
    
    dataset_cell_dict = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{cell_id_dict}')"))
      )  |> 
      filter(file_id_cellNexus_single_cell == my_file_id)
    # |> 
    #   filter(dataset_id == my_dataset_id) |> 
    #   select(cell_id, dataset_id, new_cell_id) 
    
    
    file_id_db = 
      file_id_db |> 
      left_join(dataset_cell_dict, by = c("file_id_cellNexus_single_cell", "cell_id"), copy=T  )
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sct_df = file_id_db |> 
      nest(cells = c(cell_id, new_cell_id)) %>% 
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sct = bplapply(
          sct_target_name,
          FUN = function(x) {
            if (is.na(x)) {
              return(NULL) # because cant get sample_id and dataset_id from NULL sct_matrix
            }
            
          tar_read_raw(x, store = my_store) |>
            select(.cell, donor_id, dataset_id, sample_id, cell_type) |>
            mutate(sample_id = as.factor(sample_id))
          
          },  # Read the raw SingleCellExperiment object
          BPPARAM = bp  # Use the defined parallel backend
        )) |> 
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sct, is.null)) |> 
      mutate(sct = map2(sct, cells, ~ {
        
        cell_map <- setNames(.y$new_cell_id, .y$cell_id) 
        
        .x |> filter(.cell %in% names(cell_map)) %>% 
          {
            colnames(.) <- cell_map[colnames(.)]
            .
          } |>
          
          # TEMPORARY FIX. NEED TO INVESTIGATE WHY THE SUFFIX HAPPENS
          mutate(sample_id = stringr::str_replace(sample_id, ".h5ad$",""))
        
      }, .progress = TRUE))

    if(nrow(sct_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
    sct_df |>
      mutate(
        sct = map(sct, \(x) {
          if (is.null(x)) return(NULL)
          SingleCellExperiment(assays = assays(x), colData = colData(x))
        })
      ) |>
      group_by(file_id_cellNexus_single_cell) |>
      summarise(
        sct = {
          scts <- compact(sct)  # drop NULLs inside each group
          
          list(
            if (length(scts) == 0) {
              NULL
            } else {
              
              # A few big samples do not return all features because it reached R limit 2^31-1 in SCTransform
              common_genes <- cellNexus:::check_gene_overlap(scts)
              
              # subset to intersection genes (and keep same order across objects)
              scts2 <- map(scts, \(z) z[common_genes, , drop = FALSE])
              
              do.call(SummarizedExperiment::cbind, scts2)
            }
          )
        },
        cells = list(do.call(rbind, cells)),
        .groups = "drop"
      )
    # mutate(
    #   sct = map(sct, ~  if (!is.null(.x)) SingleCellExperiment(assay = assays(.x), colData = colData(.x)) )
    # ) |> 
    # group_by(file_id_cellNexus_single_cell) |> 
    # summarise( sct =  list(do.call(cbind, args = sct) ),
    #            
    #            # A steo to check missing cells 
    #            cells = list(do.call(rbind, args = cells))) 
    # mutate(sce = map(sce,
    #                  ~ { .x = 
    #                    .x  |> 
    #                    left_join(file_id_db, by = join_by(.cell==cell_id, dataset_id==dataset_id, sample_id==sample_id)) 
    #                  .x |> 
    #                    HPCell:::splitColData(colData(.x)$file_id_cellNexus_single_cell) |>  # Split 'sce' by 'cell_type'
    #                    enframe(name = "file_id_cellNexus_single_cell", value = "sce")  # Convert to tibble with 'cell_type' and 'sce' columns
    #                  })) |> 
    # Step 8: Unnest the list of 'sce' objects to have one row per 'cell_type'
    # unnest_single_cell_experiment(sce) 
    
    
  }
  
  # Because they have an inconsistent failure. If I start the pipeline again they might work. Strange.
  insistent_cbind_sct_by_dataset_id <- purrr::insistently(cbind_sct_by_dataset_id, rate = purrr::rate_delay(pause = 60, max_times = 3), quiet = FALSE)
  
  
  get_dataset_id = function(target_name, my_store){
    # Try reading the target safely (for some failing targets)
    sce = tryCatch(
      tar_read_raw(target_name, store = my_store),
      error = function(e) return(NULL)
    )
    # sce = tar_read_raw(target_name, store = my_store)
    
    # Still need to catch target_name
    if(sce |> is.null()) return(tibble(sample_id = NA_character_, 
                                       dataset_id= NA_character_, 
                                       target_name= !!target_name))
    
    sce |> 
      
      # TEMPORARY FIX. NEED TO INVESTIGATE WHY THE SUFFIX HAPPENS
      mutate(sample_id = stringr::str_replace(sample_id, ".h5ad$","")) |> 
      
      distinct(sample_id, dataset_id) |> mutate(target_name = !!target_name)
  }
  
  create_chunks_for_reading_and_saving = function(dataset_id_sample_id, cell_metadata){
    
    # Solve sample_id mismatches because some end with .h5ad suffix while others dont 
    dataset_id_sample_id |> 
      
      # TEMPORARY FIX. NEED TO INVESTIGATE WHY THE SUFFIX HAPPENS
      mutate(sample_id = stringr::str_replace(sample_id, ".h5ad$", "")) |>
      
      left_join(
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))
        )   |> 
          distinct(dataset_id, sample_id, sample_chunk, cell_chunk, file_id_cellNexus_single_cell) |> 
          as_tibble(), 
        copy=T
      )
  }
  
  
  cbind_sce_by_dataset_id_get_missing_cells = function(dataset_id_sce){
    
    dataset_id_sce |>
      mutate(
        missing_cells = map2(
          sce, 
          cells, 
          ~{
            cells_in_sce <- .x |> colnames() |> sort()
            
            cells_in_query <- .y$new_cell_id |> unique() |> sort()
            
            # Find differences
            tibble(cell_id = setdiff(cells_in_query, cells_in_sce))
          }
        )
      ) |> 
      select(file_id_cellNexus_single_cell, missing_cells)
    
  }
  
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/cellNexus/2024-07-01/process_samples_hpcell_target_store", deployment = "main"), # MODIFY HERE: HPCell targets store to read SCEs from
    tar_target(cache_directory, "/vast/scratch/users/shen.m/cellNexus/cellxgene/01-07-2024", deployment = "main"), # MODIFY HERE: output cache directory for saved anndata files
    # This is the store for retrieving missing cells between cellnexus metadata and sce. A different store as it was done separately
    #tar_target(cache_directory, "/vast/scratch/users/shen.m/debug2/cellxgene/19-12-2024", deployment = "main"),
    tar_target(
      cell_metadata,
      "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_3_2_mengyuan.parquet", # MODIFY HERE: final metadata parquet (should match the COPY TO output above)
      packages = c( "arrow","dplyr","duckdb")
      
    ),
    
    tar_target(
      cell_id_dict,
      "/vast/projects//cellxgene_curated/metadata_cellxgene_mengyuan/file_id_cell_id_dict_v1_0_0_Jul_2024.parquet", # MODIFY HERE: cell_id dictionary parquet
      packages = c( "arrow","dplyr","duckdb")
      
    ),
   
    
    # pre-calculated counts
    tar_target(
      target_name,
      tar_meta(
        starts_with("sce_transformed_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    tar_target(
      dataset_id_sample_id,
      get_dataset_id(target_name, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic")
      )
    ),
    
    # pre-calculated sct
    tar_target(
      sct_target_name,
      tar_meta(
        starts_with("sct_matrix_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    tar_target(
      sct_dataset_id_sample_id,
      get_dataset_id(sct_target_name, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(sct_target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic")
      )
    ),
    
    # join
    tar_target(
      dataset_id_sample_id_target_names,
      dataset_id_sample_id |> left_join(sct_dataset_id_sample_id, by = c("sample_id", "dataset_id"), copy=T) |>
        dplyr::rename(sce_target_name = target_name.x, 
                      sct_target_name = target_name.y),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic")
      )
    ),
    
    tar_target(
      target_name_grouped_by_dataset_id,
      create_chunks_for_reading_and_saving(dataset_id_sample_id_target_names, cell_metadata) |> 
        
        # # FOR TESTING PURPOSE ONLY
        # filter(file_id_cellNexus_single_cell %in% c("3cef5b6aa0f5772485bb710f71e69456___1.h5ad",
        #                                             "cd2caa6de850f73af4ca78a2ea307dd4___1.h5ad")) |>
        
        group_by(dataset_id, sample_chunk, cell_chunk, file_id_cellNexus_single_cell) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic")
      ), 
      packages = c("arrow", "duckdb", "dplyr", "glue", "targets")
      
    ),
    
    tar_target(
      dataset_id_sce,
      cbind_sce_by_dataset_id(target_name_grouped_by_dataset_id, cell_metadata, cell_id_dict, my_store = my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),
    
    tar_target(
      dataset_id_sct,
      cbind_sct_by_dataset_id(target_name_grouped_by_dataset_id, cell_metadata, cell_id_dict, my_store = my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),
    
    # This target was run for retrieving missing cells analysis only
    tar_target(
      missing_cells_tbl,
      cbind_sce_by_dataset_id_get_missing_cells(dataset_id_sce),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly", "purrr"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),
    
    
    tar_target(
      save_anndata,
      insistent_save_anndata(dataset_id_sce, paste0(cache_directory, "/counts")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),

    tar_target(
      saved_dataset_cpm,
      insistent_save_anndata_cpm(dataset_id_sce, paste0(cache_directory, "/cpm")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),

    tar_target(
      saved_dataset_rank,
      insistent_save_rank_per_cell(dataset_id_sce, paste0(cache_directory, "/rank")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),
    
    tar_target(
      saved_sct,
      insistent_save_anndata_sct(dataset_id_sct, paste0(cache_directory, "/sct")),
      pattern = map(dataset_id_sct),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    )
  )
  
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(store_file_cellNexus, "_target_script.R"), 
    store = store_file_cellNexus, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

missing_cells_tbl = tar_read(missing_cells_tbl, store = store_file_cellNexus)
missing_cells_tbl <- map(missing_cells_tbl$missing_cells, ~ {.x}) |> bind_rows()
missing_cells <- missing_cells_tbl |> pull(cell_id)

cell_metadata |> filter(!cell_id %in% missing_cells) |> 
  
  # This method of save parquet to parquet is faster 
  cellNexus:::duckdb_write_parquet(path = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_3_2_filtered_missing_cells_mengyuan.parquet") # MODIFY HERE: output parquet after filtering missing cells

