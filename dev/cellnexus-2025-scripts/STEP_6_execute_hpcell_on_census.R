# Step 2
library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(stringr)
library(HPCell)
library(arrow)
library(CuratedAtlasQueryR)
library(targets)
library(crew)
library(crew.cluster)
library(duckdb)

version <- "2025-11-08"
directory = glue::glue("/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/{version}/")

downloaded_samples_tbl <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/2025-11-08_census_samples_to_download.parquet')")
  ) |>
  collect() |> 
  mutate(file_name = file.path(directory, paste0(sample_id, ".h5ad")) |> as.character())

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/"

sample_meta <- tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset.parquet')")
) |> collect()
  
sample_tbl = downloaded_samples_tbl |> left_join(CuratedAtlasQueryR::get_metadata(cache_directory = tempdir()) |> dplyr::select(dataset_id, contains("norm")) |>
                                                   distinct() |> dplyr::filter(!is.na(x_normalization)) |>
                                                   as_tibble(), by = "dataset_id")

sample_tbl <- sample_tbl |> left_join(sample_meta |> distinct(dataset_id, assay, x_approximate_distribution), 
                                      by = c("dataset_id", "assay")) |> 
  distinct(file_name, dataset_id, sample_id, assay, x_normalization, x_approximate_distribution) |>
  mutate(transform_method = case_when(str_like(x_normalization, "C%") ~ "log",
                                      x_normalization == "none" ~ "log",
                                      x_normalization == "normalized" ~ "log",
                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "log",
                                      is.na(x_normalization) & x_approximate_distribution == "NORMAL" ~ "NORMAL",
                                      is.na(x_normalization) & x_approximate_distribution == "COUNT" ~ "COUNT",
                                      str_like(x_normalization, "%canpy%") ~ "log1p",
                                      TRUE ~ x_normalization)) |>
  
  mutate(method_to_apply =  case_when(transform_method %in% c("log","LogNormalization","LogNormalize","log-normalization") ~ "exp",
                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "exp",
                                      str_like(transform_method, "Counts%") ~ "exp",
                                      str_like(transform_method, "%log2%") ~ "exp",
                                      transform_method %in% c("log1p", "log1p, base e", "Scanpy",
                                                              "scanpy.api.pp.normalize_per_cell method, scaling factor 10000") ~ "expm1",
                                      transform_method == "log1p, base 2" ~ "expm1",
                                      transform_method == "NORMAL" ~ "exp",
                                      transform_method == "COUNT" ~ "identity",
                                      is.na(transform_method) ~ "identity"
  ) ) |>
  mutate(comment = case_when(str_like(x_normalization, "Counts%")  ~ "a checkpoint for max value of Assay must <= 50",
                             is.na(x_normalization) & is.na(x_approximate_distribution) ~ "round negative value to 0",
                             x_normalization == "normalized" ~ "round negative value to 0"
  ))

sample_tbl <- saveRDS("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/sample_tbl.rds")


sample_tbl <- readRDS("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/sample_tbl.rds")
sample_summary_df <- arrow::read_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/sample_distribution_summary.parquet")
#sample_tbl = sample_tbl |> head()
# Infer missing x_approximate_distribution 
impute_x_approximate_distribution <- function(df) {
  df |> mutate(
    inferred_distribution = case_when(
      !is.na(x_approximate_distribution) ~ x_approximate_distribution,
      
      # 1) No negatives, no large values, no integers, no floating
      is.na(x_approximate_distribution) & !has_negative & !max_gt_20 & !all_integer & !has_floating ~ "log1p",
      
      # 2) No negatives, has large values
      is.na(x_approximate_distribution) & !has_negative &  max_gt_20 & !all_integer & !has_floating  ~ "raw",
      
      # 3) No negatives, large values, all integer, has floating
      is.na(x_approximate_distribution) & !has_negative &  max_gt_20 & all_integer & !has_floating ~ "raw",
      
      # 4) Has negatives, no large values, no integer, no floating
      is.na(x_approximate_distribution) & has_negative & !max_gt_20 & !all_integer & !has_floating ~ "log",
      
      # 5) Has negatives and large values
      is.na(x_approximate_distribution) & has_negative &  max_gt_20 & !all_integer & !has_floating ~ "raw"
    )
  )
}

sample_tbl = sample_tbl |> left_join(sample_summary_df |> mutate(sample_id = str_remove(sample_id, ".h5ad"))) |> 
  impute_x_approximate_distribution() |> 
  mutate(count_upper_bound = case_when(
    is.na(x_approximate_distribution) & has_negative & !max_gt_20 & !all_integer & !has_floating ~ 10,
    TRUE ~ 20
  )) |> 
  # base our filtering on % of expressed genes for panel technologies 500/20K : ?/462
  mutate(feature_thresh = ifelse(assay == "BD Rhapsody Targeted mRNA", 11, 200)) |>
  
  # Inverse distribution
  mutate(method_to_apply = case_when(inferred_distribution == "log" ~ "exp",
                                        inferred_distribution == "log1p" ~ "expm1",
                                        inferred_distribution == "raw" ~ "identity"))


# Datasets missing colData are excluded.
dataset_to_exclude <- sample_tbl |> dplyr::filter(inferred_distribution |> is.na(), 
                                                  is.na(x_approximate_distribution),  
                                                  is.na(has_negative),  is.na(max_gt_20),  
                                                  is.na(all_integer), is.na(has_floating)) |>
  pull(dataset_id) |> unique()

# setwd(glue("{store}"))
sliced_sample_tbl = 
  sample_tbl |> 
  filter(!dataset_id %in% dataset_to_exclude) |> 
  dplyr::select(file_name, dataset_id, sample_id, method_to_apply, assay, feature_thresh, count_upper_bound, list_length)


sliced_sample_tbl <- saveRDS("/vast/scratch/users/shen.m/cellNexus_run/sliced_sample_tbl.rds")

sliced_sample_tbl <- readRDS("/vast/scratch/users/shen.m/cellNexus_run/sliced_sample_tbl.rds")
# sliced_sample_tbl = sliced_sample_tbl |> filter(list_length < 1000) |> head()
  #head(2) |> bind_rows(sliced_sample_tbl |> filter(count_upper_bound == 10) |> head(2))

# Enable sample_names.rds to store sample names for the input
sample_names <-
  sliced_sample_tbl |> 
  pull(file_name) |> 
  str_replace("/home/users/allstaff/shen.m/scratch", "/vast/scratch/users/shen.m") |> 
  set_names(sliced_sample_tbl |> pull(sample_id))
functions = sliced_sample_tbl |> pull(method_to_apply)
feature_thresh = sliced_sample_tbl |> pull(feature_thresh)
count_upper_bound = sliced_sample_tbl |> pull(count_upper_bound)

# sample_names <- saveRDS("/vast/scratch/users/shen.m/cellNexus_run/sample_names_filtered_by_mengyuan_apr_2024.rds")
# functions <- saveRDS("/vast/scratch/users/shen.m/cellNexus_run/functions.rds")
# feature_thresh <- saveRDS("/vast/scratch/users/shen.m/cellNexus_run/feature_thresh.rds")
# count_upper_bound <- saveRDS("/vast/scratch/users/shen.m/cellNexus_run/count_upper_bound.rds")

my_store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")
job::job({
  
  library(HPCell)
  
  sample_names |>
    initialise_hpc(
      store = my_store,
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      # tier = tiers, # WE DON"T NEED AS WE HAVE ELASTIC RESOURCES NOW
      computing_resources = list(
        
        crew.cluster::crew_controller_slurm(
          name = "elastic",
          workers = 300,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_error = 10,
          options_cluster = crew.cluster::crew_options_slurm(
            memory_gigabytes_required = c(60, 80, 100, 120, 150, 200), 
            cpus_per_task = c(2, 2, 5, 6, 7, 10),
            time_minutes = c(60*24, 60*24, 60*24, 60*24, 60*24,60*24),
            verbose = T
          )
        )
        
      ),
      verbosity = "summary",
      #update = "never", 
      update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = functions, target_output = "sce_transformed", scale_max = count_upper_bound) |>
    
    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = feature_thresh) |>
    
    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>

    # Cell type harmonisation
    celltype_consensus_constructor(target_input = "sce_transformed",
                                   target_output = "cell_type_concensus_tbl") |>

    # Alive identification
    remove_dead_scuttle(target_input = "sce_transformed", target_annotation = "cell_type_concensus_tbl",
                        group_by = "cell_type_unified_ensemble") |>

    # Doublets identification
    remove_doublets_scDblFinder(target_input = "sce_transformed") |>

    # Pseudobulk
    calculate_pseudobulk(target_input = "sce_transformed",
                         group_by = "cell_type_unified_ensemble") |>

    # metacell
    cluster_metacell(target_input = "sce_transformed",  group_by = "cell_type_unified_ensemble") |>
    # 
    # # Cell Chat
    # ligand_receptor_cellchat(target_input = "sce_transformed",
    #                          group_by = "cell_type_unified_ensemble") |>
    
    print()
  
  
})
#tar_invalidate(names = metacell_tbl, store = my_store )
tar_meta(store = my_store) |> filter(!is.na(error)) |>  arrange(desc(time)) |> View()
tar_meta(store = my_store) |> filter(!is.na(error)) |> distinct(name, error)
tar_meta(starts_with("annotation_tbl_"), store = "/vast/scratch/users/shen.m/cellNexus_target_store") |> 
  filter(!data |> is.na()) |> arrange(desc(time)) |> select(error, name)

# Debug cellchat
tar_workspace(pseudobulk_se_iterated_3b741ea964ed386d, store = my_store, 
              script = paste0(my_store,".R"))
debugonce(create_pseudobulk)
create_pseudobulk(sce_transformed, sample_names, empty_tbl, alive_tbl, NULL, annotation_tbl, 
                  cell_type_concensus_tbl, doublet_tbl, "cell_type_unified_ensemble", "~/scratch/cache_temp",
                  NULL, data_container_type)


alive_identification(sce_transformed,empty_droplets_tbl = empty_tbl, cell_type_column = "cell_type_unified_ensemble", feature_nomenclature = gene_nomenclature)
# cell_communication(sce_transformed, empty_droplets_tbl = empty_tbl, alive_identification_tbl = alive_tbl,
#                    doublet_identification_tbl = doublet_tbl, cell_type_tbl = cell_type_concensus_tbl,
#                    cell_type_column = "cell_type_unified_ensemble",
#                    feature_nomenclature = gene_nomenclature)


#' Pipeline for Lightening Annotations in High-Performance Computing Environment
#' 
#' This pipeline is designed to read, process, and "lighten" large annotation tables in an HPC environment.
#' It uses the `targets` package for reproducibility and `crew` for efficient job scheduling on a Slurm cluster.
#' The `lighten_annotation` function selects and processes specific columns from large tables to reduce memory usage.
#' 
#' The pipeline consists of:
#' - **Crew Controllers**: Four tiers of Slurm controllers with varying memory allocations to optimize resource usage.
#' - **Targets**:
#'   - `my_store`: Defines the path to the target storage directory, ensuring all targets use the correct storage location.
#'   - `target_name`: Retrieves metadata to identify branch targets for annotation.
#'   - `annotation_tbl_light`: Applies `lighten_annotation` to process each target name, optimally running with `tier_1` resources.
#' 
#' @libraries:
#'   - `dplyr`, `magrittr`, `tibble`, `targets`, `tarchetypes` for data manipulation and pipeline structure.
#'   - `crew`, `crew.cluster` for parallel computation and cluster scheduling in an HPC environment.
#' 
#' @options:
#'   - Memory settings, garbage collection frequency, and error handling are set to handle large data efficiently.
#'   - The `cue` option is set to `never` for forced target updates if needed.
#'   - `controller` is a group of Slurm controllers to manage computation across memory tiers.
#' 
#' @function `lighten_annotation`: Processes each annotation table target, unnesting and selecting specific columns to reduce data size.
#'
#' @example Usage:
#'   The pipeline script is saved as `/vast/scratch/users/shen.m/lighten_annotation_tbl_target.R` by tar_script and can be run using `tar_make()`.
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
    #debug = "annotation_tbl_light", 
    cue = tar_cue(mode = "never"), 
    controller = crew_controller_group(
      list(
        crew_controller_slurm(
          name = "tier_1", 
          script_lines = "#SBATCH --mem 8G",
          slurm_cpus_per_task = 1, 
          workers = 200, 
          tasks_max = 10,
          verbose = T,
          #launch_max = 5, 
          seconds_idle = 30,
          slurm_time_minutes = 480
        ),
        
        crew_controller_slurm(
          name = "tier_2",
          script_lines = "#SBATCH --mem 10G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          #launch_max = 5, 
          seconds_idle = 30,
          slurm_time_minutes = 480
        ),
        crew_controller_slurm(
          name = "tier_3",
          script_lines = "#SBATCH --mem 15G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          #launch_max = 5, 
          seconds_idle = 30,
          slurm_time_minutes = 480
        ),
        crew_controller_slurm(
          name = "tier_4",
          script_lines = "#SBATCH --mem 50G",
          slurm_cpus_per_task = 1,
          workers = 30,
          tasks_max = 10,
          verbose = T,
          #launch_max = 5, 
          seconds_idle = 30,
          slurm_time_minutes = 480
        )
      )
    ), 
    trust_object_timestamps = TRUE
    # workspaces = "annotation_tbl_light_ffcd3d5a64bedf1f"
  )
  
  lighten_annotation = function(target_name, my_store ){
    annotation_tbl = tar_read_raw( target_name,  store = my_store )
    if(annotation_tbl |> is.null()) { 
      warning("this annotation is null -> ", target_name)
      return(NULL) 
    }
    
    annotation_tbl |> 
      unnest(blueprint_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
      unnest(monaco_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |> 
      rename(cell_ = .cell)
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/scratch/users/shen.m/cellNexus_target_store_2025-11-08", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("annotation_tbl_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    )    ,
    
    tar_target(
      annotation_tbl_light,
      lighten_annotation(target_name, my_store),
      packages = c("dplyr", "tidyr"),
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_1")
      )
    )
  )
  
  
}, script = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}.R"),
    store = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}"), 
    reporter = "summary"
  )
  
})

# Sample metadata
library(arrow)
library(dplyr)
library(duckdb)
library(targets)
library(cellNexus)

cellNexus:::duckdb_write_parquet = function(data_tbl, output_parquet, compression = "gzip") {
  
  # Establish connection to DuckDB in-memory database
  con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Register `data_tbl` within the DuckDB connection (this doesn't load it into memory)
  duckdb::duckdb_register(con_write, "data_tbl_view", data_tbl)
  
  # Use DuckDB's COPY command to write `data_tbl` directly to Parquet with compression
  copy_query <- paste0("
  COPY data_tbl_view TO '", output_parquet, "' (FORMAT PARQUET, COMPRESSION '", compression, "');
  ")
  
  # Execute the COPY command
  dbExecute(con_write, copy_query)
  
  # Unregister the temporary view
  duckdb::duckdb_unregister(con_write, "data_tbl_view")
  
  # Disconnect from the database
  dbDisconnect(con_write, shutdown = TRUE)
}


# Write annotation light
cell_metadata <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet')")
  ) |>
  mutate(cell_ = paste0(cell_, "___", dataset_id)) |> 
  select(cell_, observation_joinid, contains("cell_type"), dataset_id,  self_reported_ethnicity, tissue, donor_id,  sample_id, is_primary_data, assay)


cell_annotation = 
  tar_read(annotation_tbl_light, store = glue::glue("/vast/scratch/users/shen.m/lighten_annotation_tbl_target_{version}")) |> 
  dplyr::rename(
    blueprint_first_labels_fine = blueprint_first.labels.fine, 
    monaco_first_labels_fine = monaco_first.labels.fine, 
    azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
  ) 

cell_annotation = cell_annotation |> mutate(
  blueprint_first_labels_fine = ifelse(is.na(blueprint_first_labels_fine), "Other", blueprint_first_labels_fine),
  monaco_first_labels_fine = ifelse(is.na(monaco_first_labels_fine), "Other", monaco_first_labels_fine),
  azimuth_predicted_celltype_l2=ifelse(is.na(azimuth_predicted_celltype_l2), "Other", azimuth_predicted_celltype_l2))

# cell_annotation |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/annotation_tbl_light.parquet",
#                                         compression = "zstd")
cell_annotation |> arrow::write_parquet("~/scratch/cache_temp/annotation_tbl_light.parquet",
                                        compression = "zstd")

empty_droplet = 
  tar_read(empty_tbl, store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")) |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

alive_cells = 
  tar_read(alive_tbl, store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")) |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

doublet_cells =
  tar_read(doublet_tbl, store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")) |>
  bind_rows() |>
  dplyr::rename(cell_ = .cell)

metacell = 
  tar_read(metacell_tbl, store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")) |> 
  bind_rows() |> 
  dplyr::rename(cell_ = cell) |> 
  dplyr::rename_with(
    ~ stringr::str_replace(.x, "^gamma", "metacell_"),
    starts_with("gamma")
  )

# Save cell type concensus tbl from HPCell output to disk
cell_type_concensus_tbl = tar_read(cell_type_concensus_tbl, store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")) |>  
  bind_rows() |> 
  dplyr::rename(cell_ = .cell)

cell_type_concensus_tbl = cell_type_concensus_tbl |> mutate(cell_type_unified_ensemble = 
                                                              ifelse(is.na(cell_type_unified_ensemble),
                                                                     "Unknown",
                                                                     cell_type_unified_ensemble)) 

# cell_type_concensus_tbl |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_type_concensus_tbl_from_hpcell.parquet",
#                                                 compression = "zstd")


cell_type_concensus_tbl |> arrow::write_parquet("~/scratch/cache_temp/cell_type_concensus_tbl_from_hpcell.parquet",
                                                compression = "zstd")
# This command needs a big memory machine
cell_metadata_joined = cell_metadata |> 
  left_join(empty_droplet, copy=TRUE) |>  
  left_join(cell_type_concensus_tbl, copy=TRUE) |>
  #left_join(cell_annotation, copy=TRUE) |>  
  left_join(alive_cells, copy=TRUE) |> 
  left_join(doublet_cells, copy=TRUE) |>
  left_join(metacell, copy=TRUE)

cell_metadata_joined |> filter(is.na(blueprint_first_labels_fine))

cell_metadata_joined2 = cell_metadata_joined |> as_tibble() |> 
  # Match to how pseudobulk annotations get parsed in HPCell/R/functions preprocessing_output()
  mutate(cell_type_unified_ensemble = ifelse(cell_type_unified_ensemble |> is.na(), "Unknown", cell_type_unified_ensemble)) |>
  mutate(data_driven_ensemble = ifelse(data_driven_ensemble |> is.na(), "Unknown", data_driven_ensemble))   |>
  mutate(blueprint_first_labels_fine = ifelse(blueprint_first_labels_fine |> is.na(), "Other", blueprint_first_labels_fine)) |> 
  mutate(monaco_first_labels_fine = ifelse(monaco_first_labels_fine |> is.na(), "Other", monaco_first_labels_fine)) |> 
  mutate(azimuth_predicted_celltype_l2 = ifelse(azimuth_predicted_celltype_l2 |> is.na(), "Other", azimuth_predicted_celltype_l2)) |> 
  mutate(azimuth = ifelse(azimuth |> is.na(), "Other", azimuth)) |> 
  mutate(blueprint = ifelse(blueprint |> is.na(), "Other", blueprint)) |> 
  mutate(monaco = ifelse(monaco |> is.na(), "Other", monaco))

# (!!!!) WHENEVER MAKE CHANGES, ALWAYS INCLUDE NEW Rhapsody DATA. Rhapsody TECH WERE UPDATED SEPARATELY in 2_1_rerun_hpcell_for_Rhapsody_targeted_panel.R. 
# cell_metadata_joined2 |>
#   arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_annotation.parquet",
#                        compression = "zstd")

cell_metadata_joined2 |>
  arrow::write_parquet("~/scratch/cache_temp/cell_annotation.parquet",
                       compression = "zstd")

# Cellchat output
ligand_receptor_tbl = tar_read(ligand_receptor_tbl, store = glue::glue("/vast/scratch/users/shen.m/cellNexus_target_store_{version}")) |> bind_rows()
# # save
# con <- dbConnect(duckdb::duckdb(), dbdir = "~/cellxgene_curated/metadata_cellxgene_mengyuan/cellNexus_lr_signaling_pathway_strength.duckdb")
# duckdb::dbWriteTable(con, "lr_pathway_table", ligand_receptor_tbl, overwrite = TRUE)
# dbDisconnect(con)



write_parquet_to_parquet = function(data_tbl, output_parquet, compression = "gzip") {
  
  # Establish connection to DuckDB in-memory database
  con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Register `data_tbl` within the DuckDB connection (this doesn't load it into memory)
  duckdb::duckdb_register(con_write, "data_tbl_view", data_tbl)
  
  # Use DuckDB's COPY command to write `data_tbl` directly to Parquet with compression
  copy_query <- paste0("
  COPY data_tbl_view TO '", output_parquet, "' (FORMAT PARQUET, COMPRESSION '", compression, "');
  ")
  
  # Execute the COPY command
  dbExecute(con_write, copy_query)
  
  # Unregister the temporary view
  duckdb::duckdb_unregister(con_write, "data_tbl_view")
  
  # Disconnect from the database
  dbDisconnect(con_write, shutdown = TRUE)
}
