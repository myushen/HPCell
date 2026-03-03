library(targets)
store = "/vast/scratch/users/shen.m/cellnexus_dataset_cell_map_Jul_2024_v1_2_0_target_store"
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
        # crew_controller_slurm(
        #   name = "elastic",
        #   workers = 300,
        #   tasks_max = 20,
        #   seconds_idle = 30,
        #   crashes_error = 10,
        #   options_cluster = crew_options_slurm(
        #     #memory_gigabytes_required = c(10, 20, 40, 80, 160), 
        #     memory_gigabytes_required = c(25, 35, 40, 80, 160), 
        #     cpus_per_task = c(2, 2, 5, 10, 20), 
        #     time_minutes = c(30, 30, 30, 60*4, 60*24),
        #     verbose = T
        #   )
        # )
        crew_controller_local(workers = 8)
      )
    ), 
    trust_object_timestamps = TRUE
  )
  
  get_unique_dataset <- function(cell_metadata){
    tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      distinct(dataset_id) |> pull()
  }
  
  create_dataset_cell_id_dict <- function(cell_metadata, datasetId) {
    
    
    tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      filter(dataset_id == datasetId) |>
      dbplyr::window_order(cell_id) |> 
      mutate(cell_index = row_number(),
             new_cell_id = paste(dataset_id, cell_index, sep = "___")) |> 
      select(cell_id, dataset_id, new_cell_id) |>
      collect()
  }
  
  list(
    tar_target(cell_metadata , "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cell_metadata_cell_type_consensus_v1_2_0_mengyuan.parquet",
               deployment = "main"),
    tar_target(
      unique_dataset,
      get_unique_dataset(cell_metadata) 
      # |> head(2)
      ,
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array")
      # resources = tar_resources(
      #   crew = tar_resources_crew(controller = "elastic")
      # )
    ),
    tar_target(
      cell_id_dict,
      create_dataset_cell_id_dict(cell_metadata, unique_dataset),
      pattern = map(unique_dataset),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly", "HDF5Array")
      # resources = tar_resources(
      #   crew = tar_resources_crew(controller = "elastic")
      # )
    )
  )
  
}, script = paste0(store, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store, "_target_script.R"), 
    store = store, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

cell_id_dict = tar_read(cell_id_dict, store = store)
cell_id_dict |> arrow::write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/dataset_cell_dict_v1_2_0_Jul_2024.parquet",
                                     compression = "zstd")


