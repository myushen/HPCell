# # If only update metacell
# #use cloud 1.0.10 and left join cell_id with new metacell membership. update the stats. recalculate metacell_2
# cloud_meta = get_metadata(cloud_metadata = get_metadata_url(databases = "metadata.1.0.12.parquet"), cache_directory = "~/scratch/cache_temp")
# cloud_meta = cloud_meta |> select(-contains("metacell"))
# metacell = 
#   tar_read(metacell_tbl, store = "/vast/scratch/users/shen.m/cellNexus_target_store") |> 
#   bind_rows() |> 
#   dplyr::rename(cell_id = cell) |> 
#   dplyr::rename(metacell_2 = gamma2,
#                 metacell_4 = gamma4,
#                 metacell_8 = gamma8,
#                 metacell_16 = gamma16,
#                 metacell_32 = gamma32,
#                 metacell_64 = gamma64,
#                 metacell_128 = gamma128,
#                 metacell_256 = gamma256,
#                 metacell_512 = gamma512,
#                 metacell_1024 = gamma1024,
#                 metacell_2048 = gamma2048,
#                 metacell_4096 = gamma4096,
#                 metacell_8192 = gamma8192,
#                 metacell_16384 = gamma16384,
#                 metacell_32768 = gamma32768, 
#                 metacell_65536 = gamma65536)
# cloud_meta |> cellNexus:::duckdb_write_parquet("/vast/scratch/users/shen.m/cache_temp/temp.metadata.1.0.12.parquet")
# metacell |> arrow::write_parquet("/vast/scratch/users/shen.m/cache_temp/metacell_no_min_cell_thresh.parquet")
# job::job({
#   con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
#   
#   # Create a view for cell_annotation in DuckDB
#   dbExecute(con, "
#   CREATE VIEW cloud_metadata AS
#   SELECT
#     *
#   FROM read_parquet('/vast/scratch/users/shen.m/cache_temp/temp.metadata.1.0.12.parquet')
# ")
#   
#   dbExecute(con, "
#   CREATE VIEW metacell AS
#   SELECT
#     *
#   FROM read_parquet('/vast/scratch/users/shen.m/cache_temp/metacell_no_min_cell_thresh.parquet')
# ")
#   
#   query <- "
#     COPY (
#   SELECT cloud_metadata.*,
#          metacell.metacell_2,
#          metacell.metacell_4,
#           metacell.metacell_8,
#           metacell.metacell_16,
#           metacell.metacell_32,
#           metacell.metacell_64,
#           metacell.metacell_128,
#           metacell.metacell_256,
#           metacell.metacell_512,
#          metacell.metacell_1024,
#           metacell.metacell_2048,
#           metacell.metacell_4096,
#           metacell.metacell_8192,
#           metacell.metacell_16384,
#           metacell.metacell_32768,
#           metacell.metacell_65536
#   FROM cloud_metadata
#   LEFT JOIN metacell
#   ON cloud_metadata.cell_id = metacell.cell_id
#   ) TO '/vast/scratch/users/shen.m/cellNexus/no_min_cell_thresh_in_metacell_metadata.1.0.12.parquet'
#   (FORMAT PARQUET, COMPRESSION 'gzip');
# "
#   
#   # Execute the final query to write the result to a Parquet file
#   dbExecute(con, query)
#   
#   # Disconnect from the database
#   dbDisconnect(con, shutdown = TRUE)
#   
#   #system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_4.parquet box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
#   
#   print("Done.")
# })

# cell_metadata =
#   tbl(
#     dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#     sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/cellNexus/metacell_updated_metadata.1.0.10.parquet')")
#   )



# aggregate metacell from metadata and save assays
store_file_cellNexus = "/vast/scratch/users/shen.m/targets_prepare_database_split_datasets_chunked_1_3_0_metacell/" # MODIFY HERE: targets store directory for this pipeline
cell_metadata_path = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.1.4.0.parquet" # MODIFY HERE: cell metadata parquet (used to dynamically derive metacell column names)
# Dynamically get metacell column names from cell_metadata; fallback if file not available (e.g. local dev)
# metacell_columns <- "metacell_65536"
# tbl(
#     dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#     sql(glue::glue("SELECT * FROM read_parquet('{cell_metadata_path}')"))
#   ) |> select(contains("metacell_")) |> colnames()

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  library(tidySingleCellExperiment)
  library(SingleCellExperiment)
  library(tidyverse)
  library(glue)
  library(digest)
  library(scater)
  library(arrow)
  library(dplyr)
  library(duckdb)
  library(cellNexus)
  library(BiocParallel)
  library(parallelly)
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    cue = tar_cue(mode = "thorough"),
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
            memory_gigabytes_required = c(20, 40, 80, 120), 
            #memory_gigabytes_required = c(100, 200, 300, 400), 
            cpus_per_task = c(2, 2, 5, 10), 
            time_minutes = c(60, 60, 60*4, 60*4),
            verbose = T
          )
        ),
        crew_controller_slurm(
          name = "tier_4",
          workers = 200,
          tasks_max = 10,
          crashes_error = 5, 
          seconds_idle = 30,
          options_cluster = crew_options_slurm(
            #memory_gigabytes_required = c(120, 150, 200, 250, 300), 
            memory_gigabytes_required = c(40, 70, 90, 150, 300), 
            cpus_per_task = c(2, 2, 5, 5, 5, 10), 
            time_minutes = c(60*24,60*24,60*24,60*24,60*24),
            verbose = T
          )
        ))),
    trust_object_timestamps = TRUE)
  
  get_ids <- function(cell_metadata, metacell_column) {
    metacell_column <- as.character(metacell_column)
    tbl(dbConnect(duckdb::duckdb(), dbdir = ":memory:"), 
        sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      filter(!is.na(.data[[metacell_column]])) |>
      distinct(file_id_cellNexus_single_cell) |> pull()
  }
  
  get_sce <- function(cell_metadata, id, metacell_column, cache) {
    metacell_column <- as.character(metacell_column)
    tbl(dbConnect(duckdb::duckdb(), 
                  dbdir = ":memory:"), sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))) |> 
      filter(!is.na(.data[[metacell_column]])) |>
      filter(file_id_cellNexus_single_cell == id) |> get_single_cell_experiment(cache_directory = cache, 
                                                                                repository = NULL # this assume SCE are not uploaded to cloud
                                                                                )
  }
  
  aggregate_metacell <- function(sce, metacell) {
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)) - 1
    bp <- MulticoreParam(workers = cores, progressbar = TRUE)
    aggregate_metacell <- aggregateAcrossCells(sce, colData(sce)[, c("sample_id", metacell)], BPPARAM = bp)
    aggregate_metacell = aggregate_metacell |> mutate(cell_id = paste(sample_id, .data[[metacell]], sep = "___"))
    # Assign cell_id to SCE metadata rownames
    rownames(colData(aggregate_metacell)) <- aggregate_metacell$cell_id
    aggregate_metacell = aggregate_metacell |> select(-contains(".1"))
    aggregate_metacell
    
  }
  
  save_anndata = function(sce, cache_directory) {
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    file_id = pull(distinct(sce, file_id_cellNexus_single_cell))
    cellNexus:::write_h5ad(sce, glue("{cache_directory}/{file_id}"))
    return(TRUE)
  }
  
  c(
    list(
      tar_target(
        cell_metadata,
        "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.1.4.0.parquet", # MODIFY HERE: cell metadata parquet (must match cell_metadata_path above)
        deployment = "main",
        packages = c("arrow", "dplyr", "duckdb")
      ),
      tar_target(
        local_cache,
        "/vast/scratch/users/shen.m/cellNexus", # MODIFY HERE: local cache directory containing the single-cell h5ad files (input to get_single_cell_experiment)
        deployment = "main"
      ),
      tar_target(
        save_cache_directory,
        "/vast/scratch/users/shen.m/cellNexus/cellxgene/01-07-2024", # MODIFY HERE: output directory where aggregated metacell anndata files are saved
        deployment = "main"
      )
    ),
    tarchetypes::tar_map(
      values = tibble(
        metacell_column = tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.1.4.0.parquet')") # MODIFY HERE: cell metadata parquet path inside SQL (must match cell_metadata_path above)
        ) |> select(contains("metacell_")) |> colnames()
      ),
      names = metacell_column,
      unlist = TRUE,
      tar_target(
        file_ids,
        get_ids(cell_metadata, metacell_column)
        # |>
        #   # TEST PURPOSE ONLY
        #   head(2)
      ),
      tar_target(
        file_id_sce,
        get_sce(cell_metadata, file_ids, metacell_column, local_cache),
        pattern = map(file_ids),
        resources = tar_resources(crew = tar_resources_crew(controller = "elastic"))
      ),
      tar_target(
        metacell,
        aggregate_metacell(file_id_sce, metacell_column),
        pattern = map(file_id_sce),
        resources = tar_resources(crew = tar_resources_crew(controller = "tier_4"))
      ),
      tar_target(
        save_metacell,
        save_anndata(metacell, paste0(save_cache_directory, "/", metacell_column, "/counts")),
        pattern = map(metacell),
        resources = tar_resources(crew = tar_resources_crew(controller = "tier_4"))
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
#tar_invalidate(names = everything(), store = store_file_cellNexus)
tar_meta(store = store_file_cellNexus) |> filter(!is.na(error)) |> distinct(name, error)
tar_errored(store = store_file_cellNexus)
# With tar_map, target names are suffixed by metacell_column, e.g. file_ids_metacell_4, metacell_metacell_4
# tar_workspace("metacell_metacell_4_<hash>", store = store_file_cellNexus, script = paste0(store_file_cellNexus, "_target_script.R"))
# debugonce(aggregate_metacell)
# aggregate_metacell(file_id_sce, "metacell_4")  # when debugging a specific branch

# Check the number of file id should be created for metacell_2
cache_dir = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/" # MODIFY HERE: directory used for verification queries below (contains metadata.1.4.0.parquet)
file_count <- function(metacell_vec){
  metacell_column <- ensym(metacell_vec)
  get_metadata(
    cache_directory = cache_dir,
    cloud_metadata = NULL,
    local_metadata = file.path(cache_dir, "metadata.1.4.0.parquet")
  ) |> 
    filter(!is.na(metacell_column),
           empty_droplet == FALSE,
           alive == TRUE,
           scDblFinder.class != "doublet") |>
    mutate(combo = paste(sample_id, metacell_column, sep = "___")) |> 
    dplyr::count(combo, file_id_cellNexus_single_cell) |> 
    distinct(file_id_cellNexus_single_cell) |> dplyr::count() |> as_tibble()
}
file_count('metacell_2')

# Unit test query lung tissue, metacell 256
lung_metacell_256 = get_metadata(
  cache_directory = cache_dir,
  cloud_metadata = NULL,
  local_metadata = file.path(cache_dir, "metadata.1.4.0.parquet")
) |> 
  filter(!is.na(metacell_256),
         empty_droplet == FALSE,
         alive == TRUE,
         scDblFinder.class != "doublet") |>
  filter(tissue  == "lung") |> 
  get_metacell(cache_directory = cache_dir, cell_aggregation = "metacell_256")

# Cell number in lung_metacell_256 should match the count below:
get_metadata(
  cache_directory = cache_dir,
  cloud_metadata = NULL,
  local_metadata = file.path(cache_dir, "metadata.1.4.0.parquet")
) |> 
  filter(!is.na(metacell_256),
         empty_droplet == FALSE,
         alive == TRUE,
         scDblFinder.class != "doublet") |>
  filter(tissue  == "lung") |> 
  distinct(sample_id, metacell_256, cell_type_unified_ensemble) |> dplyr::count()


