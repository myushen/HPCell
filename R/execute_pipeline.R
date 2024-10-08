#' Run Targets Pipeline for HPCell
#'
#' @description
#' This function sets up and executes a `targets` pipeline for HPCell. It saves input data and configurations, 
#' writes a pipeline script, and runs the pipeline using the 'targets' package.
#'
#' @param input_data Input data for the pipeline.
#' @param store Directory path for storing the pipeline files.
#' @param input_reference Optional reference data.
#' @param tissue Tissue type for the analysis.
#' @param computing_resources Configuration for computing resources.
#' @param debug_step Optional step for debugging.
#' @param filter_empty_droplets Flag to indicate if input filtering is needed.
#' @param RNA_assay_name Name of the RNA assay.
#' @param sample_column Column name for sample identification.
#' @param cell_type_annotation_column Column name for cell type annotation in input data
#' @param data_container_type A character vector of length one specifies the input data type.
#' @param profiler Optional step for profilling. Default is FALSE
#' data type can be one of the following: anndata for annotated data mainly used in python.
#' sce_rds and seurat_rds for `SingleCellExperiment` and `Seurat` RDS format representively
#' seurat_rds for `Seurat` RDS format.
#' sce_hdf5 for `SingleCellExperiment` HDF5 format 
#' seurat_hdf5 for `Seurat` HDF5 format 
#'
#' @return The output of the `targets` pipeline, typically a pre-processed data set.
#'
#' @importFrom glue glue
#' @importFrom targets tar_script
#' @import crew.cluster
#' @import tarchetypes
#' @import targets
#' @import broom
#' @import ggplot2
#' @import ggupset
#' @import here
#' @import qs
#' @import crew
#' @importFrom future tweak
#' @import crew 
#' @import crew.cluster
#' @export
run_targets_pipeline <- function(
    input_data, 
    store =  "./", 
    input_reference = NULL,
    tissue,
    computing_resources = crew_controller_local(workers = 1), 
    debug_step = NULL,
    filter_empty_droplets = NULL, 
    RNA_assay_name = "RNA", 
    sample_column = "sample", 
    cell_type_annotation_column = "Cell_type_in_each_tissue",
    data_container_type
){
  
  # Fix GCHECKS 
  data_object <- NULL
  reference_file <- NULL
  tissue_file <- NULL
  filtered_file <- NULL
  sample_column_file <- NULL
  cell_type_annotation_column_file <- NULL
  reference_label_coarse <- NULL
  reference_label_fine <- NULL
  input_read <- NULL
  unique_tissues <- NULL
  reference_read <- NULL
  empty_droplets_tbl <- NULL
  cell_cycle_score_tbl <- NULL
  annotation_label_transfer_tbl <- NULL
  alive_identification_tbl <- NULL
  doublet_identification_tbl <- NULL
  non_batch_variation_removal_S <- NULL
  preprocessing_output_S <- NULL
  create_pseudobulk_sample <- NULL
  sampleName <- NULL
  cellAnno <- NULL
  pseudobulk_merge_all_samples <- NULL
  calc_UMAP_dbl_report <- NULL
  variable_gene_list <- NULL
  tar_render <- NULL
  empty_droplets_report <- NULL
  doublet_identification_report <- NULL
  Technical_variation_report <- NULL
  pseudobulk_processing_report <- NULL
  
  sample_column = enquo(sample_column)
  # cell_type_annotation_column = enquo(cell_type_annotation_column)
  
  # Save inputs for passing to targets pipeline 
  # input_data |> CHANGE_ASSAY |> saveRDS("input_file.rds")
  input_data |> saveRDS("input_file.rds")
  input_reference |> saveRDS("input_reference.rds")
  tissue |> saveRDS("tissue.rds")
  computing_resources |> saveRDS("temp_computing_resources.rds")
  filter_empty_droplets |> saveRDS("filter_empty_droplets.rds")
  sample_column |> saveRDS("sample_column.rds")
  cell_type_annotation_column |> saveRDS("cell_type_annotation_column.rds")
  data_container_type |> saveRDS("data_container_type.rds")
  debug_step |> saveRDS("debug_step_param.rds")
  # Write pipeline to a file
  tar_script({
    # library(targets)
    # library(tarchetypes)
    # library(crew)
    # library(crew.cluster)
    
    computing_resources = readRDS("temp_computing_resources.rds")
    debug_step = readRDS("debug_step_param.rds")
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set(
      packages = c(
        "HPCell",
        "readr",
        "dplyr",
        "tidyr",
        "ggplot2",
        "purrr",
        "Seurat",
        "tidyseurat",
        "glue",
        "scater",
        "DropletUtils",
        "EnsDb.Hsapiens.v86",
        "here",
        "stringr",
        "readr",
        "rlang",
        "scuttle",
        "scDblFinder",
        "ggupset",
        "tidySummarizedExperiment",
        "broom",
        "tarchetypes",
        "SeuratObject",
        "SingleCellExperiment", 
        "SingleR", 
        "celldex", 
        "tidySingleCellExperiment", 
        "tibble", 
        "magrittr",
        "qs", 
        "S4Vectors",
        "tarprof",
        "zellkonverter"
      ),
      memory = "transient",
      garbage_collection = TRUE,
      #trust_object_timestamps = TRUE,
      storage = "worker", 
      retrieval = "worker", 
      #error = "continue",         
      format = "qs", 
      debug = debug_step, # Set the target you want to debug.
      # cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
      controller = computing_resources
    )
    
    #-----------------------#
    # Future SLURM
    #-----------------------#
    
    # library(future)
    # library("future.batchtools")
    # slurm <-
    #     `batchtools_slurm` |>
    #     future::tweak( template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
    #                                  resources=list(
    #                                      ncpus = 20,
    #                                      memory = 6000,
    #                                      walltime = 172800
    #                                  )
    #     )
    # plan(slurm)
    
    # small_slurm =
    #   tar_resources(
    #     future = tar_resources_future(
    #       plan = tweak(
    #         batchtools_slurm,
    #         template = "dev/slurm_batchtools.tmpl",
    #         resources = list(
    #           ncpus = 2,
    #           memory = 40000,
    #           walltime = 172800
    #         )
    #       )
    #     )
    #   )
    #
    # big_slurm =
    #   tar_resources(
    #     future = tar_resources_future(
    #       plan = tweak(
    #         batchtools_slurm,
    #         template = "dev/slurm_batchtools.tmpl",
    #         resources = list(
    #           ncpus = 19,
    #           memory = 6000,
    #           walltime = 172800
    #         )
    #       )
    #     )
    #   )
    
    target_list = list(
      tar_target(file, "input_file.rds", format = "rds"), 
      tar_target(data_object, readRDS("input_file.rds")),
      #tar_target(reference_file, "input_reference.rds", format = "rds"), 
      tar_target(reference_file, readRDS("input_reference.rds")), 
      tar_target(tissue_file, readRDS("tissue.rds")), 
      tar_target(filtered_file, readRDS("filter_empty_droplets.rds")), 
      tar_target(sample_column_file, readRDS("sample_column.rds")), 
      tar_target(cell_type_annotation_column_file, readRDS("cell_type_annotation_column.rds")),
      tar_target(data_container_type_file, readRDS("data_container_type.rds")))
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    target_list|> c(list(
      
      # Define input files
      # tarchetypes::tar_files(name= input_track, 
      #                        data_object, 
      #                        deployment = "main"),
      # tarchetypes::tar_files(name= reference_track,
      #                        read_reference_file, 
      #                        deployment = "main"),
      tar_target(filter_empty_droplets, filtered_file, deployment = "main"),
      tar_target(tissue, tissue_file, deployment = "main", ),
      tar_target(sample_column, sample_column_file, deployment = "main"),
      tar_target(cell_type_annotation_column, cell_type_annotation_column_file, deployment = "main"),
      tar_target(reference_label_coarse, reference_label_coarse_id(tissue), deployment = "main"), 
      tar_target(reference_label_fine, reference_label_fine_id(tissue), deployment = "main"), 
      # Reading input files
      tar_target(file_path, data_object, pattern = map(data_object), format = "file", deployment = "main"),
      tar_target(unique_tissues,
                 get_unique_tissues(read_data_container(file_path, container_type = data_container_type_file), sample_column |> quo_name()),
                 pattern = map(file_path),
                 iteration = "list"),
      # tar_target(
      #   tissue_subsets,
      #   input_read, split.by = "Tissue"),
      #   pattern = map(input_read),
      #   iteration = "list"
      # ),
      tar_target(reference_read, reference_file, deployment = "main"),
      
      # Identifying empty droplets
      tar_target(empty_droplets_tbl,
                 empty_droplet_id(read_data_container(file_path, container_type = data_container_type_file), filter_empty_droplets),
                 pattern = map(file_path),
                 iteration = "list"),
      
      # Cell cycle scoring
      tar_target(cell_cycle_score_tbl, cell_cycle_scoring(read_data_container(file_path, container_type = data_container_type_file ),
                                                          empty_droplets_tbl),
                 pattern = map(file_path,
                               empty_droplets_tbl),
                 iteration = "list"),
      
      # Annotation label transfer
      tar_target(annotation_label_transfer_tbl,
                 annotation_label_transfer(read_data_container(file_path, container_type = data_container_type_file),
                                           empty_droplets_tbl,
                                           reference_read),
                 pattern = map(file_path,
                               empty_droplets_tbl),
                 iteration = "list"),
      
      # Alive identification
      tar_target(alive_identification_tbl, alive_identification(read_data_container(file_path, container_type = data_container_type_file),
                                                                empty_droplets_tbl,
                                                                annotation_label_transfer_tbl),
                 pattern = map(file_path,
                               empty_droplets_tbl,
                               annotation_label_transfer_tbl),
                 iteration = "list"),
      
      # Doublet identification
      tar_target(doublet_identification_tbl, doublet_identification(read_data_container(file_path, container_type = data_container_type_file),
                                                                    empty_droplets_tbl,
                                                                    alive_identification_tbl,
                                                                    annotation_label_transfer_tbl,
                                                                    reference_label_fine),
                 pattern = map(file_path,
                               empty_droplets_tbl,
                               alive_identification_tbl,
                               annotation_label_transfer_tbl),
                 iteration = "list"),
      
      # Non-batch variation removal
      tar_target(non_batch_variation_removal_S, non_batch_variation_removal(read_data_container(file_path, container_type = data_container_type_file),
                                                                            empty_droplets_tbl,
                                                                            alive_identification_tbl,
                                                                            cell_cycle_score_tbl),
                 pattern = map(file_path,
                               empty_droplets_tbl,
                               alive_identification_tbl,
                               cell_cycle_score_tbl),
                 iteration = "list"),
      
      # Pre-processing output
      tar_target(preprocessing_output_S, preprocessing_output(tissue,
                                                              non_batch_variation_removal_S,
                                                              alive_identification_tbl,
                                                              cell_cycle_score_tbl,
                                                              annotation_label_transfer_tbl,
                                                              doublet_identification_tbl),
                 pattern = map(non_batch_variation_removal_S,
                               alive_identification_tbl,
                               cell_cycle_score_tbl,
                               annotation_label_transfer_tbl,
                               doublet_identification_tbl),
                 iteration = "list")
      
      # pseudobulk preprocessing for each sample 
      # tar_target(create_pseudobulk_sample, create_pseudobulk(preprocessing_output_S,
      #                                                        assays = "SCT",
      #                                                        cell_type_annotation_column,
      #                                                        x = c(sampleName, cellAnno)),
      #            pattern = map(preprocessing_output_S),
      #            iteration = "list"),
      # 
      # tar_target(pseudobulk_merge_all_samples, pseudobulk_merge(create_pseudobulk_sample,
      #                                                           assays = "RNA",
      #                                                           x = c(sampleName)),
      #            iteration = "list"),
      # 
      # tar_target(calc_UMAP_dbl_report, calc_UMAP(input_read),
      #            pattern = map(input_read),
      #            iteration = "list"),
      # tar_target(variable_gene_list, find_variable_genes(input_read,
      #                                                    empty_droplets_tbl),
      #            pattern = map(input_read, empty_droplets_tbl),
      #            iteration = "list")
      
      # tar_render(
      #   name = empty_droplets_report, # The name of the target
      #   path =  paste0(system.file(package = "HPCell"), "/rmd/Empty_droplet_report.Rmd"),
      #   params = list(x1 = tar_read(input_read, store = store),
      #                 x2 = tar_read(empty_droplets_tbl, store = store),
      #                 x3 = tar_read(annotation_label_transfer_tbl, store = store),
      #                 x4 = tar_read(unique_tissues, store = store),
      #                 x5 = sample_column |> quo_name())
      # ),
      # tar_render(
      #   name = doublet_identification_report,
      #   path = paste0(system.file(package = "HPCell"), "/rmd/Doublet_identification_report.Rmd"),
      #   params = list(x1 = input_read,
      #                 x2 = calc_UMAP_dbl_report,
      #                 x3 = doublet_identification_tbl,
      #                 x4 = annotation_label_transfer_tbl,
      #                 x5 = sample_column |> quo_name(),
      #                 x6 = cell_type_annotation_column |> quo_name())
      # ),
      # tar_render(
      #   name = Technical_variation_report,
      #   path =  paste0(system.file(package = "HPCell"), "/rmd/Technical_variation_report.Rmd"),
      #   params = list(x1= input_read,
      #                 x2= empty_droplets_tbl,
      #                 x3 = variable_gene_list,
      #                 x4 = calc_UMAP_dbl_report, 
      #                 x5 = sample_column |> quo_name())
      # ),
      # tar_render(
      #   name = pseudobulk_processing_report,
      #   path = paste0(system.file(package = "HPCell"), "/rmd/pseudobulk_analysis_report.Rmd"),
      #   params = list(x1 = pseudobulk_merge_all_samples, 
      #                 x2 = sample_column |> quo_name(), 
      #                 x3 = cell_type_annotation_column |> quo_name())
    )
    )
  }, script = glue("{store}.R"), ask = FALSE)
  
  #Running targets 
  # input_files<- c("CB150T04X__batch14.rds","CB291T01X__batch8.rds")
  # run_targets <- function(input_files){
  #   tar_make(
  #     script = glue("{store}.R"),
  #     store = store
  #   )
  # }
  # run_targets(input_files)
  
  tar_make(
    script = glue("{store}.R"),
    store = store,
    callr_function = NULL
  )
  # tar_make_future(
  #   script = glue("{store}.R"),
  #   store = store, 
  #   workers = 200, 
  #   garbage_collection = TRUE
  # )
  
  message(glue("HPCell says: you can read your output executing tar_read(preprocessing_output_S, store = \"{store}\") "))
  #tar_meta_download(store = store)
  metadata<- tar_meta(store = store)
  return(metadata)
  #tar_read(preprocessing_output_S, store = store)
  
}

## my_results = run_targets_pipeline(..)