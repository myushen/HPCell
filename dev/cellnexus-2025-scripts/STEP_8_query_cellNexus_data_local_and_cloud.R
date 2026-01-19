# This scripts test cellNexus API with new data generated in ~/git_control/cellNexus/dev/STEP_7_unify_cell_metadata.R
library(dplyr)
library(cellNexus)
library(stringr)

cache = "~/scratch/cache_temp"

x = get_metadata(cache_directory = cache, cloud_metadata = NULL, local_metadata = "~/scratch/cache_temp/metadata.1.0.0.parquet")
sce = x |>
  get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL)

# Check the number of cells per dataset
x |> dplyr::count(dataset_id)

# Check the number of cells per sample_id
x |> dplyr::count(sample_id) |> dplyr::count(n>10)

x |> filter(dataset_id == "601c7454-17ff-432c-8e0d-902bfddb8013") |> select(cell_id) |> arrange(cell_id)

