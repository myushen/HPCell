# This scripts test cellNexus API with new data generated in ~/git_control/cellNexus/dev/STEP_7_unify_cell_metadata.R
library(dplyr)
library(cellNexus)
library(stringr)

cache = "~/scratch/cache_temp"

x = get_metadata(cache_directory = cache, cloud_metadata = NULL, local_metadata = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/metadata.1.4.0.parquet")
x = x |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |> stringr::str_like("%10x%") & 
      tissue == "lung parenchyma" &
      cell_type |> stringr::str_like("%CD4%")
  ) |> filter(empty_droplet == FALSE, alive == TRUE, scDblFinder.class!="doublet")

# Test SCE
sce = x |> 
  get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL)

# TEST CPM
cpm = x |> 
  get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL, assays = "cpm")

# TEST RANK
rank = x |> 
  get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL, assays = "rank")

# TEST SCT
sct = x |> get_single_cell_experiment(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL, assays = "sct")

# TEST Pseudobulk
pseudobulk = x |> get_pseudobulk(cache_directory = "/vast/scratch/users/shen.m/cellNexus", repository = NULL)

# TEST metacell_2
metacell = x |> filter(!is.na(metacell_2)) |> get_metacell(cache_directory = "/vast/scratch/users/shen.m/cellNexus", 
                                                             repository = NULL, cell_aggregation = "metacell_2")


# Check the number of cells per dataset
x |> dplyr::count(dataset_id)

# Check the number of cells per sample_id
x |> dplyr::count(sample_id) |> dplyr::count(n>10)

# Check whether cell_id strategt is implemented
x |> select(cell_id) |> arrange(cell_id)

# Check QC metrics
x |> dplyr::count(empty_droplet, alive, scDblFinder.class)

# Check empty droplet ratio
x |> dplyr::count(empty_droplet) |>
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

# Check alive ratio
x |> filter(!empty_droplet) |> dplyr::count(alive) |> 
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

# Check doublet ratio
x |> filter(!empty_droplet, alive) |> dplyr::count(scDblFinder.class) |> 
  collect() |> mutate(n_cells = sum(n), pct = n / n_cells * 100)

