# ------------------------------------------------------------------------------
# Script Name: shared_cognates.R
# Run : local
# Description: Analyzes shared cognates between subgroups in phylogenetic datasets
# -----------------------------------------------------------------------------------------

library(ape)
library(xml2)
library(here)
library(phangorn)
library(tidyverse)

# load data
# select the path repository of your analysis 
path_repository <- here("data/simulated-2025-07-28")

# compute the path for the csv output
output_path <- here("output/results/shared_cognate_summary_table.csv")

age_init_sim <- 1:17


# data -------------------------------------------------------------------------
# load the true trees
true_trees <- list.files(path_repository, full.names = TRUE, recursive = F) |>
  keep(~ str_detect(.x, "beast-data-sim")) |>
  tibble(final = _) |>
  expand_grid(age_init_sim = age_init_sim) |>
  mutate(
    num =  as.numeric(str_extract(final, "(\\d+)$")),
    path = str_glue("{final}/beast-data-sim-{num}-{age_init_sim}/tree-sim-{num}-{age_init_sim}.tree")
  ) |>
  pull(path)

# load data files 
paths_data <- list.files(path_repository, full.names = TRUE, recursive = F) |>
  keep(~ str_detect(.x, "beast-data-sim")) |>
  tibble(final = _) |>
  expand_grid(age_init_sim = age_init_sim) |>  # Créer toutes les combinaisons
  mutate(
    num = as.numeric(str_extract(final, "(\\d+)$")),
    path = str_glue("{final}/beast-data-sim-{num}-{age_init_sim}/beast-simulated-seq-{num}-{age_init_sim}.xml")
  ) |>
  pull(path)

# function to analyze the dataset
analyze_dataset <- function(tree_path, data_path) {
  
  tree_age <- as.integer(str_extract(tree_path, "(?<=-)(\\d+)(?=\\.tree)"))
  tree_simulation_number <- as.numeric(str_match(tree_path, "beast-data-sim-(\\d+)-\\d+")[, 2])
  
  # Load data
  tree <- read.tree(tree_path)
  data <- read_xml(data_path)
  
  # Get subgroups
  first_split <- Descendants(tree, tree$Nnode + 2, type = "children")
  subgroups <- lapply(first_split, function(node) {
    descendants <- tree$tip.label[unlist(allDescendants(tree)[node])]
    descendants[!is.na(descendants)]
  })
  
  # Create dataframe with cognate data
  df_subgroups <- map_df(seq_along(subgroups), function(i) {
    tibble(
      taxon = subgroups[[i]],
      subgroup = i
    )}) |>
    mutate(value = map_chr(taxon, function(taxon_name) {
      xpath_query <- str_glue("//sequence[@taxon='{taxon_name}']")
      node <- xml_find_first(data, xpath_query)
      
      if (length(node) > 0) {
        xml_attr(node, "value")
      } else {
        NA_character_
      }
    })) |> 
    mutate(positions_1 = map(value, ~which(str_split(.x, "")[[1]] == "1")))
  
  # Create similarity matrix
  similarity_df <- expand_grid(
    taxon1 = df_subgroups$taxon,
    taxon2 = df_subgroups$taxon
  ) |> 
    left_join(df_subgroups |> select(taxon, positions_1), 
              by = c("taxon1" = "taxon")) |>
    rename(positions_1_taxon1 = positions_1) |>
    left_join(df_subgroups |> select(taxon, positions_1), 
              by = c("taxon2" = "taxon")) |>
    rename(positions_1_taxon2 = positions_1) |>
    mutate(n_shared = map2_int(positions_1_taxon1, positions_1_taxon2, 
                               ~length(intersect(.x, .y)))) |>
    select(taxon1, taxon2, n_shared) |>
    pivot_wider(names_from = taxon2, values_from = n_shared) |> 
    column_to_rownames("taxon1")
  
  # Calculate shared cognates between subgroups
  n_shared <- df_subgroups |> 
    group_by(subgroup) |>
    summarise(all = list(reduce(positions_1, union))) |> 
    pull(all) 
  
  n_shared <- length(intersect(n_shared[[1]], n_shared[[2]]))/ length(union(n_shared[[1]], n_shared[[2]]))
  
  # Return results
  list(
    age = tree_age, 
    simu = tree_simulation_number,
    df_subgroups = df_subgroups,
    similarity_matrix = similarity_df,
    n_shared_cognates = n_shared,
    tree = tree
  )
}

# define datasets
datasets <- tibble(
  tree_path = true_trees,
  data_path = paths_data
)

# analyze all datasets
results <- datasets |> 
  mutate(analysis = pmap(list(tree_path, data_path), analyze_dataset))


# summary table 
summary_table <- results |> 
  mutate(
    simu = map_int(analysis, ~.x$simu),
    age = map_int(analysis, ~.x$age),
    n_shared = map_dbl(analysis, ~.x$n_shared_cognates),
    n_group1 = map_int(analysis, ~sum(.x$df_subgroups$subgroup == 1)),
    n_group2 = map_int(analysis, ~sum(.x$df_subgroups$subgroup == 2))
  ) |>
  select(simu, age, n_shared, n_group1, n_group2) |>
  group_by(age) |> 
  summarise(n_shared_mean = mean(n_shared, na.rm=T))


# save results
saveRDS(results, here("output/results/shared_cognates.rds"))

# save results
write.csv(summary_table, 
          output_path, 
          row.names = FALSE)



