# ------------------------------------------------------------------------------
# Script Name: 02_compute_number_of_nodes_resumed.R
# Description: This script computes the number of internal nodes in the 
#              Maximum Clade Credibility (MCC) trees and consensus trees 
#              for each simulation and age, based on a directory of simulated 
#              phylogenetic trees. It performs the following steps:
#                - Locates paths to true, consensus, and MCC trees
#                - Loads true topologies to extract the deepest nodes
#                - Iterates over all consensus/MCC trees to count internal nodes
#                - Saves the resulting data as a .csv file for downstream analysis
# ------------------------------------------------------------------------------

library(here)
library(stringr)
library(tibble)
library(dplyr)
library(purrr)
library(ape)

source(".Rprofile")

# load the true trees
true_trees <- list.files(here("data/simulated-2025-05-13"), full.names = TRUE, recursive = F) |>
  keep(~ str_detect(.x, "beast-data-sim")) |>
  tibble(final = _) |>
  mutate(
    num = as.numeric(str_extract(final, "(\\d+)$")),
    path = str_glue("{final}/beast-data-sim-{num}-1/tree-sim-{num}-1.tree")
  ) |>
  pull(path)

# load true topologies 
true_topologies <- lapply(true_trees, function(path) read.tree(path))

# compute the deepest nodes
deepest_nodes <- lapply(true_topologies, function(tree) getNodesByDepth(tree)[2:11])

# paths of the mcc and the consensus trees
dir_path <- here("data/simulated-2025-05-13")
paths_consensus <- list.files(dir_path, full.names = T, recursive = T) |>
  keep(~ str_detect(.x, "consensus-"))

paths_mcc <- list.files(dir_path, full.names = T, recursive = T) |>
  keep(~ str_detect(.x, "mcc-"))

# compute the number of nodes as a function of the age of the tree -------------
df_number_of_nodes <- tibble(
  age = numeric(),
  simulation = numeric(),
  n_mcc = numeric(),
  n_consensus = numeric()
)

for (i in 1:length(paths_consensus)) {
  # read paths
  path_consensus <- paths_consensus[i]
  path_mcc <- paths_mcc[i]
  
  # information about the file
  tree_age <- as.integer(str_extract(path_mcc, "(?<=-)(\\d+)(?=\\.tree)"))
  tree_simulation_number <- as.numeric(str_match(path_mcc, "beast-data-sim-(\\d+)-\\d+")[, 2])
  
  # load trees
  tree_cs <- read.tree(path_consensus)
  tree_mcc <- read.tree(path_mcc)
  
  # implement the dataframe
  df_number_of_nodes <- bind_rows(
    df_number_of_nodes,
    tibble(
      age = tree_age,
      simulation = tree_simulation_number,
      n_mcc = tree_mcc$Nnode,
      n_consensus = tree_cs$Nnode
    )
  )
}

# save the dataframe
write.csv(
  df_number_of_nodes, 
  paste0(getwd(),"/output/results/number_nodes_mcc_cs.csv"),
  row.names = F,
)


