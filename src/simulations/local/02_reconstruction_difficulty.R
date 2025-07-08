# ------------------------------------------------------------------------------
# Script Name: first_split_analysis.R
# Description: Analysis of first divergence ages (root split) in simulated 
#              phylogenetic trees. This script extracts for each reference tree
#              (tree-sim) the root age and the age of the first speciation event
#              (first divergence). Results include simulation ID, tree age, root
#              node, root age, first split node, and corresponding split age.
#              Data are exported to a CSV file for downstream analyses.
# Input:       Simulated phylogenetic trees (.tree format) in the directory
#              specified by path_repository
# Output:      CSV file containing divergence metrics for each analyzed 
#              simulation 
# ------------------------------------------------------------------------------

library(here)
library(ape)
library(stringr)
library(castor)
library(parallel)
library(phangorn)

# functions --------------------------------------------------------------------
# select the path repository of your analysis 
#path_repository <-here("data/simulated-2025-05-13")
#path_repository <- here("data/simulated-2025-07-02-6000")
path_repository <- here("data/simulated-2025-07-04-12000")


# number of different ages per simulation
#N_ages <- 17
N_ages <- 1

# compute the path for the csv output
#output_path <- here("output/results/first_split_age.csv")
#output_path <- here("output/results/first_split_age_6000.csv")
output_path <- here("output/results/first_split_age_12000.csv")


# function to get the outgroup of the tree 
first_split <- function(path){
  tree = read.tree(path)
  
  tree_age <- as.integer(str_extract(path, "(?<=-)(\\d+)(?=\\.tree)"))
  tree_simulation_number <- as.numeric(str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2])
  
  # root of the tree
  root = find_root(tree)
  # children of the root
  children = tree$edge[which(tree$edge[,1]==root),2]
  
  t = max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[children]
  
  return(list(
    age = tree_age,
    simulation = tree_simulation_number,
    root = root,
    root_age = max(node.depth.edgelength(tree)),
    first_split = children[which.max(t)],  
    root_split_age = t[which.max(t)]
  ))  
}


# load variables ---------------------------------------------------------------
# general path
path_phylo <- list.files(
  path_repository, 
  full.names = TRUE, 
  recursive = TRUE
)

# paths to true topologies, posterior and mcc  
path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)]
path_trees_true <- path_trees_true[seq(1, length(path_trees_true), by = N_ages)]

# compute for each age, simulation the root and first split ages
deepest_nodes_mcc <- lapply(path_trees_true, function(path) first_split(path))
deepest_nodes_mcc <- do.call(rbind, deepest_nodes_mcc)

# write the file 
write.csv(deepest_nodes_mcc, output_path)
