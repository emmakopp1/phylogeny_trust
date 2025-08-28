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
path_repository <-here("data/simulated-2025-07-28")

# compute the path for the csv output
output_path_mcc <- here("output/results/first_split_age_mcc.csv")
output_path_cs <- here("output/results/first_split_age_cs.csv")


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
  
  return(data.frame(
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
path_mcc <- path_phylo[grepl("mcc-", path_phylo)]
path_cs <- path_phylo[grepl("consensus-", path_phylo)]

# compute for each age, simulation the root and first split ages
# for the mcc trees
deepest_nodes_mcc <- lapply(path_mcc, function(path) first_split(path))
deepest_nodes_mcc <- do.call(rbind, deepest_nodes_mcc)

# for the consensus trees
deepest_nodes_cs <- lapply(path_cs, function(path) first_split(path))
deepest_nodes_cs <- do.call(rbind, deepest_nodes_cs)

# write the file 
write.csv(deepest_nodes_mcc, output_path_mcc)
write.csv(deepest_nodes_cs, output_path_cs)



