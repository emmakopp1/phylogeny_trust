# ------------------------------------------------------------------------------
# Script Name: 02_marginal_prob_first_split_resumed.R
# Description: This script identifies the first split (deepest node) in MCC and 
#              consensus trees inferred from posterior phylogenies. For each tree:
#                - It locates the outgroup (first split)
#                - Retrieves the posterior support (clade frequency) of that split
#                - Extracts simulation metadata (age, simulation number)
#                - Outputs a CSV summarizing these values for all simulations
#              Processing is parallelized across available CPU cores.
# ------------------------------------------------------------------------------

library(here)
library(ape)
library(castor)
library(stringr)
library(parallel)
library(tibble)
library(treeio)


# choose the paths adapted to the simulation -----------------------------------
path_repository <- here("data/simulated-2025-07-28")
#path_repository <- here("data/simulated-2025-07-22-1500")
#path_repository <- here("data/simulated-2025-07-22-6000")
#path_repository <- here("data/simulated-2025-07-22-12000")

# compute the path for the csv output
output_path <- here("output/results/marginal_prob_first_split_mcc_consensus_hipstr.csv")
#output_path <- here("output/results/marginal_prob_first_split_mcc_consensus_1500.csv")
#output_path <- here("output/results/marginal_prob_first_split_mcc_consensus_6000.csv")
#output_path <- here("output/results/marginal_prob_first_split_mcc_consensus_12000.csv")

# functions --------------------------------------------------------------------
# function to get the outgroup of the tree 
first_split <- function(tree){
  # root of the tree
  root = find_root(tree)
  # children of the root
  children = tree$edge[which(tree$edge[,1]==root),2]
  # fisr split : children with the maximum distance to the tips
  t = max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[children]
  
  return(children[which.max(t)])
}


# load variables ---------------------------------------------------------------

#cluster_directory <- getwd()
path_phylo <- list.files(
  path_repository, 
  full.names = TRUE, 
  recursive = TRUE
)

# paths to true topologies, posterior and mcc  
path_trees_mcc <- path_phylo[grepl("mcc-", path_phylo)]
path_trees_cs <- path_phylo[grepl("consensus-", path_phylo)]
path_trees_hipstr <- path_phylo[grepl("hipstr-", path_phylo)]

# load the data
trees_cs <- lapply(path_trees_cs, read.tree)
trees_mcc <- lapply(path_trees_mcc, read.tree)
trees_hipstr <- lapply(path_trees_hipstr, read.beast)

deepest_nodes_mcc <- sapply(trees_mcc, function(tree) first_split(tree))
deepest_nodes_cs <- sapply(trees_cs, function(tree) first_split(tree))
deepest_nodes_hipstr <- sapply(trees_hipstr, function(tree) first_split(tree@phylo))

# data frame of results 
df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(df) <- c('age', 
                  'simulation', 
                  'node_cs',
                  'node_mcc',
                  'node_hipstr',
                  'cs_prob', 
                  'mcc_prob',
                  'hipstr_prob')

# initialization of the file
write.csv(
  df, 
  output_path,
  row.names = FALSE,
)

# main function taking as input an increment from 1 to length(path_trees_phylo)
process_file <- function(i){
  
  tryCatch({
    
    # load the mcc, the posterior and the deepest node of the outgroup of the true treee
    cs <- trees_cs[[i]]
    mcc <- trees_mcc[[i]]
    hipstr <- trees_hipstr[[i]]@phylo
    hipstr_node_label = trees_hipstr[[i]]@data |> 
      as_tibble() |> 
      select(node,posterior) |> 
      filter(!is.na(posterior))
    
    node_cs <- deepest_nodes_cs[i]
    node_mcc <- deepest_nodes_mcc[i]
    node_hipstr <- deepest_nodes_hipstr[i]
    
    # load the paths, its simulation number and its age
    path <- path_trees_mcc[[i]]
    tree_simulation_number <- as.numeric(str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2])
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))
    
    # posterior probability of the outgroup for the consensus tree
    internal_nodes_cs <- (length(cs$tip.label) + 1):(length(cs$tip.label) + cs$Nnode)
    internal_nodes_mcc <- (length(mcc$tip.label) + 1):(length(mcc$tip.label) + mcc$Nnode)
    internal_nodes_hipstr <- (length(hipstr$tip.label) + 1):(length(hipstr$tip.label) + hipstr$Nnode)
    
    prob_cs <- cs$node.label[match(node_cs,internal_nodes_cs)]
    prob_mcc <- mcc$node.label[match(node_mcc,internal_nodes_mcc)]
    prob_hipstr <- hipstr_node_label |> 
      filter(node == node_hipstr) |> 
      select(posterior) |> 
      mutate(posterior = round(as.numeric(posterior), 3)) |> 
      as.character()
    
    # report age, simulation, the first split, and the clade frequencies of the first split 
    row <- data.frame(
      age = tree_age,
      simulation = tree_simulation_number,
      node_cs = node_cs,
      node_mcc = node_mcc,
      node_hipstr = node_hipstr,
      cs_prob = round(as.numeric(prob_cs), 3),
      mcc_prob = round(as.numeric(prob_mcc), 3),
      hipstr_prob = round(as.numeric(prob_hipstr), 3)
    )
    
    # write the results
    write.table(
      row,
      file = output_path,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(output_path),
      append = TRUE,
      quote = FALSE
    )
    
    
  }, error = function(e) {})
}

# local
res_list <- lapply(seq_along(path_trees_cs), process_file)


