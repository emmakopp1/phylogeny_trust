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
library(stringr)
library(parallel)

# functions --------------------------------------------------------------------
# function to get the outgroup of the tree 
first_split <- function(tree){
  # root of the tree
  root = find_root(tree)
  # children of the root
  children = tree$edge[which(tree$edge[,1]==root),2]
  # fisr split : children with the maximum distance to the tips
  #first_split = children[which.max(tree$edge.length[children])]
  
  t = max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[children]
  #children[which.max(t)]
  
  return(children[which.max(t)])
}

# arguements -------------------------------------------------------------------
start <- 1
end <- 850
cat("Traitement des fichiers de", start, "à", end, "\n")


# load variables ---------------------------------------------------------------
# general path
cluster_directory <- getwd()
path_phylo <- list.files(
  paste0(cluster_directory,"/data/simulated-2025-05-13"), 
  full.names = TRUE, 
  recursive = TRUE
)

# paths to true topologies, posterior and mcc  
path_trees_mcc <- path_phylo[grepl("mcc", path_phylo)][start:end] 
path_trees_cs <- path_phylo[grepl("consensus", path_phylo)][start:end] 

# load the data
trees_cs <- lapply(path_trees_cs, read.tree)
trees_mcc <- lapply(path_trees_mcc, read.tree)
deepest_nodes_mcc <- sapply(trees_mcc, function(tree) first_split(tree))
deepest_nodes_cs <- sapply(trees_cs, function(tree) first_split(tree))

# data frame of results 
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c('age', 'simulation', 'node_cs','node_mcc', 'cs_prob', 'mcc_prob')

# initialization of the file
write.csv(
  df, 
  paste0(cluster_directory,
         sprintf("/output/results/marginal_prob_first_split_mcc_consensus_%d_%d.csv", start, end)),
  row.names = FALSE,
)

# main function taking as input an increment from 1 to length(path_trees_phylo)
process_file <- function(i){
  
  output_file <- paste0(cluster_directory,
                        sprintf("/output/results/marginal_prob_first_split_mcc_consensus_%d_%d.csv", start, end))
  
  tryCatch({
    
    # load the mcc, the posterior and the deepest node of the outgroup of the true treee
    cs <- trees_cs[[i]]
    mcc <- trees_mcc[[i]]
    node_cs <- deepest_nodes_cs[i]
    node_mcc <- deepest_nodes_mcc[i]
    
    # load the paths, its simulation number and its age
    path <- path_trees_mcc[[i]]
    tree_simulation_number <- as.numeric(str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2])
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))
    
    # posterior probability of the outgroup for the consensus tree
    internal_nodes_cs <- (length(cs$tip.label) + 1):(length(cs$tip.label) + cs$Nnode)
    internal_nodes_mcc <- (length(mcc$tip.label) + 1):(length(mcc$tip.label) + mcc$Nnode)
    
    prob_cs <- cs$node.label[match(node_cs,internal_nodes_cs)]
    prob_mcc <- mcc$node.label[match(node_mcc,internal_nodes_mcc)]
    
    #plot(cs, show.node.label = T)
    #nodelabels(cex=0.6, frame='circle',  adj = c(0.5, 0.5))
    #plot(mcc, show.node.label = T)
    #nodelabels(cex=0.6, frame='circle',  adj = c(0.5, 0.5))
    
    # report age, simulation, the first split, and the clade frequencies of the first split 
    row <- data.frame(
      age = tree_age,
      simulation = tree_simulation_number,
      node_cs = node_cs,
      node_mcc = node_mcc,
      cs_prob = round(as.numeric(prob_cs), 3),
      mcc_prob = round(as.numeric(prob_mcc), 3)
    )
    
    # write the results
    write.table(
      row,
      file = output_file,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(output_file),
      append = TRUE,
      quote = FALSE
    )
    
    
  }, error = function(e) {})
}

# cluster initialisation
ncl <- detectCores() - 1
cl <- makeCluster(ncl, type="FORK")
clusterSetRNGStream(cl)


# Traitement en parallèle
res_list <- parLapply(cl, 1:length(path_trees_cs), process_file)
stopCluster(cl)
