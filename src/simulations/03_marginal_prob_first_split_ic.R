# ------------------------------------------------------------------------------
# Script Name: marginal_prob_first_split_ic.R
# Run on cluster
# Description: # This file calculate the marginale probabiliity that the inferred reconstructed
#               correctly the outgroup with the credibility interval
# -----------------------------------------------------------------------------------------
library(here)
library(ape)
library(phangorn)
library(Matrix)
library(castor)
library(gridExtra)
library(purrr)
library(phytools)
library(parallel)
library(adephylo)
library(reshape2)
library(stringr)
library(stats)

# function to get nodes by depths
getNodesByDepth <- function(tree) {
  # Recursive function
  recursiveTraversal <- function(node, result) {
    result[[length(result) + 1]] <- list(node = node, depth = distRoot(tree, node)[[1]])
    
    if (node %in% 1:(tree$Nnode + 1)) {
      return(result)
    } else {
      children <- tree$edge[tree$edge[, 1] == node, 2]
      
      for (child in children) {
        depth <- distRoot(tree, child)
        result <- recursiveTraversal(child, result)
      }
      return(result)
    }
  }
  
  nodes <- recursiveTraversal(castor::find_root(tree), list())
  nodes <- as.data.frame(do.call(rbind, nodes))
  nodes$depth <- unlist(nodes$depth)
  nodes$node <- unlist(nodes$node)
  return(nodes[order(-nodes$depth, decreasing = T), 1])
}


# load data
#cluster_directory <- gsub("simulation_analysis$", "simulated-2025-07-28", getwd()) # on the cluster 
cluster_directory <- here("data/simulated-2025-07-28") # in local 

# paths of the true tree and the phylogeny samples
path_phylo <- list.files(cluster_directory, full.names = TRUE, recursive = TRUE)
path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)]
path_trees_phylo <- path_phylo[grepl("trees$", path_phylo)]

# load trees and phylogenies 
trees_true <- lapply(path_trees_true, read.tree)
phylogenies <- lapply(path_trees_phylo, read.nexus)
deepest_nodes <- lapply(trees_true, function(tree) getNodesByDepth(tree)[2])

burnin <- 0.1

# data frame of results 
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c("node", "tree_age", "tree_simulation_number", "prob_mean", "prob_inf", "prob_sup")

output_path <- paste0(cluster_directory, "/marginal_prob_first_split_ic.csv")

write.table(
  df,
  file = output_path,
  sep = ",",
  col.names = TRUE
)


# parralel treatment
process_file <- function(i) {
  # true tree, its deepest node and posterior
  tree_true <- trees_true[[i]]
  phylo <- phylogenies[[i]]
  deepest_node <- deepest_nodes[[i]]
  
  # identify the files
  path <- path_trees_phylo[[i]]
  tree_simulation_number <- as.numeric(
    str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2]
  )
  tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))
  
  # find the smallest group of the first split of the true tree
  A <-  tree_true$tip.label[Descendants(tree_true, deepest_node, type = "tips")[[1]]] 
  B <- setdiff(tree_true$tip.label, A)
  outgroup <- if (length(A) <= length(B)) A else B
  
  # posterior thin-in
  M <- length(phylo)
  phylo <- phylo[seq(burnin * M, M, length = 200)]
  
  res <- sapply(phylo, function(t) is.monophyletic(t, outgroup))
  ic <- prop.test(sum(res), length(res), conf.level = 0.95)$conf.int
  ic_inf <- ic[1]
  ic_sup <- ic[2]
  
  row <- data.frame(
    node = deepest_node,
    tree_age = tree_age,
    tree_simulation_number = tree_simulation_number,
    prob_mean = mean(res),
    prob_inf = round(ic_inf, 3),
    prob_sup = round(ic_sup, 3)
  )
  
  
  write.table(
    row, 
    paste0(cluster_directory, "/marginal_prob_first_split_ic.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE     
  )
  
  return(row)
}

# cluster initialisation
ncl <- 40
cl <- makeCluster(ncl, type="FORK")
clusterSetRNGStream(cl)
res_list <- parLapply(cl, seq_along(path_trees_phylo), process_file)
stopCluster(cl)

