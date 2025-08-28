# ------------------------------------------------------------------------------
# Script Name: marginal_prob_first_split_ic.R
# Run : on the cluster
# Description: # This file calculate the marginale probabiliity that the inferred reconstructed
#               correctly the outgroup with the credibility interval
#               This file was run on a cluster is computationnaly costly. 
#               We strongly recommand to only make tests. 
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

# to be referenced by the user
phylo_length_test <- 100
# phylo_length_test <- 850

# path of the simulation folder you want to analyse
path_repository <- here("data/simulated-2025-05-13")

# compute the number of traits in the analysis
N_traits <- as.numeric(str_extract(path_repository, "\\d+$"))
# if N_traits is not 6 or 12 thousands, then it is the main study and N_traits = 3000
N_traits <- ifelse(N_traits %in% c(12000, 6000), N_traits, "")

# functions -------------------------------------------------------------------
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

# main function which compute for each tree, node, age and simulation the 
# the posterior of truthiness of the node by looking at the corresponding true tree 
# we also compute IC 
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
  
  # results
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
    file_path, 
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE     
  )
  
  return(row)
}

# path to access true trees and phylogenies
path_phylo <- list.files(path_repository, full.names = TRUE, recursive = TRUE)
path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)]
path_trees_phylo <- path_phylo[grepl("trees$", path_phylo)]

# load trees and phylogenies
trees_true <- lapply(path_trees_true[1:phylo_length_test], read.tree)
phylogenies <- lapply(path_trees_phylo[1:phylo_length_test], read.nexus)
deepest_nodes <- lapply(trees_true[1:phylo_length_test], function(tree) getNodesByDepth(tree)[2])

burnin <- 0.1

# dataframe of results 
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c("node", "tree_age", "tree_simulation_number", "prob_mean", "prob_inf", "prob_sup")


file_path <- ifelse(
  N_traits == "", 
  here("output/results/marginal_probability_first_split_ic.csv"), 
  here("output/results/marginal_probability_first_split_ic_%d.csv", N_traits))

write.csv(
  df, 
  file_path,
  row.names = FALSE,
)

# cluster initialisation
ncl <- 40
cl <- makeCluster(ncl, type="FORK")
clusterSetRNGStream(cl)


# parallelisation
res_list <- parLapply(cl, seq_along(path_trees_phylo[1:phylo_length_test]), process_file)
stopCluster(cl)

