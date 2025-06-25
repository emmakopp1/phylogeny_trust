# ------------------------------------------------------------------------------
# Script Name: marginal_prob_first_split_ic.R
# Description: # This file calculate the marginale probabiliity that the inferred reconstructed
#               correctly the outgroup with the credibility interval
# -----------------------------------------------------------------------------------------
library(ape)
library(castor)
library(parallel)
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

# arguments 
args <- commandArgs(trailingOnly = TRUE)
start <- as.numeric(args[1])
end <- as.numeric(args[2])
cat("Traitement des fichiers de", start, "à", end, "\n")


# Fichiers
cluster_directory <- getwd()
path_phylo <- list.files(cluster_directory, full.names = TRUE, recursive = TRUE)

path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)][start:end]
#path_trees_phylo <- path_phylo[grepl("trees$", path_phylo)][(start + 24):(end + 24)] # pas sur le cluster
path_trees_phylo <- path_phylo[grepl("trees$", path_phylo)][start:end] # sur le cluster


# passer en mode parLapply(cl, path_trees_true, read.tree)
trees_true <- lapply(path_trees_true, read.tree)
phylogenies <- lapply(path_trees_phylo, read.nexus)
deepest_nodes <- lapply(trees_true, function(tree) getNodesByDepth(tree)[2])

burnin <- 0.1

# init of results
# dataframe of results 
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c("node", "tree_age", "tree_simulation_number", "prob_mean", "prob_inf", "prob_sup")

header_written <- FALSE


write.table(
  df,
  file = paste0(getwd(), sprintf("/marginal_prob_first_split_ic_%d_%d.csv", start, end)),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)




# Fonction de traitement parallèle
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

    # find the outgroup of the true tree
    A <-  tree_true$tip.label[Descendants(tree_true, deepest_node, type = "tips")[[1]]] 
    B <- setdiff(tree_true$tip.label, A)
    outgroup <- if (length(A) <= length(B)) A else B

    # posterior thin-in
    M <- length(phylo)
    phylo <- phylo[seq(burnin * M, M, length = 200)]


    res <- sapply(phylo,function(t) is.monophyletic(t, outgroup))
    ic <- prop.test(sum(res), length(res), conf.level = 0.95)$conf.int
    ic_inf <- ic[1]
    ic_sup <- ic[2]

    row <- data.frame(
        node = deepest_node,
        tree_age = tree_age,
        tree_simulation_number = tree_simulation_number,
        prob = mean(res),
        prob_inf = round(ic_inf,3),
        prob_sup = round(ic_sup,3)
    )
        

    write.table(
        row, 
        paste0(getwd(), sprintf("/marginal_prob_first_split_ic_%d_%d.csv", start, end)), 
        sep = ",",
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE     
    )

  return()
}

# cluster initialisation
ncl <- 20
cl <- makeCluster(ncl, type="FORK")
clusterSetRNGStream(cl)


# Traitement en parallèle
res_list <- parLapply(cl, 1:length(path_trees_phylo), process_file)
stopCluster(cl)
#res_list <- lapply(1:length(path_trees_phylo), function(i) process_file(i))

# results as a dataframe
res_df <- do.call(rbind, res_list) 
res_df <- as.data.frame(res_df)
colnames(res_df) <- c("node", "tree_age", "tree_simulation_number", "prob")


write.csv2(
  res_df,
  paste0(getwd(), sprintf("/marginal_prob_first_split_%d_%d.csv", start, end)),
  row.names = FALSE
)