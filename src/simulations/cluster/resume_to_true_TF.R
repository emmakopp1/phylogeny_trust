# ------------------------------------------------------------------------------
# Script Name: resume_to_true_TF.R
# Description: #This file calculates whether a node present in the real tree is 
#                true, false or uncertain in the summary tree.
# -----------------------------------------------------------------------------------------

library(here)
library(ape)
library(phangorn)
library(Matrix)
library(castor)
library(here)
library(gridExtra)
library(purrr)
library(viridisLite)
library(viridis)
library(patchwork)
library(phytools)
library(adephylo)
library(stringr)
library(ggplot2)
library(reshape2)

# functions --------------------------------------------------------------------
# given a node (in the consensus tree), a true tree and a consensus tree, return if
# the node exists in the true tree
exist_node <- function(node, tree_true, tree_est) {
    x <- tree_est$tip.label[Descendants(tree_est, node)[[1]]]    
    return(list(res = is.monophyletic(tree_true, x), N_nodes =  tree_est$Nnode ))
}

# load data -------------
#simulation_folder <- "simulated-2025-07-08-12000"
#cluster_directory <- gsub("simulation_analysis$", simulation_folder, getwd()) # on the cluster

cluster_directory <- here("data/simulated-2025-07-08-12000")
path_phylo <- list.files(cluster_directory, full.names = TRUE, recursive = TRUE)

N_traits <- as.numeric(str_extract(cluster_directory, "\\d+$"))
# if N_traits is not 6 or 12 thousands, then it is the main study and N_traits = 3000
N_traits <- ifelse(N_traits %in% c(12000, 6000), N_traits, "")

path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)]
path_trees_cs <- path_phylo[grepl("consensus", path_phylo)]
path_trees_mcc <- path_phylo[grepl("mcc-", path_phylo)]

# passer en mode parLapply(cl, path_trees_true, read.tree)
trees_true <- lapply(path_trees_true, read.tree)
trees_cs <- lapply(path_trees_cs, read.tree)
trees_mcc <- lapply(path_trees_mcc, read.tree)


# dataframe of results 
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c('type', 'age', 'simulation', 'node', 'exist', 'N_nodes')


# initialisation of the file
write.csv(
    df, 
    paste0(cluster_directory, sprintf("/resume_to_true_TF_%d.csv", N_traits)),
    row.names = FALSE,
    )


# for the consensus trees
for (t in seq_along(trees_true)){
    # consensus and true tree
    tt = trees_true[[t]]
    cs = trees_cs[[t]]

    path <- path_trees_true[[t]]
    tree_simulation_number <- as.numeric(
      str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2]
    )
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))

    N_tip = length(cs$tip.label)

    for (node in seq(N_tip + 1, N_tip + cs$Nnode )){
        exist <- exist_node(node, tt, cs)
        row <- data.frame(
            type = 'consensus',
            age = tree_age,
            simulation = tree_simulation_number,
            node = node,
            state = exist$N_nodes,
            result = exist$res
        )
        
        # Écrire avec ou sans en-têtes selon si c'est la première fois
        write.table(
            row, 
            paste0(cluster_directory, sprintf("/resume_to_true_TF_%d.csv", N_traits)), 
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE       # Ajouter après la première fois
        )
    }
}


# for mcc tree 
for (t in seq_along(trees_true)){
    # consensus and true tree
    tt = trees_true[[t]]
    mcc = trees_mcc[[t]]

    path <- path_trees_true[[t]]
    tree_simulation_number <- as.numeric(
      str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2]
    )
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))

    N_tip = length(mcc$tip.label)
    
    for (node in seq(N_tip + 1, N_tip + mcc$Nnode)){
        exist <- exist_node(node, tt, mcc)
        row <- data.frame(
            type = 'mcc',
            age = tree_age,
            simulation = tree_simulation_number,
            node = node,
            state = exist$N_nodes,
            result = exist$res
        )
        
        # Écrire avec ou sans en-têtes selon si c'est la première fois
        write.table(
            row, 
            paste0(cluster_directory, sprintf("/resume_to_true_TF_%d.csv", N_traits)), 
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE       # Ajouter après la première fois
        )
    }
}
