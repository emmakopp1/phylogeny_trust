# ------------------------------------------------------------------------------
# Script Name: resume_to_true_TF.R
# Description: #This file calculates whether a node present in the real tree is 
#                true, false or uncertain in the summary tree.
# -----------------------------------------------------------------------------------------

library(ape)
library(castor)
library(stringr)

# functions --------------------------------------------------------------------
# given a node (in the consensus tree), a true tree and a consensus tree, return if
# the node exists in the true tree
exist_node <- function(node, tree_true, tree_est) {
    x <- tree_est$tip.label[Descendants(tree_est, node)[[1]]]    
    return(list(res = is.monophyletic(tree_true, x), N_nodes =  tree_est$Nnode ))
}


# arguments 
args <- commandArgs(trailingOnly = TRUE)
start <- as.numeric(args[1])
end <- as.numeric(args[2])
cat("Traitement des fichiers de", start, "à", end, "\n")


# load data -------------
cluster_directory <- getwd()
path_phylo <- list.files(cluster_directory, full.names = TRUE, recursive = TRUE)

path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)][start:end]
path_trees_cs <- path_phylo[grepl("consensus", path_phylo)][start:end]
path_trees_mcc <- path_phylo[grepl("mcc", path_phylo)][start:end]

# passer en mode parLapply(cl, path_trees_true, read.tree)
trees_true <- lapply(path_trees_true, read.tree)
trees_cs <- lapply(path_trees_cs, read.tree)
trees_mcc <- lapply(path_trees_mcc, read.tree)


# dataframe of results 
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c('type', 'age', 'simulation', 'node', 'exist', 'N_nodes')

header_written <- FALSE

# initialisation of the file
write.csv(
    df, 
    paste0(getwd(), sprintf("/resume_to_tue_TF_%d_%d.csv", start, end)),
    #"/home/users/kopp/work/simulated-2025-05-13/true_false_uncertain_nodes.csv",
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,  # En-têtes seulement la première fois
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
            paste0(getwd(), sprintf("/resume_to_true_TF_%d_%d.csv", start, end)), 
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
            paste0(getwd(), sprintf("/resume_to_tue_TF_%d_%d.csv", start, end)), 
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE       # Ajouter après la première fois
        )
    }
}
