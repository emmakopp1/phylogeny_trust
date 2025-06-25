# ------------------------------------------------------------------------------
# Script Name: true_false_uncertain.R
# Description: Checks whether a node present in the true tree is also found 
#              (as true, false, or uncertain) in the summary/consensus tree.
# -----------------------------------------------------------------------------------------

library(ape)
library(stringr)

# functions --------------------------------------------------------------------
# given a node (in the true tree), a true tree and a consensus tree, return if
# the node is plausible in the consensus tree. 

# to compute that we compare for 

is_plausible <- function(node, tree_true, tree_cs) {
    # descendant of node in the true tree
    x <- tree_true$tip.label[Descendants(tree_true, node)[[1]]]
    # mrca of x in the consensus tree
    mrca_cs <- mrca.phylo(tree_cs, x)
    # descendant of the mrca in the consensus tree
    desc_cs <- tree_cs$edge[tree_cs$edge[, 1] == mrca_cs, 2]

    # if its a regular node with two childrens
    if (length(desc_cs) == 2) {
        return(list(state = "regular", res = is.monophyletic(tree_cs, x)))
        } 
    # if the node is of a rake shape
    else {
        # if some descendant are tip, change their label to the tip label
        desc_cs_updated <- desc_cs
        for (i in seq_along(desc_cs)) {
            desc <- desc_cs[i]
            if (desc %in% 1:length(tree_cs$tip.label)) {
                desc_cs_updated[i] <- tree_cs$tip.label[desc]
            }
        }
    desc_cs <- desc_cs_updated

    plausible_node <- c()
    for (i in seq_along(desc_cs)) {
      desc <- desc_cs[i]
      # condition : if the descendant in the consensus tree is a node
      if (sum(grep("t", desc)) == 0) {
        # descendant of the descendant in the consensus tree
        desc_desc <- tree_cs$tip.label[Descendants(tree_cs, as.numeric(desc))[[1]]]
      }
      # if the descendant in the consensus tree is a tip
      else if ((sum(grep("t", desc)) == 1)) {
        # descendant of the descendant in the consensus tree
        desc_desc <- tree_cs$tip.label[Descendants(tree_cs, desc)[[1]]]
      }
    
    # condition1 : are the descendants in the consensus and the true tree the same ?
    condition1 <- length(intersect(desc_desc, x)) == length(desc_desc)
    # condition2 : is the desncendant set empty ? 
    condition2 <- is_empty(intersect(desc_desc, x))

    # if condition1 or condition2 is true the naude is plausible
    plausible_node <- c(plausible_node, condition1 || condition2)
    }

    return(list(state = "rateau", res = all(plausible_node)))
  }
}

# arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
start <- as.numeric(args[1])
end <- as.numeric(args[2])
cat("Traitement des fichiers de", start, "à", end, "\n")


# load data --------------------------------------------------------------------
# paths 
cluster_directory <- getwd()
path_phylo <- list.files(cluster_directory, full.names = TRUE, recursive = TRUE)
path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)][start:end]
path_trees_cs <- path_phylo[grepl("consensus", path_phylo)][start:end]

# load the true and consensus trees
trees_true <- lapply(path_trees_true, read.tree)
trees_cs <- lapply(path_trees_cs, read.tree)

# dataframe of results 
df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c('tree_age', 'tree_simulation_number', 'node', 'state', 'T_F_U')

write.csv(
    df, 
    paste0(getwd(), sprintf("/true_false_uncertain_nodes_%d_%d.csv", start, end)),
    #"/home/users/kopp/work/simulated-2025-05-13/true_false_uncertain_nodes.csv",
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,  # En-têtes seulement la première fois
    )


for (t in seq_along(trees_true)){
    # consensus and true tree
    tt <- trees_true[[t]]
    cs <- trees_cs[[t]]

    #load the paths, its simulation number and its age
    path <- path_trees_true[[t]]
    tree_simulation_number <- as.numeric(
      str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2]
    )
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))

     
    for (node in seq(tt$Nnode + 2, 2 * tt$Nnode + 1)) {
        plausible <- is_plausible(node, tt, cs)
        row <- data.frame(
            tree_age = tree_age,
            tree_simulation_number = tree_simulation_number,
            node = node,
            state = plausible$state,
            result = plausible$res
        )
        
        # Écrire avec ou sans en-têtes selon si c'est la première fois
        write.table(
            row, 
            paste0(getwd(), sprintf("/true_false_uncertain_nodes_%d_%d.csv", start, end)), 
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE       # Ajouter après la première fois
        )
    }
}



