# ------------------------------------------------------------------------------
# Script Name: true_false_uncertain.R
# Run : on the cluster
# Description: Checks whether a node present in the true tree is also found 
#              (as true, false, or uncertain) in the summary/consensus tree.
# -----------------------------------------------------------------------------------------

library(here)
library(ape)
library(phangorn)
library(Matrix)
library(castor)
library(here)
library(phytools)
library(adephylo)
library(stringr)
library(reshape2)

# to be referenced by the user -------------------------------------------------
# path of the simulation folder you want to analyse
#path_repository <- here("data/simulated-2025-07-22-1500")
#path_repository <- here("data/simulated-2025-07-28")
#path_repository <- here("data/simulated-2025-07-22-6000")
path_repository <- here("data/simulated-2025-07-22-12000")

# set the number of traits 
# if N_traits is not 6 or 12 thousands, then it is the main study and N_traits = 3000
N_traits <- as.numeric(str_extract(path_repository, "\\d+$"))
N_traits <- ifelse(N_traits %in% c(12000, 6000, 1500), N_traits, "")


# functions --------------------------------------------------------------------

# given a node (in the true tree), a true tree and a consensus tree, return if
# the node is plausible in the consensus tree. 
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

# exclude calibration nodes
get_excluded_nodes <- function(tree, tips) {
  mrca <- getMRCA(tree, tips)
  desc <- Descendants(tree, mrca, type = 'all')
  setdiff(c(mrca, desc), seq_len(Ntip(tree)))
}

# load data --------------------------------------------------------------------
path_phylo <- list.files(path_repository, full.names = TRUE, recursive = TRUE)
path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)]
path_trees_cs <- path_phylo[grepl("consensus-", path_phylo)]

# load the true and consensus trees
trees_true <- lapply(path_trees_true, read.tree)
trees_cs <- lapply(path_trees_cs, read.tree)

# dataframe of results 
df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c('tree_age', 'tree_simulation_number', 'node', 'state', 'T_F_U')

file_path <- ifelse(
  N_traits == "", 
  here("output/results/true_false_uncertain_nodes.csv"), 
  here(sprintf("output/results/true_false_uncertain_nodes_%d.csv", N_traits)))

write.csv(
  df, 
  file_path,
  row.names = FALSE,
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

  # exclude calibration nodes 
  calib_chinese <- tt$tip.label[grep('Sinitic', tt$tip.label)]
  calib_tibetan <- tt$tip.label[grep('Tibetan', tt$tip.label)]
  calib_burmish <- c("BurmishOldBurmese", "BurmishRangoon")
  
  nodes_to_exclude <- unlist(lapply(
    list(calib_chinese, calib_tibetan, calib_burmish),
    get_excluded_nodes,
    tree = tt
  ))
  
  node_for_loop = setdiff(seq(tt$Nnode + 2, 2 * tt$Nnode + 1),nodes_to_exclude)
  
  for (node in node_for_loop){
    plausible <- is_plausible(node, tt, cs)
    row <- data.frame(
      tree_age = tree_age,
      tree_simulation_number = tree_simulation_number,
      node = node,
      state = plausible$state,
      result = plausible$res
    )
    
    
    write.table(
      row, 
      file_path, 
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE      
    )
  }
}



