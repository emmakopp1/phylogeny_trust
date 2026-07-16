library(here)
library(phangorn)
library(ape)
library(TreeDist)

cluster_directory <- here("data/simulated-2025-07-28") 

# paths of the true tree and the phylogeny samples
path_phylo <- list.files(cluster_directory, full.names = TRUE, recursive = TRUE)
path_trees_true <- path_phylo[grepl("tree-sim", path_phylo)]
path_trees_phylo <- path_phylo[grepl("trees$", path_phylo)]

# load trees and phylogenies 
trees_true <- lapply(path_trees_true, read.tree)
phylogenies <- lapply(path_trees_phylo, read.nexus)

burnin <- 0.1

# data frame of results 
df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c("tree_age", "tree_simulation_number", "RF-mean", "RF-median", "RF-min", "RF-max", "RF-sd")

#output_path <- paste0(cluster_directory, "/rf_values.csv")
output_path <- here("output/results/rf_values_simu.csv")

write.table(
  df,
  file = output_path,
  sep = ",",
  col.names = TRUE
)


process_file <- function(i) {
  # true tree, its deepest node and posterior
  tree_true <- trees_true[[i]]
  phylo <- phylogenies[[i]]
  
  # identify the files
  path <- path_trees_phylo[[i]]
  tree_simulation_number <- as.numeric(
    str_match(path, "beast-data-sim-(\\d+)-\\d+")[, 2]
  )
  tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))
  
  # posterior thin-in
  M <- length(phylo)
  phylo <- phylo[seq(burnin * M, M, length = 200)]
  
  # RF- computations
  n <- Ntip(tree_true)
  rf_max <- 2 * (n - 1)
  rf_raw_all <- RF.dist(tree_true, phylo, rooted = TRUE)
  res <- rf_raw_all / rf_max
  
  row <- data.frame(
    tree_age = tree_age,
    tree_simulation_number = tree_simulation_number,
    RF-mean = round(mean(res), 3),
    RF-median = round(median(res), 3),
    RF-min = round(range(res), 3)[1], 
    RF-max = round(range(res), 3)[2],
    RF-sd = round(sd(res), 3)
  )
  
  
  write.table(
    row, 
    here("output/results/rf_values_simu.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE     
  )
  
  return(row)
}

# test
tree_true[[1]]
#lapply(cl, seq_along(path_trees_phylo[1:17]), process_file)
