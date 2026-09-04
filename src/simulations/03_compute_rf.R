library(here)
library(phangorn)
library(ape)
library(here)
library(ape)
library(phangorn)
library(castor)
library(parallel)
library(stringr)
library(stats)

# to be referenced by the user
# path of the simulation folder you want to analyse
cluster_directory <- here("data/simulated-2025-07-28") 
#cluster_directory <- here("data/simulated-2025-07-22-1500")
#cluster_directory <- here("data/simulated-2025-07-22-6000")
#cluster_directory <- here("data/simulated-2025-07-22-12000")

# set the number of traits 
# if N_traits is not 6 or 12 thousands, then it is the main study and N_traits = 3000
N_traits <- as.numeric(str_extract(cluster_directory, "\\d+$"))
N_traits <- ifelse(N_traits %in% c(12000, 6000, 1500), N_traits, "")


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
colnames(df) <- c('tree_age', 'tree_simulation_number', 'RF_mean', 'RF_median', 'RF_min', 'RF_max', 'RF_sd')


output_path <- ifelse(
  N_traits == "",
  paste0(here("output/results"), "/rf_values.csv"),
  paste0(here("output/results"), sprintf("/rf_values_%s.csv", N_traits))
)

write.table(
  df,
  file = output_path,
  sep = ",",
  col.names = TRUE
)


rf <- function(i) {
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
    RF_mean = round(mean(res), 3),
    RF_median = round(median(res), 3),
    RF_min = round(range(res), 3)[1], 
    RF_max = round(range(res), 3)[2],
    RF_sd = round(sd(res), 3)
  )
  
  
  write.table(
    row, 
    output_path,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE     
  )
  
  return()
}

lapply(seq_along(path_trees_phylo), rf)
