# ------------------------------------------------------------------------------
# Script Name: compute_mcc.R
# Description: compute the Maximum-Clade-Credibility (MCC) for all ages, simulations.
#             Change simulation_folder to choose another simulation study
# -----------------------------------------------------------------------------------------
library(here)
library(ape)
library(phangorn)
library(Matrix)
library(castor)
library(here)
library(gridExtra)
library(purrr)
library(phytools)
library(adephylo)
library(reshape2)
library(stringr)

# to be referenced by the user
# path of the simulation folder you want to analyse
simulation_folder <- here("data/simulated-2025-05-13")
#simulation_folder <- here("data/simulated-2025-07-02-6000")
#simulation_folder <- here("data/simulated-2025-07-08-12000")

# functions --------------------------------------------------------------------
# function which take a path to a posterior phylogeny and write the mcc tree
compute_mcc_phylo <- function(path) {
    burnin = 0.1
    phylogeny <- read.nexus(path)
  
    # burnin and thin-in
    M <- length(phylogeny)
    phylogeny <- phylogeny[seq(burnin * M, M, length = 200)] 

    # compute mcc 
    tree_mcc <- maxCladeCred(phylogeny)

    # identification of the age of the tree
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))
    path_mcc <- str_replace(path, "ctmc-strict-bd-(\\d+)\\.trees", paste0("mcc-", tree_age, ".tree"))

    # write the tree
    write.tree(tree_mcc, path_mcc)
    return(invisible(NULL))
}

# prepare data
path_phylo <- list.files(simulation_folder, full.names = TRUE, recursive = TRUE)
path_phylo <- path_phylo[grepl("\\.trees$", path_phylo)]

purrr::map(path_phylo, compute_mcc_phylo)

