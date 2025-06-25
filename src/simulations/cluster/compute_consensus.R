# ------------------------------------------------------------------------------
# Script Name: compute_consensus.R
# Description: compute the consensus tree for all ages, simulations.
# -----------------------------------------------------------------------------------------
library(here)
library(ape)
library(phangorn)
library(purrr)
library(stringr)

# start and end
#args <- commandArgs(trailingOnly = TRUE)
#start <- as.numeric(args[1])
#end <- as.numeric(args[2])

#cat("Traitement des fichiers de", start, "à", end, "\n")

# functions --------------------------------------------------------------------
# function which take a path to a posterior phylogeny and write the consensus tree
process_phylo <- function(path) {
    burnin = 0.1
    phylogeny <- read.nexus(path)
  
    # burnin and thin-in
    M <- length(phylogeny)
    phylogeny <- phylogeny[seq(burnin * M, M, length = 200)] # thin-in
  
    # compute consensus
    phylogeny_cs <- consensus(phylogeny, p = 0.5, rooted = TRUE)
    phylogeny_cs <- consensus.edges(phylogeny,
        consensus.tree = phylogeny_cs,
        rooted = TRUE)
  
  
    # check that its rooted
    if (!is.rooted(phylogeny_cs)) {
        phylogeny_cs$root.edge.length <- 0
        }
  
    # identification of the age of the tree
    tree_age <- as.numeric(str_extract(path, "(\\d+)(?=\\.tree)"))
    path_cs <- str_replace(path, "ctmc-strict-bd-(\\d+)\\.trees", paste0("consensus-", tree_age, ".tree"))
    
    # write the tree
    write.tree(phylogeny_cs, path_cs)
    return(invisible(NULL))
}

# prepare data
path_phylo <- list.files(getwd(), full.names = TRUE, recursive = TRUE)
path_phylo <- path_phylo[grepl("\\.trees$", path_phylo)]

# application
purrr::map(path_phylo, process_phylo)

