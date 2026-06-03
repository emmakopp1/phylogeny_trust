# ------------------------------------------------------------------------------
# Script Name: visualization.R
# Run : local
# Description: Creates heatmap visualizations for shared cognate analysis results
# -----------------------------------------------------------------------------------------

library(here)
library(pheatmap)
library(tidyverse)
library(phytools); library(xml2); library(ape); library(phangorn)
library(castor); library(adephylo)

# Partie 1 - Reel données proportion de cognats partagé entre les deux sous arbre 
# Données réel et arbre consensus 

tree_cs <- read.tree(here("data/real/st_ctmc-strict-fbd-uni/st_consensus.tree"))
dat     <- read_xml(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.xml"))

# Sequences ordonnées selon les tip labels
seq_nodes <- xml_find_all(dat, ".//sequence")
seqs <- setNames(xml_attr(seq_nodes, "value"), xml_attr(seq_nodes, "taxon"))
seqs <- seqs[tree_cs$tip.label]

# Tips des deux sous-arbres enfants de la racine
root_kids <- tree_cs$edge[tree_cs$edge[,1] == Ntip(tree_cs) + 1, 2]
get_tips  <- \(n) tree_cs$tip.label[Descendants(tree_cs, n, type = "tips")[[1]]]

# Proportion de sites partagés (colonne "1" présente dans les deux sous-arbres)
site_matrix <- \(tips) do.call(rbind, strsplit(seqs[tips], ""))
has_one     <- \(tips) colSums(site_matrix(tips) == "1") > 0

shared <- has_one(get_tips(root_kids[1])) & has_one(get_tips(root_kids[2]))
mean(shared)

# load data --------------------------------------------------------------------
shared_cognates <- readRDS(here("output/results/shared_cognates.rds"))
prop_shared_tip_pair <- readRDS(here("output/results/shared_cognate_tip_pair.pdf"))




