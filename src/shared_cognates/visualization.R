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

# paramètres de simulations 
N_traits   <- 3000
pi1        <- 0.94305
pi0        <- 0.05695
clock_rate <- 0.018
Q <- matrix(
  c(-clock_rate/(2*pi0),  clock_rate/(2*pi0),
    clock_rate/(2*pi1), -clock_rate/(2*pi1)),
  nrow = 2, byrow = TRUE
)

# load data --------------------------------------------------------------------
shared_cognates <- readRDS(here("output/results/shared_cognates.rds"))
prop_shared_tip_pair <- readRDS(here("output/results/shared_cognate_tip_pair.pdf"))

plot(shared_cognates$tree_age, shared_cognates$value,
     type = "b",
     pch  = 16,
     ylim = range(c(shared_cognates$inf, shared_cognates$sup)),
     xlab = "Tree age",
     ylab = "Proportion of shared cognates",
     main = "Shared cognates with 95% confidence interval")
arrows(x0     = shared_cognates$tree_age,
       y0     = shared_cognates$inf,
       y1     = shared_cognates$sup,
       angle  = 90,
       code   = 3,
       length = 0.05,
       col    = "steelblue")

# cognats partagé entre une pair de langue (une de l'outgroup et une de l'ingroup)
prop_shared_tip_pair <- readRDS(here("output/results/shared_cognate_tip_pair.pdf"))
tip_A <- unique(prop_shared_tip_pair$tipA)
tip_B <- unique(prop_shared_tip_pair$tipB)
plot(1:17, prop_shared_tip_pair$prop,
     pch  = 16,
     col  = "steelblue",
     xlab = "Âge (millénaires)",
     ylab = "Proportion de cognats partagés",
     main = paste0("Paire : ", tip_A, " — ", tip_B),
     ylim = c(0, 1))
# Courbe théorique de Swadeash
curve(exp(-2 * Q[1,2] * x), from = 1, to = 17, add = TRUE, col = "red", lty = 2, lwd = 2)
