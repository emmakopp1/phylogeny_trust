library(rwty)
library(here)
library(patchwork)
library(ape)
library(TreeTools)
library(tracerer)
library(phytools)
library(HDInterval)
library(ggtree)
library(tidyverse)
library(phangorn)
library(adephylo)


# Chargement des données
trees <- load.trees(
  file    = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/sino-tibet-covarion-fbd-strict.trees",
  logfile = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/sino-tibet-covarion-fbd-strict.log"
)

# --- originFBD ---
pp_param <- makeplot.param(trees, parameter = "TreeHeight.t.tree", burnin = 20)

# --- Topologie ---
pp_topo <- makeplot.topology(trees, burnin = 20)

# --- Plot trace of parameters ---
combined <- (pp_param$trace.plot | pp_topo$trace.plot) 
ggsave("/Users/kopp/Documents/these_overleaf/ch2/fig/convergence_mcmc.png", plot = combined, width = 12, dpi = 300)


# Compare summary trees
burnin=0.2
trees = read.nexus(file = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/sino-tibet-covarion-fbd-strict.trees")
trees = trees[ceiling(length(trees) * burnin):length(trees)]
log_tree = parse_beast_tracelog_file("/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/sino-tibet-covarion-fbd-strict.log")
log_tree = log_tree[ceiling(length(log_tree) * burnin):length(log_tree),]

# ------ Compute summary trees --------

# MAP
map_index = which.max(log_tree[,"posterior"])
tree_map = trees[[map_index]]

write.tree(
  tree_map,
  "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/map.tree"
)

# MRC 
consensus <- trees[ceiling(length(trees) * burnin):length(trees)]
consensus <- consensus(trees, p = .5, rooted = TRUE)
consensus <- consensus.edges(trees,
                             consensus.tree = consensus,
                             rooted = TRUE)

if (!is.rooted(consensus)) {
  consensus$root.edge.length <- 0
}

write.tree(
  consensus,
  "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/consensus.tree"
)


# ---- Plot summary trees ------
mcc = read.nexus(file = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/mcc.tree")
map = read.tree(file = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/map.tree")
hipstr = read.nexus(file = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/hipstr.tree")
consensus = read.tree( file = "/Users/kopp/Documents/phylogeny_trust_thesis_redaction/data/real/st_covarion-strict-fbd/consensus.tree")


# MCC
pdf("/Users/kopp/Documents/these_overleaf/ch2/fig/st_mcc_tree.pdf", width = 12, height = 8)
plot(mcc, 
     show.node.label = TRUE, 
     cex = 0.7,           
     no.margin = FALSE)
dev.off()

# MAP
pdf(here("/Users/kopp/Documents/these_overleaf/ch2/fig/st_map_tree.pdf"), width = 12, height = 8)
plot(map, 
     show.node.label = TRUE, 
     cex = 0.7,           
     no.margin = F)
dev.off()

# MRC
consensus$node.label = round(as.numeric(consensus$node.label),2)
pdf(here("/Users/kopp/Documents/these_overleaf/ch2/fig/st_consensus_tree.pdf"), width = 10, height = 8)
plot(consensus, 
     show.node.label = TRUE, 
     cex = 0.7,           
     no.margin = F)
dev.off()

# Hipstr
pdf(here("/Users/kopp/Documents/these_overleaf/ch2/fig/st_hipstr_tree.pdf"), width = 10, height = 8)
plot(hipstr, 
     show.node.label = TRUE, 
     cex = 0.7,           
     no.margin = F)
dev.off()

