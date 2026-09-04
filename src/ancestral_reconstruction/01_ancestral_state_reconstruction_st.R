# ==============================================================================
# SCRIPT: ancestral_reconstruction.R
# DESCRIPTION:
# This script performs ancestral state reconstruction on linguistic traits
# (cognates) for Indo-European and Sino-Tibetan languages using a posterior 
#distribution of trees inferred from BEAST.
#
# For each linguistic "meaning" (semantic category), and each sampled tree,
# the script identifies which internal nodes have a reconstructed probability
# of presence > 0.5 for each trait. It then computes the maximum depth
# (distance from the root) at which a trait can be reconstructed with high
# confidence.
#
# The final output is a dataset linking meanings to their maximum reconstruction
# depth, which can then be visualized or analyzed further.
#
# Main steps:
# - Load a set of BEAST trees (after burn-in and thinning)
# - Load linguistic data and meaning boundaries
# - For each meaning and each tree:
#     - Drop tips with missing data
#     - Reconstruct ancestral states using a Markov model
#     - Record the maximum depth for nodes with posterior prob > 0.5
# - Save results to CSV
#
# This files don't use the tidyverse library as it was not available on the 
# cluster.
# ==============================================================================

# ----- Load required libraries ------------------------------------------------
library(parallel)
library(here) 
library(ape)
library(phytools)
library(castor)
library(phangorn)
library(adephylo)

# arguments 
# for testing, this will lower the execution time 
#length_phylo <- 2
# to reproduce the exact results fix 200
length_phylo <- 200


# ----- Main processing function -----------------------------------------------
process_k <- function(k, param) {
  
  Y <- param$Y
  meanings_sets <- param$meanings_sets
  pi <- param$pi
  tree <- param$tree
  I_k <- param$I_k
  M <- param$M
  Q <- param$Q
  path_out <- param$path_out
  
  
  meaning_k <- meanings_sets$meaning[k]
  I <- I_k[k]
  start <- meanings_sets$start[k]
  end <- meanings_sets$end[k]
  Y_pruned <- t(Y[start:end, , drop = FALSE])
  
  
  for (t in 1:M) {
    tree_pruned <- tree[[t]]
    
    for (trait in 1:ncol(Y_pruned)) {
      y <- Y_pruned[, trait]
      ii <- which(is.na(y))
      
      if (length(ii) > 0) {
        y <- y[-ii]
        tree_pruned <- drop.tip(tree_pruned, names(ii))
        
      }

      if (sum(y, na.rm = TRUE) > 1 & length(unique(y)) > 1) {
        rec <- ancr(fitMk(tree_pruned, y, fixedQ = Q, "ARD", pi = as.numeric(pi)))
        indice <- which(rec[["ace"]][, 2] > 0.5)
        nodes <- as.integer(names(indice))
        
        if (length(nodes) > 0) {
            
          all_dists <- distRoot(tree_pruned, 1:max(tree_pruned$edge))
          values    <- as.numeric(max(distRoot(tree_pruned)) - all_dists)
          
          # Take the max only among the reconstructed nodes
          node  <- nodes[which.max(values[nodes])]
          value <- values[node]
          
          row <- data.frame(
            value = value,
            sens = meaning_k,
            tree = t,
            trait = trait + start - 1,
            node = node,
            root_age = max(all_dists),
            stringsAsFactors = FALSE
          )
          
          write.table(
            row,
            file = path_out,
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE
          )
        }
      }
    }
  }
}

# Sino-Tibetan  ---------------------------------------------------------------
# Load BEAST posterior trees and thin sample 
phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees"))
M_st <- length(phylo_st)
phylo_st <- phylo_st[seq(0.8 * M_st, M_st, length = length_phylo)]
M_st <- length(phylo_st)  # update M to number of retained trees

# Load and clean linguistic data 
Y_st <- read.nexus.data(here("data/real/st_ctmc-strict-fbd-uni/st.nex"))
Y_st <- lapply(Y_st, function(col) replace(col, col == "?", NA))
Y_st <- lapply(Y_st, as.numeric)
Y_st <- as.data.frame(Y_st)

# Load meaning boundaries and filter parameters
meanings_sets_st <- read.csv(here("output/results/meanings_sets_st.csv"))
bounds_real_tb_by_sens <- read.csv(here("output/results/tracelog_summary.csv"))
bounds_real_tb_st <- bounds_real_tb_by_sens[bounds_real_tb_by_sens$family == "ST", ]

# Define constants and substitution model 
K_st <- length(meanings_sets_st$meaning)
I_k_st <- meanings_sets_st$end - meanings_sets_st$start + 1

pi_st <- bounds_real_tb_st[, c("pi0", "pi1")]
clock_rate <-  0.018
lambda_st <- clock_rate / (2 * pi_st$pi0)
mu_st <- clock_rate / (2 * pi_st$pi1)
Q_st <- cbind(c(-lambda_st, mu_st), c(lambda_st, -mu_st))

# Initialize output CSV
path_out_st <- here("output/results/ancestral_reconstruction_st.csv")
write.csv(
  x = data.frame(
    value = character(),
    sens = character(),
    tree = character(), 
    trait = character(), 
    node = character(),
    root_age =  character(),
    stringsAsFactors = FALSE),
  file = path_out_st,
  row.names = FALSE
)

# Parameter for the main function 
param_st = list(tree = phylo_st, 
                Y = Y_st,
                meanings_sets = meanings_sets_st,
                I_k = I_k_st, 
                M = M_st, 
                pi = pi_st, 
                Q = Q_st, 
                path_out = path_out_st
)

# Execution for water -----------------------------------------------------
lapply(1:K_st, function(i) process_k(i, param_st))





