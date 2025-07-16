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
length_phylo <- 2 
# to reproduce the exact results fix 200
#length_phylo <- 200


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
        rec <- ancr(fitMk(tree_pruned, y, "ARD", fittedQ = Q, pi = as.numeric(pi)))
        indice <- which(rec[["ace"]][, 2] > 0.5)
        nodes <- as.integer(names(indice))
        
        if (length(nodes) > 0) {
          value <- max(as.numeric(max(distRoot(tree_pruned)) - distRoot(tree_pruned, nodes)))
          
          row <- data.frame(
            value = value,
            sens = meaning_k,
            tree = t,
            trait = trait + start - 1,
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

# Indo-European  ---------------------------------------------------------------
# Load BEAST posterior trees and thin sample 
phylo_ie <- read.nexus(here("data/real/iecor_ctmc-strict-M1/IECoR_M1_CTMC_Gamma_1_Rate_For_All_Mgs_combined.trees"))
M_ie <- length(phylo_ie)
phylo_ie <- phylo_ie[seq(0.8 * M_ie, M_ie, length = length_phylo)]
M_ie <- length(phylo_ie)  # update M to number of retained trees

# Load and clean linguistic data 
Y_ie <- read.nexus.data(here("data/real/iecor_ctmc-strict-M1/iecor.nex"))
Y_ie <- lapply(Y_ie, function(col) replace(col, col == "?", NA))
Y_ie <- lapply(Y_ie, as.numeric)
Y_ie <- as.data.frame(Y_ie)

# Load meaning boundaries and filter parameters
meanings_sets_ie <- read.csv(here("output/results/meanings_sets_ie.csv"))
bounds_real_tb_by_sens_ie <- read.csv(here("output/results/tracelog_summary.csv"))
bounds_real_tb_ie <- bounds_real_tb_by_sens_ie[bounds_real_tb_by_sens_ie$family == "IE", ]

# Define constants and substitution model 
K_ie <- length(meanings_sets_ie$meaning)
I_k_ie <- meanings_sets_ie$end - meanings_sets_ie$start + 1

pi_ie <- bounds_real_tb_ie[, c("pi0", "pi1")]
lambda_ie <- 1 / (2 * pi_ie$pi0)
mu_ie <- 1 / (2 * pi_ie$pi1)
Q_ie <- cbind(c(-lambda_ie, mu_ie), c(lambda_ie, -mu_ie))

# Initialize output CSV
path_out_ie <- here("output/results/ancestral_reconstruction_ie.csv")
write.csv(
  x = data.frame(
    value = character(), 
    sens = character(), 
    tree = character(), 
    trait = character(),
    stringsAsFactors = FALSE),
  file = path_out_ie,
  row.names = FALSE
)

# Parameter for the main function 
param_ie = list(tree = phylo_ie, 
                Y = Y_ie,
                meanings_sets = meanings_sets_ie,
                I_k = I_k_ie, 
                M = M_ie, 
                pi = pi_ie, 
                Q = Q_ie, 
                path_out = path_out_ie
)


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
lambda_st <- 1 / (2 * pi_st$pi0)
mu_st <- 1 / (2 * pi_st$pi1)
Q_st <- cbind(c(-lambda_st, mu_st), c(lambda_st, -mu_st))

# Initialize output CSV
path_out_st <- here("output/results/ancestral_reconstruction_st.csv")
write.csv(
  x = data.frame(
    value = character(),
    sens = character(),
    tree = character(), 
    trait = character(), 
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



# ----- Parallel execution -----------------------------------------------------
ncl <- detectCores() - 1
cl <- makeCluster(ncl, type = "FORK")
clusterSetRNGStream(cl)

# Run the rest in parallel
clusterExport(cl, varlist = c(
  "process_k", "param_st", "param_ie"
))

parLapply(cl, 1:K_st, function(k) {
  process_k(k, param_st)
})


parLapply(cl, 1:K_ie, function(k) {
  process_k(k, param_ie)
})

stopCluster(cl)
