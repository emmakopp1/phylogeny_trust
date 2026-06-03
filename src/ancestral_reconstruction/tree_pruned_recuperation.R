library(ape)
library(dplyr)
library(here)
library(phytools)
library(castor)
library(phangorn)
library(adephylo)
library(glue)

process_k_save_trees <- function(k, param) {
  
  Y <- param$Y
  meanings_sets <- param$meanings_sets
  pi <- param$pi
  tree <- param$tree
  I_k <- param$I_k
  M <- param$M
  Q <- param$Q
  path_out <- param$path_out
  df_ancr = param$df_ancr
  
  # Créer la clé une seule fois
  ancr_keys <- paste(df_ancr$tree, df_ancr$trait, sep = "_")
  
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
      
      trait_global <- trait + start - 1
      
      if (paste(t, trait_global, sep = "_") %in% ancr_keys) {
        saveRDS(
          tree_pruned,
          file = here(paste0(param$dir_out, "/tree_", meaning_k, "_t", t, "_trait", trait_global, ".rds"))
        )
      }
    }
  }
}

# Sino-Tibetan  ---------------------------------------------------------------
length_phylo <- 200
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

dir_out_st <- here("output/trees/ancr/tree_pruned_st/")

ancr_st = read.csv(here("output/results/ancestral_reconstruction_st.csv"))

# Parameter for the main function 
param_st = list(tree = phylo_st, 
                Y = Y_st,
                meanings_sets = meanings_sets_st,
                I_k = I_k_st, 
                M = M_st, 
                pi = pi_st, 
                Q = Q_st,
                dir_out = dir_out_st,
                df_ancr = ancr_st
)

dir.create(param_st$dir_out, showWarnings = FALSE, recursive = TRUE)

# Indo-Europeen ----------------------------------------------------------------
# Sino-Tibetan  ---------------------------------------------------------------

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

dir_out_ie <- here("output/trees/ancr/tree_pruned_ie")

# results from ancestral reconstruction
ancr_ie = read.csv(here("output/results/ancestral_reconstruction_ie.csv"))

# Parameter for the main function 
param_ie = list(tree = phylo_ie, 
                Y = Y_ie,
                meanings_sets = meanings_sets_ie,
                I_k = I_k_ie, 
                M = M_ie, 
                pi = pi_ie, 
                Q = Q_ie,
                dir_out = dir_out_ie,
                df_ancr = ancr_ie
)

dir.create(param_ie$dir_out, showWarnings = FALSE, recursive = TRUE)

# Execute code
# Sino-Tibétain 
for (k in 1:K_st) {
  process_k_save_trees(k, param_st)
}

# Indo-Européen
for (k in 1:K_ie) {
  process_k_save_trees(k, param_ie)
}

# Test 
length(list.files("/Users/kopp/Documents/phylogeny_trust/output/trees/ancr/tree_pruned_ie"))
dim(param_ie$df_ancr)[1]










