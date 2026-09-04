# Libraries 
library(here) 
library(ape)
library(phytools)
library(castor)
library(phangorn)
library(adephylo)

# Import data ------------------------------------------------------------------
length_phylo <- 200
# Load BEAST posterior trees and thin sample 
phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees"))
M_st <- length(phylo_st)
phylo_st <- phylo_st[seq(0.8 * M_st, M_st, length = length_phylo)]
M_st <- length(phylo_st)
root_age <- mean(sapply(phylo_st, function(tree) max(distRoot(tree))))

# Load meaning boundaries and filter parameters
meanings_sets_st <- read.csv(here("output/results/meanings_sets_st.csv"))
bounds_real_tb_by_sens <- read.csv(here("output/results/tracelog_summary.csv"))
bounds_real_tb_st <- bounds_real_tb_by_sens[bounds_real_tb_by_sens$family == "ST", ]

# Functions --------------------------------------------------------------------
S <- function(node, tree, mu_st){
  
  if (node %in% 1:Ntip(tree)){return(1)}
  else{
    children = Descendants(tree, node, type = "children")
    branch_length_chidren_index = sapply(children, function(i) which(tree$edge[,1] == node & tree$edge[,2] == i))
    branch_length_children = tree$edge.length[branch_length_chidren_index]
    exp_branch_length_children = exp(-mu_st * branch_length_children)
    
    return(
      1 - (1 - exp_branch_length_children[1] + exp_branch_length_children[1]*(1 - S(children[1], tree, mu_st)))
      * (1 - exp_branch_length_children[2] + exp_branch_length_children[2]*(1 - S(children[2], tree, mu_st)))
    )
  }
}

T_1 <- function(node, tr, mu_st, lambda_st) {
  if (node %in% 1:Ntip(tr)) { return(1) }
  
  children <- Descendants(tr, node, type = "children")
  idx      <- sapply(children, function(i) which(tr$edge[,1] == node & tr$edge[,2] == i))
  bl       <- tr$edge.length[idx]
  
  p_loss_1 <- (1 - exp(-mu_st * bl[1])) * (1 - T_1(children[1], tr, mu_st, lambda_st)) +
    exp(-mu_st * bl[1])       * (1 - T_0(children[1], tr, mu_st, lambda_st))
  
  p_loss_2 <- (1 - exp(-mu_st * bl[2])) * (1 - T_1(children[2], tr, mu_st, lambda_st)) +
    exp(-mu_st * bl[2])       * (1 - T_0(children[2], tr, mu_st, lambda_st))
  
  return(1 - p_loss_1 * p_loss_2)
}

T_0 <- function(node, tr, mu_st, lambda_st) {
  if (node %in% 1:Ntip(tr)) { return(1) }
  
  children <- Descendants(tr, node, type = "children")
  idx      <- sapply(children, function(i) which(tr$edge[,1] == node & tr$edge[,2] == i))
  bl       <- tr$edge.length[idx]
  
  p_loss_1 <- (1 - exp(-lambda_st * bl[1])) * (1 - T_0(children[1], tr, mu_st, lambda_st)) +
    exp(-lambda_st * bl[1])        * (1 - T_1(children[1], tr, mu_st, lambda_st))
  
  p_loss_2 <- (1 - exp(-lambda_st * bl[2])) * (1 - T_0(children[2], tr, mu_st, lambda_st)) +
    exp(-lambda_st * bl[2])        * (1 - T_1(children[2], tr, mu_st, lambda_st))
  
  return(1 - p_loss_1 * p_loss_2)
}

# Compute parameters -----------------------------------------------------------
pi_st <- bounds_real_tb_st[, c("pi0", "pi1")]
clock_rate_st <- 0.018
lambda_st <- clock_rate_st / (2 * pi_st$pi0)
mu_st     <- clock_rate_st / (2 * pi_st$pi1)

# Compute S_root over all trees ------------------------------------------------
S_root_all <- numeric(M_st)

for (i in seq_len(M_st)) {
  tree <- phylo_st[[i]]
  
  root <- find_root(tree)
  
  first_split <- Descendants(tree, root, type = "children")
  branch_length_chidren_index <- sapply(first_split, function(j) which(tree$edge[,1] == root & tree$edge[,2] == j))
  branch_length_children      <- tree$edge.length[branch_length_chidren_index]
  exp_branch_length_children  <- exp(-mu_st * branch_length_children)
  
  S_root_all[i] <- ( exp_branch_length_children[1]) * S(first_split[1], tree, mu_st) *
    ( exp_branch_length_children[2]) * S(first_split[2], tree, mu_st)
}

# Summary
summary(S_root_all)
hist(S_root_all, main = "Distribution de S_root sur les arbres postérieurs",
     xlab = "S_root", col = "steelblue", border = "white")
mean(S_root_all)

# Apply T(root)
T_1(first_split[1], tree, mu_st, lambda_st) * (1 - exp(-mu_st * branch_length_children[1])) + 
  T_1(first_split[2], tree, mu_st, lambda_st) * (1 - exp(-mu_st * branch_length_children[2])) - 
  T_1(first_split[1], tree, mu_st, lambda_st) * (1 - exp(-mu_st * branch_length_children[1])) *
  T_1(first_split[2], tree, mu_st, lambda_st) * (1 - exp(-mu_st * branch_length_children[2]))


# Sensitivity of S(root) to tree age -------------------------------------------
target_ages <- seq(1, 17, by = 0.1)

S_matrix <- matrix(NA, nrow = M_st, ncol = length(target_ages))

for (k in seq_len(M_st)) {
  tree_k     <- phylo_st[[k]]
  root_k     <- find_root(tree_k)
  root_age_k <- max(distRoot(tree_k))
  
  for (i in seq_along(target_ages)) {
    
    coef_k <- target_ages[i] / root_age_k
    
    tree_scaled             <- tree_k
    tree_scaled$edge.length <- tree_k$edge.length * coef_k
    
    first_split_k <- Descendants(tree_scaled, root_k, type = "children")
    idx_root_k    <- sapply(first_split_k, function(j)
      which(tree_scaled$edge[,1] == root_k & tree_scaled$edge[,2] == j))
    bl_root_k     <- tree_scaled$edge.length[idx_root_k]
    
    S_matrix[k, i] <- (exp(-mu_st * bl_root_k[1])) * S(first_split_k[1], tree_scaled, mu_st) *
      ( exp(-mu_st * bl_root_k[2])) * S(first_split_k[2], tree_scaled, mu_st)
  }
}

# Mean over all trees, vector of length: length(target_ages)
S_mean <- colMeans(S_matrix)

results_S <- data.frame(
  tree_age = target_ages,
  S_root   = S_mean
)

# Write files
write.csv(results_S, here("output/results/shared_cognate_thq_no_homoplasie.csv"))


