rm(list=ls())
library(ape)
library(phangorn)
library(Matrix)
library(matrixStats)
library(castor)
library(here)
library(purrr)
library(tidyverse)

source(".Rprofile")
# paramter ---------------------------------------------------------------------
K_tot <-3963

# binary covarion --------------------------------------------------------------
#alpha <- 6.371e-2 # slow mode
#s <- 0.201 # fast mode
#pi0 <- 0.995
#pi1 <- 1- pi0
#c <- 1.928E-3
#mu <- 1
#lambda <- (c * mu)/(2 * pi0 * pi1)


#Q <- cbind(
#  c(-pi0 - s*pi0, pi0, s*pi0, 0),
#  c(pi1, -pi1 - s*pi1, 0, s* pi1),
#  c(s*pi0, 0, -s*pi0 - alpha*pi0, alpha*pi0),
#  c(0, s*pi1, alpha*pi1, -s*pi1 - alpha*pi1)
#) 

path_bi_cov <- here("data/real/st_covarion-strict-fbd/st_covarion-strict-fbd.trees")
phylo_st_bcov <- read.nexus(path_bi_cov)
# ctmc -------------------------------------------------------------------------
pi0 <- 0.943
pi1 <- 5.695e-2
mu <- 1 
c <- 1.807e-2
lambda <- c * mu 
q <-  (pi0 + pi1)/(2 * pi0 * pi1)

Q <- cbind(
  c(-pi0, pi0),
  c(pi1, -pi1)
) 

path_ctmc_ht <- '/Users/kopp/Documents/phylogeny_trust/data/real/st_covarion-strict-fbd-ht/st_covarion-strict-fbd-heterogene.trees'
path_ctmc_uni <- here('data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees')
phylo_ctmc <- read.nexus(path_ctmc_uni)

# ∆S ---------------------------------------------------------------------------
# ∆S for binary covarion model one tree, for one trait
ub_DS_bicov_tree <- function(tree, lambda, Q, s_node){
  n_row <- dim(Q)[1]
  somme <- 0
  for (l in 1:length(tree$tip.label)){
    #print(l)
    ancestors <- nodepath(tree, l, s_node)
    edges_root_to_tip <- head(tree$edge.length[match(ancestors,tree$edge[,2])],-1)
    edges_root_to_tip <- edges_root_to_tip[!is.na(edges_root_to_tip)]
    
    prod <- 1
    for (edge in edges_root_to_tip){
      prod <- prod * (
        1 - sum(colMins(
        as.matrix(expm(lambda * edge * Q) +  diag(Inf,n_row))))
        )}
    somme <- somme + prod
  }
  return(somme)
}

#library(ape)
#tree <-read.tree('/Users/kopp/Desktop/test.tree')

# Test 2 -----------------------------------------------------------------------
# calcul de la borne pour un arbre et plusieurs cuttings points (les 20 noeuds les plus profonds)
tree <- phylo_ctmc[[999]]
s_set <- getNodesByDepth(tree)
K <- sum(parameters_st$n_cogsets)
sapply(s_set[1:1], function(s) K*ub_DS_bicov_tree(tree, lambda, Q, s))


# MARCHE PLUS
# ∆S for binary covarion model phylogeny
ub_DS_bicov_phylo <- function(phylo, lambda, Q, burnin){
  M <- length(phylo)
  mean(sapply(phylo[ceiling(burnin*M):M], function(tree) ub_DS_bicov_tree(tree, lambda, Q)))
}


# pour un modèle binary covarion, la seule borne qui existe est moins fine, et donc
# peu utile en pratique. pk jai dis ça ?


# ∆S for ctmc model for one tree
ub_DS_ctmc_tree <- function(tree, lambda, q){
  
  root <- find_root(tree)
  somme <- 0
  
  for (l in 1:length(tree$tip.label)){
    ancestors <- nodepath(tree, l, root)
    edges_root_to_tip <- head(tree$edge.length[match(ancestors,tree$edge[,2])],-1)
    t <- sum(edges_root_to_tip) 
    
    cat('t=',t,'epx=',exp(-lambda * q * t),'\n')
    somme <- somme + exp(-lambda * q * t)
  }
  return(somme)
}

ub_DS_ctmc_tree(tree,lambda,q*mu)

#  ∆S for ctmc model for a phylogeny
ub_DS_ctmc_phylo <- function(phylo, lambda, q, burnin){
  M <- length(phylo)
  mean(sapply(phylo[ceiling(burnin*M):M], function(tree) ub_DS_ctmc_tree(tree, lambda, q)))
}

# ∆T ---------------------------------------------------------------------------

# add cutting point
calibrations <- read_csv(here("output/results/calibration.csv")) |>
  select(family, calibration, tip, s)

# tipages summary
tipages_summary <- read_csv(here("output/results/tipages_summary.csv")) |>
  full_join(calibrations, by = c("family", "tip")) |>
  mutate(s = if_else(is.na(s), depth, s)) |>
  group_by(family, tip) |>
  filter(s == min(s)) |>
  mutate(s = ifelse(depth > 0, s + 5, s))|>
  ungroup()


# ∆T for ctmc model for on leaf, cutting points and parameters
ub_DT_ctmc_leaf <- function(l, tree, lambda, q, s){
  # t(l)
  root <- find_root(tree)
  ancestors <- nodepath(tree, l, root)
  edges_root_to_tip <- head(tree$edge.length[match(ancestors,tree$edge[,2])],-1)
  # age of the root (t_R) - t(l)
  d <- max(node.depth.edgelength(tree)) - sum(edges_root_to_tip)
  #cat(tree$tip.label[l], round(s-d,3), "\n")
  return(exp(-lambda * q * (s - d)))
  
}

# ∆T for ctmc model for one tree, cutting points, parameters
ub_DT_ctmc_tree <- function(tree, lambda, q, cutting_points, K){
  
  root <- find_root(tree)
  somme <- 0
  
  
  for (l in 1:length(tree$tip.label)){
    s <- cutting_points[l]
    somme <- somme + ub_DT_ctmc_leaf(l, tree, lambda, q, s)
  }
  return(K * somme)
}

# test1 ------------------------------------------------------------------------
#s_set
s_node <- s_set[1]
s <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[s_node]
ss <- rep(s, length(phylo_ctmc[[1000]]$tip.label))
ub_DT_ctmc_tree(phylo_ctmc[[1000]], lambda, q, ss, sum(parameters_st$n_cogsets))

# ∆T for ctmc model, one tree, parameters
ub_DT_ctmc_auto_tree <- function(tree, lambda, q, K, coef, s_option = 'constant'){
  N <- tree$Nnode
  internal_nodes <- (N + 2):(2 * N + 1)
  internal_nodes_depth <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[internal_nodes]
  internal_nodes_depth <- as.numeric(internal_nodes_depth) * as.numeric(coef)

  # compute ∆ for s constant taking internal nodes depth values
  map_df(internal_nodes_depth, function(node_depth) {
    cutting_points <- rep(node_depth, length(tree$tip.label))
    tibble(s= node_depth, ub_DT = ub_DT_ctmc_tree(tree, lambda, q, cutting_points, K))
  })
}

#ub_DT_ctmc_auto_tree(phylo_ctmc[[150,]])

#∆T for ctmc model for one phylogeny cutting points and parametrs
ub_DT_ctmc_phylo <- function(phylo, lambda, q, cutting_points, K, burnin){
  M <- length(phylo)
  mean(sapply(phylo[ceiling(burnin*M):M], function(tree) ub_DT_ctmc_tree(tree, lambda, q, cutting_points, K)))
}

# uniform tree prior -----------------------------------------------------------
# uniform prior for a tree given a number of leaf 
uniform_tree_prior_s <- function(S){
  1/(factorial(2*S - 3) / (2^(S - 2) * factorial(S -2 )))
}


# uniform tree prior for a tree given a tree and a cutting point
uniform_tree_prior <- function(tree, s){
  # distance to the tips
  nodes_depth <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)
  nodes_depth <-nodes_depth[-(1:length(tree$tip.label))] #internal nodes
  S <- sum(nodes_depth >= s) + 1 # at the right of the internal node s
  uniform_tree_prior_s(S)
}

# uniform tree prior given all the possible cutting points (internal nodes)
uniform_tree_prior_auto <- function(tree){
  
  N <- tree$Nnode
  internal_nodes <- (N + 2):(2 * N + 1)
  internal_nodes_depth <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[internal_nodes]
  
  # compute ∆ for s constant taking internal nodes depth values
  map_dbl(internal_nodes_depth, function(s) {
    uniform_tree_prior(tree,s)
  })
}


# ∆T < prior + model_contribution ----------------------------------------------
ub_DT_prior_model <- function(tree, lambda, q, K, coef){
  N <- tree$Nnode
  internal_nodes <- (N + 2):(2 * N + 1)
  internal_nodes_depth <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[internal_nodes]
  internal_nodes_depth <- coef * internal_nodes_depth
  
  ub_DT_node_depth <- tibble(
    node = internal_nodes,  
    s = internal_nodes_depth,
    model_contrib = ub_DT_ctmc_auto_tree(tree, lambda, q, K, coef)$ub_DT,
    uniform_prior = uniform_tree_prior_auto(tree)
  ) |> 
    mutate(ub_DT = coalesce(uniform_prior, 0) + model_contrib) |>
    mutate(uniform_prior = round(uniform_prior,3)) |>
    arrange(ub_DT) 
  
  return(ub_DT_node_depth)
}

#∆T minimizer over cutting points
ub_DT_min_tree <- function(tree, lambda, q, K, coef, s_option = 'constant'){
  ub_DT_prior_model(tree, lambda, q, K, coef) |> 
    filter(ub_DT == min(ub_DT)) 
}


# heterogene parameter, process ------------------------------------------------

# tracelog summary
tracelog_summary <- read_csv(here("output/results/tracelog_summary.csv")) |>
  filter(!(family == "ST_by_sens" & concept == "the_name")) |>
  filter(!(family == "ST_by_sens" & concept == "four")) 

parameters_st <- tracelog_summary |>
  filter(family=='ST_by_sens') |>
  select(pi0, pi1,mu, clock_rate, t_R, concept, n_cogsets) |> 
  relocate(concept,.before = pi0)

pi0 <- parameters_st$pi0[1]
pi1 <- parameters_st$pi1[1]
mu <- parameters_st$mu[1]
c <- parameters_st$clock_rate[1]
K <- parameters_st$n_cogsets[1]
lambda <- c * mu 
q <-  (pi0 + pi1)/(2 * pi0 * pi1)

## heteorgene parameter, test of ub_DT ----------------------------------------
M <- length(phylo_ctmc)
#results <- tibble()
M_seq <- seq(1,M,length=10)


coef <- 12

for (n_tree in M_seq){
  for (i in 1:nrow(parameters_st)) {
    
    # initialize parameters
    pi0 <- parameters_st$pi0[i]
    pi1 <- parameters_st$pi1[i]
    mu <- parameters_st$mu[i]
    c <- parameters_st$clock_rate[i]
    K <- parameters_st$n_cogsets[i]
    concept <- parameters_st$concept[i]
    lambda <- (c * mu)/(2 * pi0 * pi1) 
    q <-  (pi0 + pi1)
    
    # compute the new row
    new_row <- ub_DT_min_tree(phylo_ctmc[[n_tree]], lambda, q, K, coef) |> 
      mutate(model_contrib = round(model_contrib,3)) |>
      mutate(concept = concept,
             pi0 = pi0, 
             lambda = lambda, 
             coef = coef,
             n_tree = n_tree) |> 
      relocate(concept, .before = node)
    
    # bind rows
    results <- bind_rows(results, new_row)
  }
}  

# Vérifier le résultat final
glimpse(results)

res
sum(res)

# results ----------------------------------------------------------------------
coef_set = 1:11
res <- c()
for (coef in coef_set){
  res <- cbind(res, ub_DT_min_tree(phylo_ctmc[[150]], lambda, q, K, coef)$ub_DT)
}



# cutting points non constant
cutting_points =  tipages_summary |> filter(family == 'ST') |> select(s) |> as.numeric()

# ∆ with cutting point as a variable 
ub_DT_ctmc_tree(
  phylo_ctmc[[2]], 
  lambda, 
  q, 
  cutting_points, 
  K)

# Tests
#s = tipages_summary |> filter(family == 'ST') |> 
#  select(s) |> 
#  unlist() |> 
#  as.numeric()

#ub_DT_ctmc_leaf(1, phylo_ctmc[[2]], lambda, q, s)
#ub_DT_ctmc_tree(tree = phylo_ctmc[[150]], lambda = lambda, q = q, cutting_points = s, K = K)
#ub_DT_ctmc_phylo(phylo_ctmc, lambda, q, s, K, 0.98)
#ub_DT_ctmc_auto_tree(phylo_ctmc[[150]], lambda, q, K, coef)  
#ub_DT_cutting_points <- ub_DT(phylo_ctmc[[150]], lambda, q, K, coef)    
#ub_DT_ctmc_phylo(phylo_ctmc, lambda, q, cutting_points, K, 0.96)
#ub_DT_ctmc_phylo(phylo_ctmc, lambda, q, cutting_points, K, 0.98)
#ub_DS_ctmc_phylo(phylo_ctmc, lambda, q, 0.9) 
#ub_DS_bicov_phylo(phylo_ctmc, lambda, Q, 0.99)









