rm(list=ls())
library(ape)
library(phangorn)
library(Matrix)
library(matrixStats)
library(castor)
library(here)
library(purrr)
library(tidyverse)
library(phytools)
library(TreeTools)
library(adephylo)

source(".Rprofile")

# probability that there exist one leaf l descendant of i such that there no mutation on <i,l>
B <- function(i, tree, M, pi0, pi1){
  pi1 * C1(i, tree, M) + pi0 * C0(i, tree, M)
}

# probability that there exist one leaf l descendant of i such that there no mutation on <i,l> condition to Y_i = 1
C1 <- function(i, tree, M){
  1 - C1_bar(i, tree, M)
}

# probability that for all path from i to desc(i) there is at least one mutation condition to Y_i = 1
C1_bar <- function(i, tree, M){
  if (i %in% 1:(tree$Nnode+1)){ 
    return(0)}
  else {
    childrens <- Descendants(tree , i , type = c("children"))
    return(
      (1 - exp(- M[2,1] * find_branch_length(tree, i, childrens[1])) + exp(- M[2,1] * find_branch_length(tree, i, childrens[1])) * C1_bar(childrens[1], tree, M)) *
        (1 - exp(- M[2,1] * find_branch_length(tree, i, childrens[2])) + exp(- M[2,1] * find_branch_length(tree, i, childrens[2])) * C1_bar(childrens[2], tree, M)) 
      )
  }
}

# probability that there exist one leaf l descendant of i such that there no mutation on <i,l> condition to Y_i = 0
C0 <- function(i, tree, M){
  1 - C0_bar(i, tree, M)
}

# probability that for all path from i to desc(i) there is at least one mutation condition to Y_i = 0
C0_bar <- function(i, tree, M){
  if (i %in% 1:(tree$Nnode+1)){ 
    return(0)}
  else {
    childrens <- Descendants(tree , i , type = c("children"))
    return(
      (1 - exp(- M[1,2] * find_branch_length(tree, i, childrens[1])) + exp(- M[1,2] * find_branch_length(tree, i, childrens[1])) * C1_bar(childrens[1], tree, M)) *
        (1 - exp(- M[1,2] * find_branch_length(tree, i, childrens[2])) + exp(- M[1,2] * find_branch_length(tree, i, childrens[2])) * C1_bar(childrens[2], tree, M)) 
    )
  }
}


# on every tree of a phylogeny, probability that there exist one leaf l descendant of i such that there no mutation on <i,l>
B_phylo <- function(phylo, pi0, pi1, clock_rate, mutation_rate, burnin = 0.9){
  
  # create stochastic matrix
  Q <- cbind(
    c(-pi0, pi0),
    c(pi1, -pi1)
  )
  # scale
  Q_scaled <- (clock_rate * mutation_rate) / (2 * pi0 * pi1) * Q
  
  # length of the phylogeny
  M <- length(phylo)
  
  # compute the probability
  mean(sapply(phylo[ceiling(burnin * M):M], function(tree) {
    root <- find_root(tree)
    B(root, tree, Q_scaled, pi0, pi1)
  }))
}


# import phylogenies
phylo_st_uni <- read.nexus( here('data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees'))
phylo_bantu <- read.nexus(here('data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.trees'))
phylo_bantu_subset <- read.nexus(here('data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.trees'))
phylo_bantu_subset2 <- read.nexus(here('data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.trees'))
#phylo_st_ht <- read.nexus(here('data/real/st_ctmc-strict-fbd-ht/st_ctmc-strict-fbd-heterogene.trees'))
#phylo_tea <- read.nexus(here('data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees'))
phylo_iecor <- read.nexus(here('/Users/kopp/Desktop/thesis/indo-euro/IECoR_Suppl_Files/02_Alternative_Models_M1_M2_M4/IECoR_M1_CTMC_Gamma_1_Rate_For_All_Mgs/IECoR_M1_CTMC_Gamma_1_Rate_For_All_Mgs_combined.trees'))

# import parameters
bounds_real_tb_by_sens <- read_csv(here("output/results/bounds_real_tb_by_sens.csv"))
foo_bantu <- filter(bounds_real_tb_by_sens, family == 'Bantu')
foo_bantu_subset <- filter(bounds_real_tb_by_sens, family == 'Bantu_subset')
foo_bantu_subset2 <- filter(bounds_real_tb_by_sens, family == 'Bantu_subset2')
foo_ie <- filter(bounds_real_tb_by_sens, family == 'IE')
foo_st <- filter(bounds_real_tb_by_sens, family == 'ST')
foo_st_by_sens <- filter(bounds_real_tb_by_sens, family == 'ST_by_sens')
foo_tea <- filter(bounds_real_tb_by_sens, family == 'TEA')

# probabilities on phylogenies 
B_phylo(phylo_st_uni, foo_st$pi0, foo_st$pi1, foo_st$clock_rate, foo_st$mu)
B_phylo(phylo_bantu, foo_bantu$pi0, foo_bantu$pi1, foo_bantu$clock_rate, foo_bantu$mu, burnin = 0.99)
B_phylo(phylo_bantu_subset, foo_bantu_subset$pi0, foo_bantu_subset$pi1, foo_bantu_subset$clock_rate, foo_bantu_subset$mu)
B_phylo(phylo_bantu_subset2, foo_bantu_subset2$pi0, foo_bantu_subset2$pi1, foo_bantu_subset2$clock_rate, foo_bantu_subset2$mu)
B_phylo(phylo_iecor, foo_ie$pi0, foo_ie$pi1, foo_ie$clock_rate, foo_ie$mu, burnin = 0.99)

# on a tree multiply the branch length by a coeficiant 
B_scaling_tree <- function(coef, tree, pi0, pi1, clock_rate, mutation_rate){
  # scale tree
  new_tree <- tree
  new_tree$edge.length <- tree$edge.length * coef
  root <- find_root(tree)
  t_R <- max(node.depth.edgelength(new_tree))
  
  # parameters
  # create stochastic matrix
  Q <- cbind(
    c(-pi0, pi0),
    c(pi1, -pi1)
  )
  # scale
  Q_scaled <- (clock_rate * mutation_rate) / (2 * pi0 * pi1) * Q
  
  B(root, new_tree, Q_scaled, pi0, pi1)
}




# sino-tibetain 
B_prob_st <- sapply(
  seq(1, 30, length.out = 30)/max(node.depth.edgelength(phylo_st_uni[[150]])), 
  function(coef) 
    B_scaling_tree(coef, phylo_st_uni[[150]], foo_st$pi0, foo_st$pi1, foo_st$clock_rate, foo_st$mu)
  )

# indo-european
B_prob_iecor <- sapply(
  seq(1, 30, length.out = 30)/max(node.depth.edgelength(phylo_iecor[[150]])), 
  function(coef) 
    B_scaling_tree(coef, phylo_iecor[[150]], foo_ie$pi0, foo_ie$pi1, foo_ie$clock_rate, foo_ie$mu)
)

# bantu
B_prob_bantu <- sapply(
  seq(1, 30, length.out = 30)/max(node.depth.edgelength(phylo_bantu[[150]])), 
  function(coef) 
    B_scaling_tree(coef, phylo_bantu[[150]], foo_bantu$pi0, foo_bantu$pi1, foo_bantu$clock_rate, foo_bantu$mu)
)


delta_prob <- 
  bind_rows(
    tibble(family = "ST", delta_prob = B_prob_st, age = 1:30),
    tibble(family = "IE", delta_prob = B_prob_iecor, age = 1:30),
    tibble(family = "Bantu", delta_prob = B_prob_bantu, age = 1:30)) |> 
  mutate(delta_prob = round(delta_prob, 3)) |>
  #pivot_wider(names_from = family, values_from = delta_prob) |> 
  pivot_wider(names_from = age, values_from = delta_prob) |> 
  as_tibble()
  #relocate(age)
  
write_csv(delta_prob, here('output/results/delta_prob.csv'))










