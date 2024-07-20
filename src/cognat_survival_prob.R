library(here)
library(tidyverse)
library(ape)
library(adephylo)
library(dplyr)
library(purrr)
library(castor)
library(phytools)
library(maps)
library(phangorn)


# Prob of a trait born in i to survive in bot subtree
Q <- function(tree, i, mu) {
  childrens <- Descendants(tree, i, type = c("children"))
  delta(tree, i, childrens[1], mu) * P(tree, childrens[1], mu) *
    delta(tree, i, childrens[2], mu) * P(tree, childrens[2], mu)
}


# This function output the branch length in a tree between two nodes
find_branch_length <- function(tree, father, children) {
  indice <- which(tree$edge[, 1] == father & tree$edge[, 2] == children)
  tree$edge.length[indice]
}

# This function outputs the probability of a cognat to survive between node i to node j 
delta <- function(tree, i, j, mu) {
  exp(-mu * find_branch_length(tree, i, j))
}

# See notations in article Gray&Nicholls 2018
# Sequence u
u <- function(tree, i, mu) {
  if (i %in% 1:(tree$Nnode + 1)) {0} 
  else {
    childrens <- Descendants(tree, i, type = c("children"))
    (1 - delta(tree, i, childrens[1], mu) + delta(tree, i, childrens[1], mu) * u(tree, childrens[1], mu)) *
      (1 - delta(tree, i, childrens[2], mu) + delta(tree, i, childrens[2], mu) * u(tree, childrens[2], mu))
    }}

# Probability of a trait born in i to survive in 1 or more leave
P <- function(tree, i, mu) {
  1 - u(tree, i, mu)
}


# Separate the tip of tree into two set regarding to the subtrees induced by the root
tips_of_subtree = function(tree){
  root = find_root(tree)
  child1 = Descendants(tree, root, type = c("children"))[[1]]
  
  tipsA = tree$tip.label[Descendants(tree, child1, type = c("tips"))[[1]]]
  tipsB = setdiff(tree$tip.label,tipsA)
  
  return(list(A=tipsA, B=tipsB))
} 


# 1. Theoretical value
# Initialize the number of meanings 
bounds_real_tb = read_csv(here("output/results/bounds_real_tb.csv")) |> 
  mutate(n_meanings = c(254,150,150,150,170,200,200)) |> 
  relocate(n_meanings, .after = n_cogsets)

# Load the trees 
tree_bantu = read.tree(here("output/results/bantu/bantu_ctmc-strict-bd_tree.nex"))
tree_bantu_subset = read.tree(here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tree.nex"))
tree_bantu_subset2 = read.tree(here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tree.nex"))
tree_ie = read.tree(here("output/results/ie/iecor_ctmc-strict-M1_tree.nex"))
tree_st = read.tree(here("output/results/st/st_ctmc-strict-fbd_tree.nex"))
tree_tea = read.tree(here("output/results/tea/tea_ctmc-strict-fbd-constrained_tree.nex"))
tree_st_bysens = read.tree(here("output/results/st_by_sens/st_ctmc-strict-fbd_by_sens_tree.nex"))

# Transition from 1 to 0 i.e death trait rate
transition = bounds_real_tb |> 
  mutate(pi10= pi1/(pi0^2+pi1^2)) |> 
  select(family,pi10)


# Functions to transform into tidy 
compute_survival_prob_by_ages <- function(tree, n_sens, mu, age_min = 1, age_max = 20) {
  t <- as.numeric(distRoot(tree, 1))
  root <- find_root(tree)
  ks <- (age_min:age_max) / t
  
  # Create 20 trees with modified edge lengths
  trees_scaled <- purrr::map(ks, function(x) {
    tree_k <- tree
    tree_k$edge.length <- as.numeric(tree$edge.length) * x
    return(tree_k)
  })
  
  # Compute theoretical survival probabilities
  tibble(
    age = age_min:age_max,
    q_theo = map2_dbl(trees_scaled, ks, ~ Q(.x, root, mu * .y) * n_sens)
  )
}


compute_survival_prob_by_ages2 = function(tree,n_sens,mu){
  t = as.numeric(distRoot(tree,1))
  root = find_root(tree)
  ks = (1:20)/t
  
  # create 20 trees
  for (k in ks) {
    tree_k = tree
    tree_k$edge.length = tree_k$edge.length * k 
    assign(paste0("tree_", k*t), tree_k)
  }
  
  q_theo =c()
  for (k in 1:20){
    tree_k = get(paste0("tree_",k))
    q_theo=c(q_theo,Q(tree,root,mu*k/t)*n_sens)
  }
  
  
  df_theo = data.frame(
    age = 1:length(q_theo),
    q_theo = q_theo
  )
  
  return(df_theo)}

# Apply the function to each family and combine the results
trees <- c(tree_tea, tree_bantu, tree_bantu_subset, tree_bantu_subset2, tree_ie, tree_st, tree_st_bysens)

qs_tb <- map_df(1:length(trees), function(i) {
  family <- bounds_real_tb$family
  n_meanings <- bounds_real_tb$n_meanings
  pi10 <- transition$pi10
  
  compute_survival_prob_by_ages2(trees[[i]], n_meanings[i], pi10[i]) |>
    mutate(family = family[i])
}) |>
  pivot_wider(names_from = age, values_from = q_theo) |>
  group_by(family)  

# Write outputs
write_csv(qs_tb, here("output/results/qs_tb.csv"))







qs_tb_min <- qs_tb |>
  pivot_longer(cols = -family, names_to = "millennium", values_to = "value") |>
  mutate(millennium = as.numeric(millennium)) %>%
  left_join(qs_tb %>% select(family, value_1 = '1' ), by = "family") |>
  filter(value < 0.05*value_1) |>
  group_by(family) |>
  summarise(min_millennium = min(millennium), .groups = 'drop')
write_csv(qs_tb_min, here("output/results/qs_tb_min.csv"))


# 2. Number of common meanings between a leaf and a groups  
# Separate the tip of tree into two set regarding to the subtrees induced by the root
tip_set = tips_of_subtree(tree_st)
tip_set$A
tip_set$B


# Convert nexus data into dataframe
trait = read.nexus.data(here("data/real/st_ctmc-strict-fbd/st.nex")) 
languages = names(trait)

trait_by_family = trait|> 
  as_tibble(.name_repair = 'unique') |>
  mutate_all(~replace(as.numeric(replace(., . == "?", NA)), is.na(.), NA)) |>
  t() |> 
  as_tibble() |> 
  mutate(Family = languages) |>
  mutate(
    Subtree = case_when(
      Family %in% tip_set$A ~ 1,
      Family %in% tip_set$B ~ 2,
      TRUE ~ NA)) |> 
  relocate(c(Family,Subtree), .before = 1) 
  

sub1 = trait_by_family |>
  filter(Subtree==1) 

sub2 = trait_by_family |>
  filter(Subtree==2)


# Count number of common trait of two languages
count_common_traits = function(trait1, trait2){
  sum(trait1==1&trait2==1, na.rm = T)
}

# Count number of common trait : 1 language vs 1 df 
count_common_traits_df = function(df, trait2){
  res = c()
  for (i in 1:nrow(df)){
    #print(i)
    row = df[i,]
    #print(trait2)
    res = rbind(res, count_common_traits(row,trait2))
  }
  res
}


count_common_traits_df2 = function(df1,df2){
  # Res matrix
  res = matrix(ncol = nrow(df2) , nrow= nrow(df1))
  colnames(res) = df2$Family
  rownames(res) = df1$Family
  
  # Just the trait 
  df1 = df1 |> select(-Family, -Subtree)
  df2 = df2 |> select(-Family, -Subtree)
  
  for (j in 1:nrow(df2)){
    row = df2[j,]
    res[,j] = count_common_traits_df(df1,row)
  }
  res
}


res_final = count_common_traits_df2(
  df1 = sub1 ,
  df2 = sub2 
  )





