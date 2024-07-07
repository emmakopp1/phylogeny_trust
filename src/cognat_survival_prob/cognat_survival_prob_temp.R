library(here)
library(tidyverse)
library(ape)
library(adephy)
library(dplyr)

# 1. Theoretical value
# Initialise the number of meanings 
bounds_real_tb = read_csv(here("output/results/bounds_real_tb.csv")) |> 
  mutate(n_meanings = c(254,150,150,150,170,200)) |> 
  relocate(n_meanings, .after = n_cogsets)

# Load the tree 
tree_bantu = read.tree(here("output/results/bantu/bantu_ctmc-strict-bd_tree.nex"))
tree_bantu_subset = read.tree(here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tree.nex"))
tree_bantu_subset2 = read.tree(here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tree.nex"))
tree_ie = read.tree(here("output/results/ie/iecor_ctmc-strict-M1_tree.nex"))
tree_st = read.tree(here("output/results/st/st_ctmc-strict-fbd_tree.nex"))
tree_tea = read.tree(here("output/results/tea/tea_ctmc-strict-fbd-constrained_tree.nex"))

# Transition from 1 to 0 i.e death trait rate
transition = bounds_real_tb |> 
  mutate(pi10= pi1/(pi0^2+pi1^2)) |> 
  select(family,pi10)


compute_survival_prob <- function(family, tree, n_meanings, pi10) {
  return(compute_survival_prob_by_ages(tree, n_meanings, pi10) |>
    mutate(family = family))
}

# Apply the function to each family and combine the results
family <- bounds_real_tb$family
tree <- c(tree_tea, tree_bantu,tree_bantu_subset,tree_bantu_subset2,tree_ie,tree_st)
n_meanings <- bounds_real_tb$n_meanings
pi10 <- transition$pi10

qs_tb <- map_df(1:length(phylo), function(i) {compute_survival_prob(family[i], tree[[i]], n_meanings[i], pi10[i])}) |>
  pivot_wider(names_from = age, values_from = q_theo) |>
  group_by(family) 
  
# Write outputs
write_csv(qs_tb, here("output/results/qs_tb.csv"))



# 2. Empirical value pas clair pour le moment - potentiel sert à rien

# Here I compute the expected number of meaning surviving from the root two subgroups of the root
# and the observed number of meaning which survived in the two subgroups 

l_theo_empirical_st = t(as.data.frame(
  unname(sapply(tree_st, 
                function(arbre) 
                  c(
                    distRoot(arbre,1),
                    as.numeric(Q(arbre,find_root(arbre))*200),
                    survival_frequency(data_st$path_cognates,arbre)
                  ))),
  row.names = c("age","theorique","empirical")
))


# 3. Number of common meanings between a leaf and a groups  
# Separate the tip of tree into two set regarding to the subtrees induced by the root
tip_set = tips_of_subtree(tree_st)
tip_set$A
tip_set$B


# Convert nexus data into dataframe
df_cognates = as_tibble(read.nexus.data(here("data/real/st_ctmc-strict-fbd/st.nex"))) |> 
  t() |>
  rename(Language = V1)

#df_cognates = data_to_df(data_st$path_cognates)

# Dataframe of Sinitc group & other groups

v1 <- df_cognates %>% filter(tip_set$B) 

mutate_all(~replace(as.numeric(replace(., . == "?", NA)), is.na(.), NA))

v1 = df_cognates[tip_set$B,] |> as.numeric()
v2 = df_cognates[tip_set$A,]


resultats = c()

for (sinitic_row in tip_set$B){
  #x = apply(v2, 1, function(row) sum(v1[sinitic_row,] == row & v1[sinitic_row,] == 1, na.rm = TRUE))
  x = sum(apply(v2, 2, max,na.rm=T)&v1[sinitic_row,],na.rm=T)
  resultats = cbind(resultats,x)
}

colnames(resultats) = tip_set$B
resultats



