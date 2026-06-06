# ------------------------------------------------------------------------------
# Script Name: 03_simulation_analyses_main.R
# Description: This script consolidates and tidies all the data required for 
#              plotting and analysis of phylogenetic reconstruction accuracy. 
#              It performs the following tasks:
#                - Loads and merges summary statistics from previous simulations
#                - Computes aggregated node reconstruction counts and proportions
#                - Saves cleaned datasets in .csv and .rds formats for plotting
#                - Computes summary statistics on marginal posterior probabilities
#                - Prepares datasets for modeling and runs logistic regressions
#              Data covers comparisons between true trees and both consensus 
#              and MCC summary trees across varying simulation ages.
# ------------------------------------------------------------------------------
library(here)
library(broom)
library(corrplot)
library(HDInterval)
library(tidyverse)
N_sim <- 50

# load data --------------------------------------------------------------------
# marginal probability of the first split with IC 
marginal_probability_first_split_ic <- read.csv(here("output/results/marginal_prob_first_split_ic.csv"))

hdi_results <- map(1:17, ~ {
  test <- marginal_probability_first_split_ic |>
    filter(tree_age == .x) |>
    select(prob_mean)
  list(
    lower = hdi(test, credMass = 0.90)[1],
    upper = hdi(test, credMass = 0.90)[2]
  )
})
hdi_lower <- map_dbl(hdi_results, "lower")
hdi_upper <- map_dbl(hdi_results, "upper")

marginal_probability_first_split_ic = marginal_probability_first_split_ic |>
  group_by(tree_age) |>
  summarise(
    prob_mean = mean(prob_mean),
    prob_inf = mean(prob_inf),
    prob_sup = mean(prob_sup)
  )
marginal_probability_first_split_ic$prob_inf <- hdi_lower
marginal_probability_first_split_ic$prob_sup <- hdi_upper
marginal_probability_first_split_ic

write.csv(marginal_probability_first_split_ic, here("output/results/marginal_prob_first_split_ic.csv"))

# frequency of good reconstruction of all the nodes in of the mcc
mcc_to_true_TF <- read.csv(here("output/results/resume_to_true_TF.csv")) |> 
  filter(type == 'mcc')

# for each mcc tree, age between the root and the first split of the true tree
first_split_age_mcc <- read.csv(here("output/results/first_split_age_mcc.csv"))
first_split_age_cs <- read.csv(here("output/results/first_split_age_cs.csv"))

# frequency of good reconstruction of all the nodes in consensus tree (true -> summary)
# the value of node represent the node in the true tree
true_false_uncertain <- read.csv(
  file = here("output/results/true_false_uncertain_nodes.csv"),
  sep = ",",
  header = T) |> 
  rename(age = tree_age, simulation = tree_simulation_number) |>
  # add 0/1/2 for false/true/uncertain nodes for the consensus tree
  mutate(
    value = case_when(
      T_F_U == FALSE ~ 0,
      T_F_U == TRUE & state == "rateau" ~ 2,
      T_F_U == TRUE & state == "regular" ~ 1
    )
  )

# number of node in the consensus and mcc tree
df_number_of_nodes <- read.csv(
  file = here("output/results/number_nodes_mcc_cs.csv"),
  sep = ",",
  header = T
)

# marginal probability of the first split in the mcc and consensus tree
prob_first_split_summary = read.csv(
  here("output/results/marginal_prob_first_split_mcc_consensus.csv")
)

# process data -----------------------------------------------------------------

# marginal probability of the first split in the mcc 
prob_first_split_mcc = prob_first_split_summary |> 
  select(- cs_prob, - node_cs, - node_mcc) |> 
  group_by(age) |> 
  summarise(mean_mcc_prob = mean(mcc_prob, na.rm=T), .groups='drop') |> 
  ungroup() |> 
  write.csv(file = here("output/results/prob_first_split_mcc.csv"), row.names = FALSE)

# number of node in the summary tree
df_number_of_nodes_avg <- df_number_of_nodes |>
  group_by(age) |>
  summarise(
    n_mcc = mean(n_mcc, na.rm = TRUE),
    n_consensus = mean(n_consensus, na.rm = TRUE)
  ) |>
  write.csv(file = here("output/results/number_of_nodes_summary.csv"), row.names = FALSE)


# for the consensus trees, count the number of true, false and uncertain nodes
# with special labels for the plot
count_true_to_cs <- true_false_uncertain |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = sum(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

write.csv(count_true_to_cs, here("output/results/count_true_to_cs.csv"))

# for the mcc tree count the number of true, false

# for each summary tree, age, simulation this dataframe indicates the proprtions 
# of true and false nodes
resume_to_true_grouped <- read_csv(here("output/results/resume_to_true_TF.csv"), col_names = T)|> 
  rename(N_node = exist, exist = N_nodes) |>
  mutate(exist = as.numeric(exist)) |>  # TRUE -> 1, FALSE -> 0
  group_by(type, age, simulation, exist, N_node) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(type, age, simulation) |>
  mutate(proportion = n / N_node) |>
  ungroup() |>
  mutate(exist = factor(exist, levels = c("0", "1")))

# number of well reconstructed node in the mcc tree
count_true_to_mcc <- resume_to_true_grouped |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(n_mean = sum(n)/N_sim, .groups = "drop") |>
  tidyr::pivot_wider(names_from = exist, values_from = n_mean, names_prefix = "exist_") |>
  mutate(
    total = exist_0 + exist_1
  ) |>
  pivot_longer(cols = starts_with("exist_"), names_prefix = "exist_", names_to = "exist", values_to = "n_mean") |>
  mutate(
    exist = factor(exist, levels = c("0", "1")),
    y_label = ifelse(exist == "1", 0, total) # alignement manuel
  ) 

write.csv(count_true_to_mcc,here("output/results/count_true_to_mcc.csv"))


# count the numer of true node from de consensus tree to the true tree
count_cs_to_true <- resume_to_true_grouped |>
  filter(type == "consensus") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n), .groups = "drop") |>
  tidyr::pivot_wider(names_from = exist, values_from = n_mean, names_prefix = "exist_") |>
  mutate(
    total = exist_0 + exist_1
  ) |>
  pivot_longer(cols = starts_with("exist_"), names_prefix = "exist_", names_to = "exist", values_to = "n_mean") |>
  mutate(
    exist = factor(exist, levels = c("0", "1")),
    y_label = ifelse(exist == "1", 0, total) 
  )

write_csv(count_cs_to_true, here("output/results/count_cs_to_true.csv"))

# proportion of true,false and uncertain nodes in the consensus tree (true -> consensus)
prop_true_to_cs <- true_false_uncertain |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = sum(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

write_csv(prop_true_to_cs, here("output/results/prop_true_to_cs.csv"))

# proportion of true, false node from the mcc to the true tree
prop_mcc_to_true <- resume_to_true_grouped |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(n_mean = sum(n)/N_sim, .groups = "drop") 

write_csv(prop_mcc_to_true, here("output/results/prop_mcc_to_true.csv"))

# proportion of true, false node from the consensus to the true tree
prop_cs_to_true <- resume_to_true_grouped |>
  filter(type == "consensus") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n), .groups = "drop") |> 
  ungroup()

write_csv(prop_cs_to_true, here("output/results/prop_cs_to_true.csv"))

# number of true nodes in the summary tree
prop_mcc_to_true <- resume_to_true_grouped |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n), .groups = "drop") |> 
  ungroup()

write_csv(prop_mcc_to_true, here("output/results/prop_mcc_to_true.csv"))

# regression on the first split ------------------------------------------------
# regression mcc
df_reg_mcc <- mcc_to_true_TF |>
  inner_join(prob_first_split_summary, by = c("age", "simulation")) |>
  filter(node == node_mcc) |>
  full_join(first_split_age_mcc, by = c("age","simulation")) |>
  select(-type, -exist, -node_cs,-node_mcc, -cs_prob, -X)|>
  rename(y = N_nodes) |> 
  mutate(root_split_age_prop = as.numeric(root_split_age)/root_age) |>
  rename(first_split_prob = mcc_prob) |> 
  mutate(y = as.factor(y))

# regression with age, probability of the first split and age of the first split
# model mcc regression 
model_mcc <- glm(y ~ age + first_split_prob + root_split_age_prop, data = df_reg_mcc, family = 'binomial')
summary(model_mcc)

# regression cs 
resume_to_true_TF_cs <- read.csv(here("output/results/resume_to_true_TF.csv")) |> 
  filter(type == 'consensus')

df_reg_cs <- resume_to_true_TF_cs |>
  inner_join(prob_first_split_summary, by = c("age", "simulation")) |>
  filter(node == node_cs ) |>
  full_join(first_split_age_cs, by = c("age","simulation")) |>
  select(-type, -exist, -node_cs,-node_mcc, -mcc_prob, -X)|>
  rename(y = N_nodes) |>
  # delete the tree for which the first split is a leaf
  filter(!is.na(y)) |> 
  mutate(root_split_age_prop = as.numeric(root_split_age)/root_age) |>
  rename(first_split_prob = cs_prob) |> 
  mutate(y = as.factor(y))

# model consensus regression
model_cs <- glm(y ~ age + first_split_prob + root_split_age_prop, data = df_reg_cs, family = 'binomial')
summary(model_cs)

# regression with age, probability of the first split
model_mcc2 <- glm(y ~ age + first_split_prob, data = df_reg_mcc, family = 'binomial')
model_cs2 <- glm(y ~ age + first_split_prob, data = df_reg_cs, family = 'binomial')

summary(model_mcc2)
summary(model_cs2)


# Sauvegarde des modèles dans des fichiers .rds
saveRDS(model_mcc2, here("output/results/model_mcc.rds"))
saveRDS(model_cs2, here("output/results/model_cs.rds"))

# correlation matrix between covariates
cor_matrix <- cor(df_reg_mcc[c("age", "first_split_prob", "root_split_age_prop")], use = "complete.obs")
#corrplot(cor_matrix, method = "circle", type = "full")



