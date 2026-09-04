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
library(purrr)
library(tidyverse)
N_sim <- 50

# load data --------------------------------------------------------------------
# marginal probability of the first split with IC 
marginal_probability_first_split_ic <- read.csv(here("output/results/marginal_prob_first_split_ic.csv"))

hdi_results <- purrr::map(1:17, ~ {
  test <- marginal_probability_first_split_ic |>
    filter(tree_age == .x) |>
    select(prob_mean)
  list(
    lower = HDInterval::hdi(test, credMass = 0.90)[1],
    upper = HDInterval::hdi(test, credMass = 0.90)[2]
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

write.csv(marginal_probability_first_split_ic, here("output/results/marginal_prob_first_split_hdi.csv"))

# for each mcc tree, age between the root and the first split of the true tree
first_split_age_mcc <- read.csv(here("output/results/first_split_age_mcc.csv"))
first_split_age_cs <- read.csv(here("output/results/first_split_age_cs.csv"))
first_split_age_hipstr <- read.csv(here("output/results/first_split_age_hipstr.csv"))

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

# posterior probability of the first split in the mcc, consensus and hipstr trees
prob_first_split_summary = read.csv(
  here("output/results/marginal_prob_first_split_mcc_consensus_hipstr.csv")
)

# process data -----------------------------------------------------------------

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

# proportion of true,false and uncertain nodes in the consensus tree (true -> consensus)
prop_true_to_cs <- true_false_uncertain |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = mean(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

write_csv(prop_true_to_cs, here("output/results/prop_true_to_cs.csv"))

# proportion of true, false node from the mcc to the true tree
prop_mcc_to_true <- resume_to_true_grouped |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n)/N_sim, .groups = "drop") 

write_csv(prop_mcc_to_true, here("output/results/prop_mcc_to_true.csv"))

# proportion of true, false node from the consensus to the true tree
prop_cs_to_true <- resume_to_true_grouped |>
  filter(type == "consensus") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n)/N_sim, .groups = "drop") |> 
  ungroup()

write_csv(prop_cs_to_true, here("output/results/prop_cs_to_true.csv"))

# number of true nodes in the hipstr tree
prop_hipstr_to_true <- resume_to_true_grouped |>
  filter(type == "hipstr") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n)/N_sim, .groups = "drop") |> 
  ungroup()

write_csv(prop_hipstr_to_true, here("output/results/prop_hipstr_to_true.csv"))

# Robinson-Foulds metrics ------------------------------------------------------
rf <- read.csv(here('output/results/rf_values.csv'))

rf_results <- purrr::map(1:17, ~ {
  test <- rf |>
    filter(tree_age == .x) |>
    select(RF_mean)
  list(
    lower = HDInterval::hdi(test, credMass = 0.90)[1],
    upper = HDInterval::hdi(test, credMass = 0.90)[2]
  )
})
rf_lower <- map_dbl(rf_results, "lower")
rf_upper <- map_dbl(rf_results, "upper")


rf_ic = rf |>
  group_by(tree_age) |>
  summarise(
    RF_mean = mean(RF_mean),
    RF_min = mean(RF_min),
    RF_max = mean(RF_max)
  )

rf_ic$inf <- rf_lower
rf_ic$sup <- rf_upper
rf_ic

write.csv(rf_ic, here("output/results/rf_hdi.csv"))

# mcc - true false reconstruction and marginale probability
mcc_to_true_TF <- read.csv(here("output/results/resume_to_true_TF.csv")) |> 
  filter(type == 'mcc')

mcc_reconstruction_proba_first_split = mcc_to_true_TF |>
  inner_join(prob_first_split_summary, by = c("age", "simulation")) |>
  filter(node == node_mcc ) |>
  full_join(first_split_age_mcc, by = c("age","simulation")) |>
  rename(y = N_nodes) |>
  select(age, simulation, node, y, mcc_prob, first_split) |> 
  # delete the tree for which the first split is a leaf
  filter(!is.na(y))

write_csv(mcc_reconstruction_proba_first_split, here("output/results/mcc_reconstruction_proba_first_split.csv"))


# consensus - true false reconstruction and marginale probability
resume_to_true_TF_cs <- read.csv(here("output/results/resume_to_true_TF.csv")) |> 
  filter(type == 'consensus')

cs_reconstruction_proba_first_split = resume_to_true_TF_cs |>
  inner_join(prob_first_split_summary, by = c("age", "simulation")) |>
  filter(node == node_cs ) |>
  full_join(first_split_age_cs, by = c("age","simulation")) |>
  rename(y = N_nodes) |>
  select(age, simulation, node, y, cs_prob, first_split) |> 
  # delete the tree for which the first split is a leaf
  filter(!is.na(y))

write_csv(cs_reconstruction_proba_first_split, here("output/results/cs_reconstruction_proba_first_split.csv"))


# hipstr - true false reconstruction and marginale probability
resume_to_true_TF_hipstr <- read.csv(here("output/results/resume_to_true_TF.csv")) |> 
  filter(type == 'hipstr')

hipstr_reconstruction_proba_first_split = resume_to_true_TF_hipstr |>
  inner_join(prob_first_split_summary, by = c("age", "simulation")) |>
  filter(node == node_hipstr ) |>
  full_join(first_split_age_hipstr, by = c("age","simulation")) |>
  rename(y = N_nodes) |>
  select(age, simulation, node, y, hipstr_prob, first_split) |> 
  # delete the tree for which the first split is a leaf
  filter(!is.na(y))

write_csv(hipstr_reconstruction_proba_first_split, here("output/results/hipstr_reconstruction_proba_first_split.csv"))


