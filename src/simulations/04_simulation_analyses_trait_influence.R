# ------------------------------------------------------------------------------
# Script Name: 03_simulation_analyses.R
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
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidyr)
library(HDInterval)
N_sim <- 50

# load the data ---------
# frequency of good reconstruction of all the nodes in consensus tree (true -> summary)
# the value of node represent the node in the true tree
true_false_uncertain_1500 <- read.csv(
  file = here("output/results/true_false_uncertain_nodes_1500.csv"),
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

true_false_uncertain_3000 <- read.csv(
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

true_false_uncertain_6000 <- read.csv(
  file = here("output/results/true_false_uncertain_nodes_6000.csv"),
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

true_false_uncertain_12000 <- read.csv(
  file = here("output/results/true_false_uncertain_nodes_12000.csv"),
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

# process the data ---------
# for the consensus trees, count the number of true, false and uncertain nodes
# averaged by the number of simulations
# obtain 3 points as age = 8
# for 1500 traits
count_true_to_cs_1500 <- true_false_uncertain_1500 |>
  count(age, simulation, value) |>
  group_by(age, simulation) |>
  mutate(prop = n / sum(n)) |>          # frequency per simulation
  group_by(age, value) |>
  summarise(prop = mean(prop), .groups = "drop") |>   # mean across simulations
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value)

saveRDS(count_true_to_cs_1500, here("output/results/prop_true_to_cs_1500.rds"))

# for 3000 traits
count_true_to_cs_3000 <- true_false_uncertain_3000 |>
  count(age, simulation, value) |>
  group_by(age, simulation) |>
  mutate(prop = n / sum(n)) |>          # frequency per simulation
  group_by(age, value) |>
  summarise(prop = mean(prop), .groups = "drop") |>   # mean across simulations
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value)

saveRDS(count_true_to_cs_3000, here("output/results/prop_true_to_cs_3000.rds"))

# for 6000 traits
count_true_to_cs_6000 <-true_false_uncertain_6000 |>
  count(age, simulation, value) |>
  group_by(age, simulation) |>
  mutate(prop = n / sum(n)) |>          # frequency per simulation
  group_by(age, value) |>
  summarise(prop = mean(prop), .groups = "drop") |>   # mean across simulations
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value)

saveRDS(count_true_to_cs_6000, here("output/results/prop_true_to_cs_6000.rds"))

# for 12000 traits
count_true_to_cs_12000 <- true_false_uncertain_12000 |>
  count(age, simulation, value) |>
  group_by(age, simulation) |>
  mutate(prop = n / sum(n)) |>          # frequency per simulation
  group_by(age, value) |>
  summarise(prop = mean(prop), .groups = "drop") |>   # mean across simulations
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value)

saveRDS(count_true_to_cs_12000, here("output/results/prop_true_to_cs_12000.rds"))

# for the mcc tree count the number of true, false
# for each summary tree, age, simulation this dataframe indicates the proprtions 
# of true and false nodes
# for 1500 traits
resume_to_true_grouped_1500 <- read_csv(here("output/results/resume_to_true_TF_1500.csv"), col_names = T)|> 
  rename(N_node = exist, exist = N_nodes) |>
  mutate(exist = as.numeric(exist)) |>  # TRUE -> 1, FALSE -> 0
  group_by(type, age, simulation, exist, N_node) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(type, age, simulation) |>
  mutate(proportion = n / N_node) |>
  ungroup() |>
  mutate(exist = factor(exist, levels = c("0", "1")))

# number of well reconstructed node in the mcc tree
count_true_to_mcc_1500 <- resume_to_true_grouped_1500 |>
  #filter(simulation==1) |> 
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(proportion = mean(proportion), .groups = "drop") |>
  tidyr::pivot_wider(names_from = exist, values_from = proportion, names_prefix = "exist_") |>
  mutate(
    total = exist_0 + exist_1
  ) |>
  pivot_longer(cols = starts_with("exist_"), names_prefix = "exist_", names_to = "exist", values_to = "n_proportion") |>
  mutate(
    exist = factor(exist, levels = c("0", "1")),
    y_label = ifelse(exist == "1", 0, total) # manual alignment
  ) 

saveRDS(count_true_to_mcc_1500,here("output/results/prop_true_to_mcc_1500.csv"))

# for 3000 traits
resume_to_true_grouped_3000 <- read_csv(here("output/results/resume_to_true_TF.csv"), col_names = T)|> 
  rename(N_node = exist, exist = N_nodes) |>
  mutate(exist = as.numeric(exist)) |>  # TRUE -> 1, FALSE -> 0
  group_by(type, age, simulation, exist, N_node) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(type, age, simulation) |>
  mutate(proportion = n / N_node) |>
  ungroup() |>
  mutate(exist = factor(exist, levels = c("0", "1")))

# number of well reconstructed node in the mcc tree
count_true_to_mcc_3000 <- resume_to_true_grouped_3000 |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(proportion = mean(proportion), .groups = "drop") |>
  tidyr::pivot_wider(names_from = exist, values_from = proportion, names_prefix = "exist_") |>
  mutate(
    total = exist_0 + exist_1
  ) |>
  pivot_longer(cols = starts_with("exist_"), names_prefix = "exist_", names_to = "exist", values_to = "n_proportion") |>
  mutate(
    exist = factor(exist, levels = c("0", "1")),
    y_label = ifelse(exist == "1", 0, total) # manual alignment
  ) 

saveRDS(count_true_to_mcc_3000,here("output/results/prop_true_to_mcc_3000.csv"))

# for 6000 traits
resume_to_true_grouped_6000 <- read_csv(here("output/results/resume_to_true_TF_6000.csv"), col_names = T)|> 
  rename(N_node = exist, exist = N_nodes) |>
  mutate(exist = as.numeric(exist)) |>  # TRUE -> 1, FALSE -> 0
  group_by(type, age, simulation, exist, N_node) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(type, age, simulation) |>
  mutate(proportion = n / N_node) |>
  ungroup() |>
  mutate(exist = factor(exist, levels = c("0", "1")))

# number of well reconstructed node in the mcc tree
count_true_to_mcc_6000 <- resume_to_true_grouped_6000 |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(proportion = mean(proportion), .groups = "drop") |>
  tidyr::pivot_wider(names_from = exist, values_from = proportion, names_prefix = "exist_") |>
  mutate(
    total = exist_0 + exist_1
  ) |>
  pivot_longer(cols = starts_with("exist_"), names_prefix = "exist_", names_to = "exist", values_to = "n_proportion") |>
  mutate(
    exist = factor(exist, levels = c("0", "1")),
    y_label = ifelse(exist == "1", 0, total) # manual alignment
  ) 

saveRDS(count_true_to_mcc_6000,here("output/results/prop_true_to_mcc_6000.csv"))

# for 12000 traits
resume_to_true_grouped_12000 <- read_csv(here("output/results/resume_to_true_TF_12000.csv"), col_names = T)|> 
  rename(N_node = exist, exist = N_nodes) |>
  mutate(exist = as.numeric(exist)) |>  # TRUE -> 1, FALSE -> 0
  group_by(type, age, simulation, exist, N_node) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(type, age, simulation) |>
  mutate(proportion = n / N_node) |>
  ungroup() |>
  mutate(exist = factor(exist, levels = c("0", "1")))

# number of well reconstructed node in the mcc tree
count_true_to_mcc_12000 <- resume_to_true_grouped_12000 |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(proportion = mean(proportion), .groups = "drop") |>
  tidyr::pivot_wider(names_from = exist, values_from = proportion, names_prefix = "exist_") |>
  mutate(
    total = exist_0 + exist_1
  ) |>
  pivot_longer(cols = starts_with("exist_"), names_prefix = "exist_", names_to = "exist", values_to = "n_proportion") |>
  mutate(
    exist = factor(exist, levels = c("0", "1")),
    y_label = ifelse(exist == "1", 0, total) # manual alignment
  ) 

saveRDS(count_true_to_mcc_12000, here("output/results/prop_true_to_mcc_12000.csv"))

# --- count_true_to_cs data ---
count_true_to_cs_data <- bind_rows(
  readRDS(here("output/results/prop_true_to_cs_12000.rds")) |> mutate(n_trait = 12000),
  readRDS(here("output/results/prop_true_to_cs_6000.rds")) |> mutate(n_trait = 6000),
  readRDS(here("output/results/prop_true_to_cs_3000.rds")) |> filter(age==8) |> mutate(n_trait = 3000),
  readRDS(here("output/results/prop_true_to_cs_1500.rds")) |> mutate(n_trait = 1500)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = prop,
    names_prefix = "prop_n_"
  ) |>
  select(age, 
         value, 
         prop_n_1500,
         prop_n_3000, 
         prop_n_6000, 
         prop_n_12000) |>
  arrange(value) |> 
  pivot_longer(
    cols = starts_with("prop_n_"),
    names_to = "n_trait_col",
    values_to = "prop_n"
  ) |>
  mutate(
    n_trait = as.numeric(gsub("prop_n_", "", n_trait_col)), # Extract the number of traits
    value = as.factor(value) # Make sure 'value' is a factor for the fill aesthetic
  )



# --- count_true_to_mcc data ---
count_true_to_mcc_data <- bind_rows(
  readRDS(here("output/results/prop_true_to_mcc_12000.csv")) |> mutate(n_trait = 12000),
  readRDS(here("output/results/prop_true_to_mcc_6000.csv")) |> mutate(n_trait = 6000),
  readRDS(here("output/results/prop_true_to_mcc_3000.csv")) |> filter(age==8) |> mutate(n_trait = 3000),
  readRDS(here("output/results/prop_true_to_mcc_1500.csv")) |> mutate(n_trait = 1500)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = n_proportion,
    names_prefix = "prop_n_"
  ) |>
  select(age, 
         exist, 
         prop_n_1500,
         prop_n_3000, 
         prop_n_6000, 
         prop_n_12000) |>
  arrange(exist) |> 
  pivot_longer(
    cols = starts_with("prop_n_"),
    names_to = "n_trait_col",
    values_to = "prop_n"
  ) |>
  mutate(
    n_trait = as.numeric(gsub("prop_n_", "", n_trait_col)), # Extract the number of traits
    value = as.factor(exist) # Make sure 'value' is a factor for the fill aesthetic
  ) |>
  select(-exist)

write_csv(count_true_to_cs_data, here("output/results/prop_true_to_cs_data_long.csv"))
write_csv(count_true_to_mcc_data, here("output/results/prop_true_to_mcc_data_long.csv"))


