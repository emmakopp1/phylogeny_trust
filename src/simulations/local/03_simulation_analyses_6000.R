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
N_sim <- 50

# load data --------------------------------------------------------------------
# marginal probability of the first split with IC 
marginal_probability_first_split_ic_6000 <- read.csv(
  here("output/results/marginal_prob_first_split_ic_1_50_6000.csv")
  )


# frequency of good reconstruction of all the nodes in of the mcc
mcc_to_true_TF_6000 <- read.csv(here("output/results/resume_to_true_TF_1_50_6000.csv")) |> 
  filter(type == 'mcc')

# for each tree simulation, age between the root and the first split
first_split_age_6000 <- read.csv(here("output/results/first_split_age_6000.csv"))

# frequency of good reconstruction of all the nodes in consensus tree (true -> summary)
# the value of node represent the node in the true tree
true_false_uncertain_6000 <- read.csv(
  file = here("output/results/true_false_uncertain_nodes_1_50_6000.csv"),
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
df_number_of_nodes_6000 <- read.csv(
  file = here("output/results/number_nodes_mcc_cs_6000.csv"),
  sep = ",",
  header = T
)

# marginal probability of the first split in the mcc and consensus tree
prob_first_split_summary_6000 = read.csv(
  here("output/results/marginal_prob_first_split_mcc_consensus_6000.csv")
  )

# process data 6000 -------------------------------------------------------------

# marginal probability of the first split in the mcc 
# obtain one point as age = 8
prob_first_split_mcc_6000 = prob_first_split_summary_6000 |> 
  select(- cs_prob, - node_cs, - node_mcc) |> 
  group_by(age) |> 
  summarise(mean_mcc_prob = mean(mcc_prob, na.rm=T), .groups='drop') |> 
  ungroup() 

write.csv(prob_first_split_mcc_6000,
          here("output/results/prob_first_split_mcc_6000.csv"), row.names = FALSE)

# number of node in the summary tree
# obtain one point as age = 8
df_number_of_nodes_avg_6000 <- df_number_of_nodes_6000 |>
  group_by(age) |>
  summarise(
    n_mcc = mean(n_mcc, na.rm = TRUE),
    n_consensus = mean(n_consensus, na.rm = TRUE)
  ) 

write.csv(
  df_number_of_nodes_avg_6000,
  here("output/results/number_of_nodes_summary_6000.csv"), row.names = FALSE)

# posterior of the first split with IC
# obtain one point as age = 8
marginal_probability_first_split_ic_6000 = marginal_probability_first_split_ic_6000|>
  rename(age = tree_age, simulation = tree_simulation_number) |>
  group_by(age) |>
  summarise(
    prob = mean(prob_mean, na.rm = TRUE),
    inf = mean(prob_inf, na.rm = TRUE),
    sup = mean(prob_sup, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  marginal_probability_first_split_ic_6000, 
  here("output/results/marginal_probability_first_split_ic_6000.csv"))

# for the consensus trees, count the number of true, false and uncertain nodes
# with special labels for the plot
# obtain 3 points as age = 8
count_true_to_cs_6000 <- true_false_uncertain_6000 |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = sum(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

saveRDS(count_true_to_cs_6000, here("output/results/count_true_to_cs_6000.rds"))

# for the mcc tree count the number of true, false
# for each summary tree, age, simulation this dataframe indicates the proprtions 
# of true and false nodes
resume_to_true_grouped_6000 <- read_csv(here("output/results/resume_to_true_TF_1_50_6000.csv"), col_names = T)|> 
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

saveRDS(count_true_to_mcc_6000,here("output/results/count_true_to_mcc_6000.csv"))


# count the numer of true node from de consensus tree to the true tree
count_cs_to_true_6000 <- resume_to_true_grouped_6000 |>
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

write_csv(count_cs_to_true_6000, here("output/results/count_cs_to_true_6000.csv"))

# proportion of true,false and uncertain nodes in the consensus tree (true -> consensus)
prop_true_to_cs_6000 <- true_false_uncertain_6000 |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = sum(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

write_csv(prop_true_to_cs_6000, here("output/results/prop_true_to_cs_6000.csv"))

# proportion of true, false node from the mcc to the true tree
prop_mcc_to_true_6000 <- resume_to_true_grouped_6000 |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(n_mean = sum(n)/N_sim, .groups = "drop") 

write_csv(prop_mcc_to_true_6000, here("output/results/prop_mcc_to_true_6000.csv"))

# proportion of true, false node from the consensus to the true tree
prop_cs_to_true_6000 <- resume_to_true_grouped_6000 |>
  filter(type == "consensus") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n), .groups = "drop") |> 
  ungroup()

write_csv(prop_cs_to_true_6000, here("output/results/prop_cs_to_true_6000.csv"))

# number of true nodes in the summary tree
prop_mcc_to_true_6000 <- resume_to_true_grouped_6000 |>
  filter(type == "mcc") |>
  group_by(age, exist) |>
  summarise(n_mean = mean(n), .groups = "drop") |> 
  ungroup()

write_csv(prop_mcc_to_true_6000, here("output/results/prop_mcc_to_true_6000.csv"))

# analysis and comparison between number of traits -----------------------------

# --- Données prob_first_split_mcc ---
prob_first_split_mcc_data <- bind_rows(
  read.csv(here("output/results/prob_first_split_mcc_6000.csv")) |> mutate(n_trait = 6000),
  read.csv(here("output/results/prob_first_split_mcc.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = mean_mcc_prob,
    names_prefix = "prob_mcc_"
  ) |>
  select(age, prob_mcc_3000, prob_mcc_6000)


# --- Données number_of_nodes_summary ---
number_of_nodes_summary_data <- bind_rows(
  read.csv(here("output/results/number_of_nodes_summary_6000.csv")) |> mutate(n_trait = 6000),
  read.csv(here("output/results/number_of_nodes_summary.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = c(n_mcc, n_consensus), # Pivoter les deux colonnes n_mcc et n_consensus
    names_prefix = "" # Pas de préfixe global, car les noms des valeurs sont déjà clairs (n_mcc_, n_consensus_)
  ) |>
  # Renommer spécifiquement pour la clarté (ex: n_mcc_3000, n_consensus_3000)
  rename(
    n_mcc_3000 = `n_mcc_3000`,
    n_mcc_6000 = `n_mcc_6000`,
    n_consensus_3000 = `n_consensus_3000`,
    n_consensus_6000 = `n_consensus_6000`
  ) |>
  select(age, n_mcc_3000, n_mcc_6000, n_consensus_3000, n_consensus_6000)


# --- Données marginal_probability_first_split_ic ---
marginal_prob_ic_data <- bind_rows(
  read_csv(here("output/results/marginal_probability_first_split_ic_6000.csv")) |> mutate(n_trait = 6000),
  read_csv(here("output/results/marginal_probability_first_split_ic.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = c(prob, inf, sup), # Pivoter prob, inf, sup
    names_prefix = ""
  ) |>
  rename(
    prob_3000 = 'prob_3000', prob_6000 = 'prob_6000',
    inf_3000 = 'inf_3000', inf_6000 = 'inf_6000',
    sup_3000 = 'sup_3000', sup_6000 = 'sup_6000'
  ) |>
  select(age, prob_3000, prob_6000, inf_3000, inf_6000, sup_3000, sup_6000)


# --- Données count_true_to_cs ---
count_true_to_cs_data <- bind_rows(
  readRDS(here("output/results/count_true_to_cs_6000.rds")) |> mutate(n_trait = 6000),
  readRDS(here("output/results/count_true_to_cs.rds")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = mean_n,
    names_prefix = "mean_n_"
  ) |>
  select(age, value, mean_n_3000, mean_n_6000) |>
  arrange(value)


# --- Données count_true_to_mcc ---
count_true_to_mcc_data <- bind_rows(
  readRDS(here("output/results/count_true_to_mcc_6000.csv")) |> mutate(n_trait = 6000),
  readRDS(here("output/results/count_true_to_mcc.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = n_mean,
    names_prefix = "n_mean_"
  ) |>
  select(age, exist, n_mean_3000, n_mean_6000) |>
  arrange(exist)


# --- Données count_cs_to_true ---
count_cs_to_true_data <- bind_rows(
  read_csv(here("output/results/count_cs_to_true_6000.csv")) |> mutate(n_trait = 6000),
  read_csv(here("output/results/count_cs_to_true.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    id_cols = c(age, exist),
    names_from = n_trait,
    values_from = c(n_mean, total), # Pivot both n_mean and total
    names_prefix = ""
  ) |>
  rename(
    n_mean_3000 = "n_mean_3000", n_mean_6000 = "n_mean_6000",
    total_3000 = "total_3000", total_6000 = "total_6000"
  ) |>
  select(age, exist, n_mean_3000, n_mean_6000, total_3000, total_6000) |>
  arrange(age, exist)

# visualization ----------------------------------------------------------------
# currently here because 12000 analysis coming 

# from true to cs 
count_true_to_cs_data_long <- count_true_to_cs_data |>
  pivot_longer(
    cols = starts_with("mean_n_"),
    names_to = "n_trait_col",
    values_to = "mean_n"
  ) |>
  mutate(
    n_trait = as.numeric(gsub("mean_n_", "", n_trait_col)), # Extraire le nombre de traits
    value = as.factor(value) # Assurez-vous que 'value' est un facteur pour l'esthétique de remplissage
  )

# from true to mcc
count_true_to_mcc_data_long <- count_true_to_mcc_data |>
  pivot_longer(
    cols = starts_with("n_mean_"),
    names_to = "n_trait_col",
    values_to = "mean_n"
  ) |>
  mutate(
    n_trait = parse_number(n_trait_col),  
    exist = as.factor(exist),
    n_trait_col = str_replace(n_trait_col, "n_mean_", "mean_n_") 
  )


p1<- ggplot(count_true_to_cs_data_long, aes(x = as.factor(n_trait), y = mean_n, fill = value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "",
    x = "number of traits",
    y = "average number of nodes",
    fill = "node category"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue", "2" = "darkorange"),
    breaks = c("0", "2", "1"),
    labels = c("false", "plausible", "true")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


p2<- ggplot(count_true_to_mcc_data_long, aes(x = as.factor(n_trait), y = mean_n, fill = exist)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "",
    x = "number of traits",
    y = "average number of nodes",
    fill = "node category"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("false", "true")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p1+p2




