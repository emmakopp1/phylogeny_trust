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

# analysis and comparison between number of traits -----------------------------

# --- Données prob_first_split_mcc ---
prob_first_split_mcc_data <- bind_rows(
  read.csv(here("output/results/prob_first_split_mcc_12000.csv")) |> mutate(n_trait = 12000),
  read.csv(here("output/results/prob_first_split_mcc_6000.csv")) |> mutate(n_trait = 6000),
  read.csv(here("output/results/prob_first_split_mcc.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = mean_mcc_prob,
    names_prefix = "prob_mcc_"
  ) |>
  select(age, prob_mcc_3000, prob_mcc_6000, prob_mcc_12000)


# --- Données number_of_nodes_summary ---
number_of_nodes_summary_data <- bind_rows(
  read.csv(here("output/results/number_of_nodes_summary_12000.csv")) |> mutate(n_trait = 12000),
  read.csv(here("output/results/number_of_nodes_summary_6000.csv")) |> mutate(n_trait = 6000),
  read.csv(here("output/results/number_of_nodes_summary.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = c(n_mcc, n_consensus), # Pivoter les deux colonnes n_mcc et n_consensus
    names_prefix = "" 
  ) |>
  rename(
    n_mcc_3000 = 'n_mcc_3000',
    n_mcc_6000 = 'n_mcc_6000',
    n_mcc_12000 = 'n_mcc_12000',
    n_consensus_3000 = 'n_consensus_3000',
    n_consensus_6000 = 'n_consensus_6000',
    n_consensus_12000 = 'n_consensus_12000'
  ) |>
  select(age, n_mcc_3000, n_mcc_6000,n_mcc_12000, n_consensus_3000, n_consensus_6000,n_consensus_12000)


# --- Données marginal_probability_first_split_ic ---
marginal_prob_ic_data <- bind_rows(
  read_csv(here("output/results/marginal_probability_first_split_ic_12000.csv")) |> mutate(n_trait = 12000),
  read_csv(here("output/results/marginal_probability_first_split_ic_6000.csv")) |> mutate(n_trait = 6000),
  read_csv(here("output/results/marginal_probability_first_split_ic.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = c(prob, inf, sup), # Pivoter prob, inf, sup
    names_prefix = ""
  ) |>
  rename(
    prob_3000 = 'prob_3000', prob_6000 = 'prob_6000', prob_12000 = 'prob_12000',
    inf_3000 = 'inf_3000', inf_6000 = 'inf_6000', inf_12000 = 'inf_12000',
    sup_3000 = 'sup_3000', sup_6000 = 'sup_6000', sup_12000 = 'sup_12000'
  ) |>
  select(age, prob_3000, prob_6000, prob_12000, inf_3000, inf_6000, inf_12000, sup_3000, sup_6000, sup_12000)


# --- Données count_true_to_cs ---
count_true_to_cs_data <- bind_rows(
  readRDS(here("output/results/count_true_to_cs_12000.rds")) |> mutate(n_trait = 12000),
  readRDS(here("output/results/count_true_to_cs_6000.rds")) |> mutate(n_trait = 6000),
  readRDS(here("output/results/count_true_to_cs.rds")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = mean_n,
    names_prefix = "mean_n_"
  ) |>
  select(age, value, mean_n_3000, mean_n_6000, mean_n_12000) |>
  arrange(value)


# --- Données count_true_to_mcc ---
count_true_to_mcc_data <- bind_rows(
  readRDS(here("output/results/count_true_to_mcc_12000.csv")) |> mutate(n_trait = 12000),
  readRDS(here("output/results/count_true_to_mcc_6000.csv")) |> mutate(n_trait = 6000),
  readRDS(here("output/results/count_true_to_mcc.csv")) |> filter(age==8) |> mutate(n_trait = 3000)
) |>
  pivot_wider(
    names_from = n_trait,
    values_from = n_mean,
    names_prefix = "n_mean_"
  ) |>
  select(age, exist, n_mean_3000, n_mean_6000, n_mean_12000) |>
  arrange(exist)


# --- Données count_cs_to_true ---
count_cs_to_true_data <- bind_rows(
  read_csv(here("output/results/count_cs_to_true_12000.csv")) |> mutate(n_trait = 12000),
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
    n_mean_3000 = "n_mean_3000", n_mean_6000 = "n_mean_6000", n_mean_12000 = "n_mean_12000",
    total_3000 = "total_3000", total_6000 = "total_6000",  total_12000 = "total_12000"
  ) |>
  select(age, exist, n_mean_3000, n_mean_6000, n_mean_12000, total_3000, total_6000, total_12000) |>
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
    title = "Consensus",
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
    title = "MCC",
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

plt_number_of_traits_influence <- p1 + p2
ggsave(plt_number_of_traits_influence, filename = here("output/figs/number_of_traits_influence.png"), 
       width = 12, height = 6, dpi = 300)




