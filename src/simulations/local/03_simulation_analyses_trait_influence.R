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

# load the data ---------
# frequency of good reconstruction of all the nodes in consensus tree (true -> summary)
# the value of node represent the node in the true tree
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
# with special labels for the plot
# obtain 3 points as age = 8
# for 6000 traits
count_true_to_cs_6000 <- true_false_uncertain_6000 |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = sum(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

saveRDS(count_true_to_cs_6000, here("output/results/count_true_to_cs_6000.rds"))

# for 12000 traits
count_true_to_cs_12000 <- true_false_uncertain_12000 |>
  count(age, simulation, value) |>
  group_by(age, value) |>
  summarise(mean_n = sum(n)/N_sim, .groups = "drop") |>
  mutate(value = factor(value, levels = c("0", "2", "1"))) |>
  arrange(age, value) 

saveRDS(count_true_to_cs_12000, here("output/results/count_true_to_cs_12000.rds"))

# for the mcc tree count the number of true, false
# for each summary tree, age, simulation this dataframe indicates the proprtions 
# of true and false nodes
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

saveRDS(count_true_to_mcc_12000, here("output/results/count_true_to_mcc_12000.csv"))


# count the numer of true node from de consensus tree to the true tree
# for 6000 traits
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

# for 12000 traits
count_cs_to_true_12000 <- resume_to_true_grouped_12000 |>
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

write_csv(count_cs_to_true_12000, here("output/results/count_cs_to_true_12000.csv"))


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
