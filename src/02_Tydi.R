library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)


# Get the number of taxa and traits from a nexus file
get_nexus_parameters <- function(file) {
  phydt <- ReadAsPhyDat(file)
  tibble(N = length(attributes(phydt)$names), k = length(attributes(phydt)$index))
}


# Get the values of pi0, pi1, the number of generated trees, and compute q
# from a BEAST .log file
get_tracerlog_parameters <- function(file, burnin = 0.2) {
  beast_log_full <- read.csv(file)
  beast_log <- remove_burn_ins(beast_log_full, burn_in_fraction = burnin)
  beast_log |>
    select(starts_with("freqParameter"), TreeHeight.t.tree) |>
    summarise(across(everything(), ~ median(.x))) |>
    mutate(nTrees = max(beast_log$Sample)) 
}


# Sino-tibetain
st_file = here("output/results/st/st_ctmc-strict-fbd_tracelog.csv")
st_par <- get_tracerlog_parameters(st_file) %>%
  select(matches("\\d$"), TreeHeight.t.tree, nTrees) %>%
  pivot_longer(cols = -c(TreeHeight.t.tree, nTrees), names_to = "Column", values_to = "Value") %>%
  mutate(Group = str_extract(Column, "\\d$")) %>%
  group_by(Group) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean_Value) %>%
  rename(pi0 = "1", pi1 = "2") %>%
  bind_cols(get_tracerlog_parameters(st_file) %>% select(TreeHeight.t.tree, nTrees)) %>%
  mutate(q = 1 / (pi0^2 + pi1^2)) %>%
  relocate(q, .after = pi1) %>%
  rename( t_R = TreeHeight.t.tree)

# Transeurasien
tea_file = here("output/results/tea/tea_ctmc-strict-fbd-constrained_tracelog.csv")
tea_par <- get_tracerlog_parameters(tea_file) %>%
  select(matches("\\d$"), TreeHeight.t.tree, nTrees) %>%
  pivot_longer(cols = -c(TreeHeight.t.tree, nTrees), names_to = "Column", values_to = "Value") %>%
  mutate(Group = str_extract(Column, "\\d$")) %>%
  group_by(Group) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean_Value) %>%
  rename(pi0 = "1", pi1 = "2") %>%
  bind_cols(get_tracerlog_parameters(tea_file) %>% select(TreeHeight.t.tree, nTrees)) %>%
  mutate(q = 1 / (pi0^2 + pi1^2)) %>%
  relocate(q, .after = pi1) %>%
  rename( t_R = TreeHeight.t.tree) %>%
  mutate(t_R = t_R*100)


# Bantu
bantu_file = here("output/results/bantu/bantu_ctmc-strict-bd_tracelog.csv")
bantu_par <- get_tracerlog_parameters(bantu_file) %>%
  select(matches("\\d$"), TreeHeight.t.tree, nTrees) %>%
  pivot_longer(cols = -c(TreeHeight.t.tree, nTrees), names_to = "Column", values_to = "Value") %>%
  mutate(Group = str_extract(Column, "\\d$")) %>%
  group_by(Group) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean_Value) %>%
  rename(pi0 = "1", pi1 = "2") %>%
  bind_cols(get_tracerlog_parameters(bantu_file) %>% select(TreeHeight.t.tree, nTrees)) %>%
  mutate(q = 1 / (pi0^2 + pi1^2)) %>%
  relocate(q, .after = pi1) %>%
  rename( t_R = TreeHeight.t.tree) %>%
  mutate(t_R = t_R*100)


# Bantu subsample
bantu_subsample_file = here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tracelog.csv")
bantu_subsample_par <- get_tracerlog_parameters(bantu_subsample_file) %>%
  select(matches("\\d$"), TreeHeight.t.tree, nTrees) %>%
  pivot_longer(cols = -c(TreeHeight.t.tree, nTrees), names_to = "Column", values_to = "Value") %>%
  mutate(Group = str_extract(Column, "\\d$")) %>%
  group_by(Group) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean_Value) %>%
  rename(pi0 = "1", pi1 = "2") %>%
  bind_cols(get_tracerlog_parameters(bantu_subsample_file) %>% select(TreeHeight.t.tree, nTrees)) %>%
  mutate(q = 1 / (pi0^2 + pi1^2)) %>%
  relocate(q, .after = pi1) %>%
  rename( t_R = TreeHeight.t.tree) %>%
  mutate(t_R = t_R*100)

# Bantu subsample 2
bantu_subsample2_file = here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tracelog.csv")
bantu_subsample2_par <- get_tracerlog_parameters(bantu_subsample2_file) %>%
  select(matches("\\d$"), TreeHeight.t.tree, nTrees) %>%
  pivot_longer(cols = -c(TreeHeight.t.tree, nTrees), names_to = "Column", values_to = "Value") %>%
  mutate(Group = str_extract(Column, "\\d$")) %>%
  group_by(Group) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean_Value) %>%
  rename(pi0 = "1", pi1 = "2") %>%
  bind_cols(get_tracerlog_parameters(bantu_subsample2_file) %>% select(TreeHeight.t.tree, nTrees)) %>%
  mutate(q = 1 / (pi0^2 + pi1^2)) %>%
  relocate(q, .after = pi1) %>%
  rename( t_R = TreeHeight.t.tree) %>%
  mutate(t_R = t_R*100)




