# Load packages
library(here)
library(tidyverse)
# library(dplyr)
# library(jsonlite)
# library("readr")

# Config
path_compute_bounds_function <- here("src/compute_bounds/functions.R")
path_to_config <- here("src/compute_bounds/config.json")
source(path_compute_bounds_function)

# Data
data_st <- get_tree_par_fun(path_to_config, "sino-tibetan", n_tree = 1)
data_iecor <- get_tree_par_fun(path_to_config, "iecor", n_tree = 1)
data_bantu <- get_tree_par_fun(path_to_config, "bantu", n_tree = 1)

# Topology function
f_topology_st <- data_st$f_topology
f_topology_bantu <- data_bantu$f_topology
f_topology_iecor <- data_iecor$f_topology

# Root function
f_root_st <- data_st$f_root
f_root_bantu <- data_bantu$f_root
f_root_iecor <- data_iecor$f_root

# Bantu subset
data_bantu_sub <- get_tree_par_fun(path_to_config, "bantu_subsample")
data_bantu_sub2 <- get_tree_par_fun(path_to_config, "bantu_subsample_2")

f_topology_bantu_sub <- data_bantu_sub$f_topology
f_topology_bantu_sub2 <- data_bantu_sub2$f_topology

f_root_bantu_sub <- data_bantu_sub$f_root
f_root_bantu_sub2 <- data_bantu_sub2$f_root

# Valeur issue de Tracer
data_st$param$t
data_iecor$param$t
data_bantu$param$t


dt_params <- list("Sino-Tibetan" = data_st$param, "Bantu" = data_bantu$param, "Indo-European" = data_iecor$param)
bounds_tb <- dt_params |>
  map_df(
    ~ tibble(
      t = .x$t, k = .x$k, N = .x$n,
      DT = compute_upper_bound_topology(t, k, .x$Q, N),
      DR = compute_upper_bound_root(t, .x$Q, N),
      inf_topo = find_t_value(k, .x$Q, N),
      inf_root = find_t_value_root(.x$Q, N)
    )
  ) |>
  mutate(
    family = names(dt_params),
    "Substitution model" = "CTMC",
    "Clock model" = "strict",
    "Tree model" = "BD",
    .before = t
  )

write_csv(bounds_tb, here("output/results/bounds_tb.csv"))

# Bounds values
# Sino tibetan
compute_upper_bound_topology(
  data_st$param$t,
  data_st$param$k,
  data_st$param$Q,
  data_st$param$n
)

# Bantu
compute_upper_bound_topology(
  data_bantu$param$t,
  data_bantu$param$k,
  data_bantu$param$Q,
  data_bantu$param$n
)

compute_upper_bound_topology(
  data_iecor$param$t,
  data_iecor$param$k,
  data_iecor$param$Q,
  data_iecor$param$n
)

# You can do the same by using the function compute_upper_bound_root

t_values <- seq(0, 20, length.out = 100)
bounds_byt_tb <- tibble(
  t = rep(t_values, 3),
  "Delta_T" = c(f_topology_st(t_values), f_topology_bantu(t_values), f_topology_iecor(t_values)),
  "Delta_R" = c(f_root_st(t_values), f_root_bantu(t_values), f_root_iecor(t_values)),
  family = rep(c("Sino-tibetan", "Bantu", "Indo-European"), each = 100)
) |>
  pivot_longer(-c(t, family), names_to = "Delta") |>
  mutate(Delta = str_remove(Delta, "Delta_"))
write_csv(bounds_byt_tb, here("output/results/bounds_byt_tb.csv"))

# ---------------- Infima --------------------------------

# For topology
# Sino-tibetain
inf_topology_st <- find_t_value(data_st$param$k, data_st$param$Q, data_st$param$n)
inf_topology_st

# Bantu
inf_topology_bantu <- find_t_value(data_bantu$param$k, data_bantu$param$Q, data_bantu$param$n)
inf_topology_bantu

# Iecor
inf_topology_iecor <- find_t_value(data_iecor$param$k, data_iecor$param$Q, data_iecor$param$n)
inf_topology_iecor


# For root
# Sino-tibetain
inf_root_st <- find_t_value_root(data_st$param$Q, data_st$param$n)
inf_root_st

# Bantu
inf_root_bantu <- find_t_value_root(data_bantu$param$Q, data_bantu$param$n)
inf_root_bantu

# Iecor
inf_root_iecor <- find_t_value_root(data_iecor$param$Q, data_iecor$param$n)
inf_root_iecor
