library(here)
library(tidyverse)
library(TreeTools)
library(kableExtra)
library(tracerer)

# Compute the upper bound of the probability of inferring the true tree topology
compute_upperbound_DT <- function(k, N, q, t) {
  k * N * exp(-q * t)
}

# Compute the upper bound of the probability of correctly inferring ancestral states
compute_upperbound_DS <- function(pi0, pi1, N, q, t) {
  max(pi0, pi1) + N * exp(-q * t)
}

# Compute the time threshold beyond which the upper bound of the probability
# of inferring the true tree topology falls below 1
compute_inf_t_DT <- function(k, N, q, t, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DT(k, N, q, t) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}

# Compute the time threshold beyond which the upper bound of the probability
# of correctly inferring ancestral states falls below 1
compute_inf_t_DS <- function(pi0, pi1, N, q, t, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DS(pi0, pi1, N, q, t) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}

# Get the number of taxa and traits from a nexus file
get_nexus_parameters <- function(file) {
  phydt <- ReadAsPhyDat(file)
  tibble(N = length(attributes(phydt)$names), k = length(attributes(phydt)$index))
}

# Get the values of pi0, pi1, the number of generated trees, and compute q
# from a BEAST .log file
get_tracerlog_parameters <- function(file, burnin = 0.2) {
  beast_log_full <- parse_beast_tracelog_file(file)
  beast_log <- remove_burn_ins(beast_log_full, burn_in_fraction = burnin)
  beast_log |>
    select(starts_with("freqParameter"), TreeHeight.t.tree) |>
    summarise(across(everything(), ~ mean(.x))) |>
    rename_all(str_replace, pattern = "freq.+(\\d)", replacement = "pi\\1") |>
    rename(pi0 = pi1, pi1 = pi2, t_R = TreeHeight.t.tree) |>
    mutate(q = 1 / (pi0^2 + pi1^2)) |>
    mutate(nTrees = max(beast_log$Sample)) |>
    relocate(q, .after = pi1)
}

# Combine all of the above
get_all_parameters <- function(logfile, nexusfile, burnin = 0.2, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  bind_cols(
    get_nexus_parameters(nexusfile),
    get_tracerlog_parameters(logfile)
  ) |>
    relocate(nTrees, .before = pi0) |>
    mutate(ub_DT = compute_upperbound_DT(k, N, q, t_R)) |>
    mutate(ub_DS = compute_upperbound_DS(pi0, pi1, N, q, t_R)) |>
    mutate(inf_t_DT = compute_inf_t_DT(k, N, q, t_R)) |>
    mutate(inf_t_DS = compute_inf_t_DS(pi0, pi1, N, q, t_R))
}

# Run get_all_parameters on the files within each folder
dt_real <- list.dirs(here("data/real"), full.names = TRUE, recursive = FALSE) |>
  map_df(function(x) {
    d <- str_remove_all(x, ".*/")
    logfile <- list.files(x, "\\.log", full.names = TRUE)
    nexusfile <- list.files(x, "\\.nex", full.names = TRUE)
    #burnin <- ifelse(str_detect(x, "^tea"), .8, .2)
    burnin <- 0.2
    if (length(logfile) > 0 & length(nexusfile) > 0) {
      bind_cols(tibble(d), get_all_parameters(logfile, nexusfile, burnin = burnin))
    } else {
      tibble(d)
    }
  }) |>
  filter(!is.na(N)) |>
  mutate(d = case_when(
    str_detect(d, "^bantu.+subsample$") ~ "Bantu subset",
    str_detect(d, "^bantu.+subsample2$") ~ "Bantu subset 2",
    str_detect(d, "^bantu") ~ "Bantu",
    str_detect(d, "^ie") ~ "Indo-European",
    str_detect(d, "^st.+tibetanPrior$") ~ "Sino-Tibetan prior",
    str_detect(d, "^st") ~ "Sino-Tibetan",
    str_detect(d, "^tea") ~ "Trans-Eurasian",
  )) |>
  rename(family = d) |> 
  group_by(family) |> 
  slice(1)
write_csv(dt_real, here("output/results/bounds_real_tb.csv"))


dt_sim <- list.dirs(here("data/simulated"), full.names = TRUE, recursive = FALSE) |>
  map_df(function(x) {
    d <- str_remove_all(x, ".*/")
    logfile <- list.files(x, "\\.log", full.names = TRUE)
    nexusfile <- here("src/preprocess_simulated_data/tree-sim.nex")
    bind_cols(tibble(d), get_all_parameters(logfile, nexusfile))
  }) |> 
  mutate(d = str_remove_all(d, "[^0-9]") |> as.integer()) |>
  rename(family = d) |> 
  arrange(family)


# simtr <- treeio::read.tree(here("data/simulated/beast-data-sim-10/tree-sim-10.tree"))
# max(castor::get_all_pairwise_distances(simtr)[, Ntip(simtr)+1], na.rm = TRUE)
# 
# read.nexus(here("data/simulated/beast-data-sim-1/ctmc-strict-bd-1.trees"))

t_values <- seq(0, 20, length.out = 101)
bounds_real_byt_tb <- dt_real |>
  group_by(family) |>
  slice(1) |>
  mutate(count = length(t_values)) |>
  uncount(count) |>
  mutate(t = t_values) |>
  rowwise() |>
  mutate(ub_DT = compute_upperbound_DT(k, N, q, t)) |>
  mutate(ub_DS = compute_upperbound_DS(pi0, pi1, N, q, t))
write_csv(bounds_real_byt_tb, here("output/results/bounds_real_byt_tb.csv"))

