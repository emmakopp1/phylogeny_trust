library(here)
library(tidyverse)
library(TreeTools)
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
dt <- list.dirs(here("data/real"), full.names = TRUE, recursive = FALSE) |>
  str_subset("/iecor_co", negate = TRUE) |>
  map_df(function(x) {
    d <- str_remove_all(x, ".*/")
    logfile <- list.files(x, "\\.log", full.names = TRUE)
    nexusfile <- list.files(x, "\\.nex", full.names = TRUE)
    if (length(logfile) > 0 & length(nexusfile) > 0) {
      bind_cols(tibble(d), get_all_parameters(logfile, nexusfile))
    } else {
      tibble(d)
    }
  })

t_values <- seq(0, 20, length.out = 101)

dt |> 
  filter(!is.na(N)) |> 
  select(d, N, k, pi0, pi1, q) |> 
  mutate(count = length(t_values)) |> 
  uncount(count) |> 
  group_by(d) |> 
  mutate(t = t_values) |> 
  mutate(ub_DT = compute_upperbound_DT(k, N, q, t)) |> 
  mutate(ub_DS = compute_upperbound_DS(pi0, pi1, N, q, t)) |> 
  filter(t == max(t))

