library(here)
library(tidyverse)


# Compute the upper bound of the probability of inferring the true tree topology for one sens
compute_upperbound_DT <- function(k, q, t, depth, s) {
  k * sum(exp(-q * (t - depth - s)))
}

# Compute the upper bound of the probability of correctly inferring ancestral states
compute_upperbound_DS <- function(pi0, pi1, q, t, depth) {
  max(as.numeric(pi0), as.numeric(pi1)) + sum(exp(-as.numeric(q) * (t - depth)))
}

# Compute the time threshold beyond which the upper bound of the probability
# of inferring the true tree topology falls below 1
compute_inf_t_DT <- function(k, q, depth, s, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DT(k, q, t, depth, s) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}

# Compute the time threshold beyond which the upper bound of the probability
# of correctly inferring ancestral states falls below 1
compute_inf_t_DS <- function(pi0, pi1, q, depth, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DS(pi0, pi1, q, t, depth) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}

# Add cutting point
calibrations <- read_csv(here("output/results/calibration.csv")) |>
  select(family, calibration, tip, s)

# Tipages summary
tipages_summary <- read_csv(here("output/results/tipages_summary.csv")) |>
  full_join(calibrations, by = c("family", "tip")) |>
  mutate(s = if_else(is.na(s), depth, s)) |>
  group_by(family, tip) |>
  filter(s == min(s)) |>
  ungroup()

# Tracelog summary
tracelog_summary <- read_csv(here("output/results/tracelog_summary.csv"))

# Bounds
bounds_real_tb <- tracelog_summary |>
  mutate(n_cogsets = if_else(is.na(n_cogsets), k, n_cogsets)) |>
  rename(familyx = family) |> 
  relocate(c(N, k), .after = n_trees) |>
  rowwise() |>
  mutate(
    ub_DS = compute_upperbound_DS(pi0, pi1, q,
      t = t_R,
      depth = filter(tipages_summary, family == family)$depth
    ),
    inf_t_DS = compute_inf_t_DS(pi0, pi1, q, depth = filter(tipages_summary, family == family)$depth),
    ub_DT = compute_upperbound_DT(k = n_cogsets, q = q, t = t_R, depth = filter(tipages_summary, family == familyx)$depth, s = filter(tipages_summary, family == familyx)$s)
  ) |> 
  ungroup() |> 
  rename(family = familyx)

compute_upperbound_DT(
  k = 3859, 
  q = 1.03, 
  t = 6.17, 
  depth = filter(tipages_summary, family == "Bantu")$depth, 
  s = filter(tipages_summary, family == "Bantu")$s)


# TEA add line with the mean of parameter for all cognates
TEA_summary <- bounds_real_tb %>%
  filter(family == "TEA") %>%
  select(-concept, -family, -n_cogsets, -ub_DT) %>%
  colMeans() %>%
  t() %>%
  as_tibble() %>%
  mutate(concept = NA, family = "TEA_all", n_cogsets = NA) %>%
  relocate(family, .before = n_trees) %>%
  relocate(c(concept, n_cogsets), .before = ub_DS) |>
  mutate(ub_DT = sum(filter(bounds_real_tb, family == "TEA")$n_cogsets * filter(bounds_real_tb, family == "TEA")$ub_DT))

bounds_real_tb <- bind_rows(TEA_summary, bounds_real_tb) |>
  mutate(
    inf_t_DT = compute_inf_t_DT(
    k=k, q=q, 
    depth = filter(tipages_summary, family == family)$depth, 
    s = filter(tipages_summary, family == family)$s
    ))

write_csv(bounds_real_tb, here("output/results/bounds_real_tb.csv"))


compute_inf_t_DT(3421,1.10,filter(tipages_summary, family == 'TEA')$depth, filter(tipages_summary, family == 'TEA')$s)

compute_inf_t_DT(57120,1.03,filter(tipages_summary, family == 'Bantu')$depth, filter(tipages_summary, family == 'Bantu')$s,
                 interval = c(-20, 0))



# dt_real_ages <- list.dirs(here("output/results"), full.names = TRUE, recursive = FALSE) %>%
#   map_df(function(x) {
#     ages <- list.files(list.dirs(here("output/results"), full.names = TRUE, recursive = FALSE), "tipages", full.names = TRUE) %>%
#       read_csv() %>%
#       select(age) %>%
#       tibble()
#     parameters <- read.csv(here("output/results/bounds_real_tb.csv"))
#   })
#
# # Find the extension of a file
# path_extension <- function(path) {
#   tolower(substr(path, nchar(path) - 3, nchar(path)))
# }
#
#
# # Liste des chemins des fichiers
# tipages_files <- list.files(list.dirs(here("output/results"), full.names = TRUE, recursive = FALSE), "tipages", full.names = TRUE)
#
# # Fonction qui lit les fichiers d'ages et sort les âges
# compute_ages <- function(x) {
#   if (path_extension(x) == ".csv") {
#     read_csv(x, col_types = cols()) %>%
#       select(age)
#   } else {
#     readRDS(x) %>%
#       select(age)
#   }
# }
#
# map_df(tipages_files, ~ compute_ages(.x)) # marche pas il ne fait pas de bind_rows()
