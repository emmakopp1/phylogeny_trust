library(here)
library(tidyverse)


# Compute the upper bound of the probability of inferring the true tree topology
compute_upperbound_DT <- function(k, N, q, t, s = 0) {
  k * N * exp(-q * (t - s))
}

# Compute the upper bound of the probability of inferring the true tree topology
compute_upperbound_DT <- function(k, N, q, t, depth, s) {
  k * N * exp(-q * (t - depth - s))
}

# Compute the upper bound of the probability of correctly inferring ancestral states version2
compute_upperbound_DS <- function(pi0, pi1, q, t, depth) {
  max(as.numeric(pi0), as.numeric(pi1)) + sum(exp(-as.numeric(q) *(t - depth)))
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
compute_inf_t_DS <- function(pi0, pi1, q, depth, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DS(pi0,pi1,q,t,depth) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}


tipages_summary <- read_csv(here("output/results/tipages_summary.csv"))

tracelog_summary <- read_csv(here("output/results/tracelog_summary.csv"))

bounds_real_tb <- tracelog_summary |>
  relocate(c(N, k), .after = n_trees) |> 
  rowwise() |>
  mutate(
    ub_DS = compute_upperbound_DS(pi0, pi1, q, t =filter(tipages_summary, family == family)$root_age, 
                                  depth = filter(tipages_summary, family == family)$depth),
    inf_t_DS = compute_inf_t_DS(pi0, pi1, q, depth = filter(tipages_summary, family == family)$depth)
  )

TEA_summary <- bounds_real_tb %>% 
  filter(family == "TEA") %>%
  select(-concept, -family, -n_cogsets) %>%
  colMeans() %>%
  as.data.frame() %>%
  t() %>%
  as_tibble()%>%
  mutate(concept = NA, family = "TEA_all", n_cogsets = NA) %>%
  relocate(family, .before=n_trees) %>%
  relocate(c(concept,n_cogsets), .before = ub_DS)

bounds_real_tb <- bind_rows(bounds_real_tb,TEA_summary)

write_csv(bounds_real_tb, here("output/results/bounds_real_tb.csv"))



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
