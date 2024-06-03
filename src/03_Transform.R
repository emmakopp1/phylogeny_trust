library(here)
library(tidyverse)


# Compute the upper bound of the probability of inferring the true tree topology
compute_upperbound_DT <- function(k, N, q, t, s = 0) {
  k * N * exp(-q * (t - s))
}

# Compute the upper bound of the probability of correctly inferring ancestral states version2
compute_upperbound_DS <- function(pi0, pi1, q, ages) {
  max(pi0, pi1) + sum(exp(-as.numeric(q) * ages))
}

compute_upperbound_DS2 <- function(pi0, pi1, q, t, N) {
  max(pi0, pi1) + N*exp(-q * t)
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
compute_inf_t_DS <- function(pi0, pi1, q, N, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DS2(pi0,pi1,q,t,N) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}


tipages_summary <- read_csv(here("output/results/tipages_summary.csv"))

tracelog_summary <- read_csv(here("output/results/tracelog_summary.csv"))

bounds_real_tb <- tracelog_summary |>
  relocate(c(N, k), .after = n_trees) |> 
  rowwise() |>
  mutate(
    #ub_DT = compute_upperbound_DT(k, N, q, t_R),
    #ub_DS = compute_upperbound_DS(pi0, pi1, q, ages = filter(tipages_summary, family == family)$age, N),
    ub_DS = compute_upperbound_DS(pi0, pi1, q, ages = filter(tipages_summary, family == family)$age),
    #inf_t_DT = compute_inf_t_DT(k, N, q, t_R),
    inf_t_DS = compute_inf_t_DS(pi0, pi1, q, N)
  )


write_csv(bounds_real_tb, here("output/results/bounds_real_tb.csv"))

tt

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
