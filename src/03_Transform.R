library(here) 
library(tidyverse)
library(ggplot2)
library(dplyr)


# Compute the upper bound of the probability of inferring the true tree topology for one sens
compute_upperbound_DT <- function(k, q, t, depth, s) {
  k * sum(exp(-q  * (t - depth - s)))
}

# Compute the time threshold beyond which the upper bound of the probability
# of inferring the true tree topology falls below 1
compute_inf_t_DT <- function(k, q, depth, s, interval = c(0, 20), tol = 1e-6, maxiter = 10000) {
  uniroot(function(t) {
    compute_upperbound_DT(k, q, t, depth, s) - 1
  }, extendInt = "yes", interval = interval, tol = tol, maxiter = maxiter)$root
}

# Compute the upper bound of the probability of correctly inferring ancestral states
compute_upperbound_TauS <- function(pi0, pi1, q, t, depth) {
  sum(exp(-as.numeric(q) * (t - depth)))/(1 - max(as.numeric(pi0), as.numeric(pi1))) 
}

compute_upperbound_DS <- function(pi0, pi1, q, t, depth) {
  sum(exp(-as.numeric(q) * (t - depth)))
}

# Compute the time threshold beyond which the upper bound of the probability
# of correctly inferring ancestral states falls below 1
compute_inf_t_DS <- function(pi0, pi1, q, depth, interval = c(0, 1), tol = 1e-6, maxiter = 10000) {
  uniroot(function(t) {
    compute_upperbound_TauS(pi0, pi1, q , t, depth) - 0.10
  }, extendInt = "yes", interval = interval, tol = tol, maxiter = maxiter)$root 
}



compute_inf_t_DS2 <- function(pi0, pi1, q, mu, t, depth, interval = c(0, q), tol = 1e-6, maxiter = 10000) {
  compute_inf_t_DS2_varyq <- function(qq) {
    compute_upperbound_TauS(pi0, pi1, qq *mu , t, depth) - 1
  }
  uniroot(Vectorize(compute_inf_t_DS2_varyq), extendInt = "yes", interval = interval, tol = tol, maxiter = maxiter)$root * t /q
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
tracelog_summary <- read_csv(here("output/results/tracelog_summary.csv")) |>
  filter(!(family == "ST_by_sens" & concept == "the_name")) |>
  filter(!(family == "ST_by_sens" & concept == "four")) |>
  mutate(t_R = ifelse(family == "ST_by_sens", 11.0, t_R)) |> 
  mutate(t_R = ifelse(family == "TEA", 12.0, t_R)) 

# Bounds
bounds_real_tb_by_sens <- tracelog_summary |>
  mutate(n_cogsets = if_else(is.na(n_cogsets), k, n_cogsets)) |>
  rename(familyx = family) |>
  relocate(c(N, k), .after = n_trees) |>
  rowwise() |>
  mutate(
    ub_DS = compute_upperbound_DS(pi0, pi1, q * mu,
      t = t_R,
      depth = filter(tipages_summary, family == familyx)$depth
    ),
    #inf_t_DS = compute_inf_t_DS2(pi0, pi1, q, mu, t, depth = filter(tipages_summary, family == familyx)$depth, interval=c(q/2,q)),
    ub_DT = compute_upperbound_DT(k = n_cogsets, q = q * mu, t = t_R, depth = filter(tipages_summary, family == familyx)$depth, s = filter(tipages_summary, family == familyx)$s),
    inf_t_DT = compute_inf_t_DT(k = k, q = q * mu, depth = filter(tipages_summary, family == familyx)$depth, s = filter(tipages_summary, family == familyx)$s)
    ) |>
  #mutate(inf_t_DS = round(inf_t_DS, 3)) |>
  mutate(inf_t_DT = round(inf_t_DT, 3)) |>
  ungroup() |>
  rename(family = familyx) 
  #relocate(inf_t_DS, .after = ub_DT)

# ST par sens 
# Filtrage des données
foo_st <- filter(bounds_real_tb_by_sens, family == 'ST_by_sens')
foo_tea <- filter(bounds_real_tb_by_sens, family == 'TEA')

# Tri des valeurs en ordre décroissant
bar_st <- sort(foo_st$ub_DT, decreasing = TRUE)
bar_tea <- sort(foo_tea$ub_DT, decreasing = TRUE)

# Calcul des sommes cumulées
df_st <- data.frame(Index = seq_along(bar_st), Cumsum = cumsum(bar_st), Family = "Sino-Tibetan")
df_tea <- data.frame(Index = seq_along(bar_tea), Cumsum = cumsum(bar_tea), Family = "Transeurasian")

# Fusion des données
df_plot <- bind_rows(df_st, df_tea)


compute_inf_t_DS2(
  foo$pi0[10], foo$pi1[10], foo$q[10], foo$mu[10], foo$t_R[10],
  filter(tipages_summary, family == 'ST_by_sens')$depth)

# ST For each sens we compute t_inf 
tau <- matrix(nrow=nrow(foo), ncol=2)
tau[,1] <- foo$concept

for (k in 1:nrow(foo)){
  tau[k,2] <- compute_inf_t_DS2(
    foo$pi0[k], foo$pi1[k], foo$q[k], foo$mu[k], foo$t_R[k],
    filter(tipages_summary, family == 'ST_by_sens')$depth)
}

tau <- as_tibble(tau)



# TEA for each sens we compute t_inf 
tau_tea <- matrix(nrow=nrow(foo_tea), ncol=2)
tau_tea[,1] <- foo_tea$concept

for (k in 1:nrow(foo_tea)){
  tau_tea[k,2] <- compute_inf_t_DS2(
    foo_tea$pi0[k], foo_tea$pi1[k], foo_tea$q[k], foo_tea$mu[k], foo_tea$t_R[k],
    filter(tipages_summary, family == 'TEA')$depth)
}

tau_tea <- as_tibble(tau_tea)

write_csv(bounds_real_tb_by_sens, here("output/results/bounds_real_tb_by_sens.csv"))

# bound 
sum(foo_tea$ub_DT)


# TEA add line with the mean of parameter for all cognates
# check
TEA_summary <- bounds_real_tb_by_sens |>
  filter(family == "TEA") |>
  select(-concept, -family, -n_cogsets, -ub_DT, -ub_DS) |>
  colMeans(na.rm=T) |>
  t() |>
  as_tibble() |>
  mutate(concept = NA, family = "TEA_all", n_cogsets = NA) |>
  mutate(n_cogsets = if_else(is.na(n_cogsets), k, n_cogsets)) |>
  relocate(family, .before = n_trees) |>
  relocate(c(concept, n_cogsets), .before = t_R) |>
  #rowwise() |>
  mutate(
    ub_DT = sum(filter(bounds_real_tb_by_sens, family == "TEA")$ub_DT),
    ub_DS = sum(filter(bounds_real_tb_by_sens, family == 'TEA')$ub_DS),
    inf_t_DT = NA,
    inf_t_DS = NA)

ST_bysens_summary <- bounds_real_tb_by_sens |>
  filter(family == "ST_by_sens") |>
  select(-concept, -family, -n_cogsets, -ub_DT) |>
  colMeans(na.rm=T) |>
  t() |>
  as_tibble() |>
  mutate(concept = NA, family = "ST_bysens", n_cogsets = NA) |>
  mutate(n_cogsets = if_else(is.na(n_cogsets), k, n_cogsets)) |>
  relocate(family, .before = n_trees) |>
  relocate(c(concept, n_cogsets), .before = ub_DS) |>
  rowwise() |>
  mutate(
    ub_DT = sum(filter(bounds_real_tb_by_sens, family == "ST_by_sens")$ub_DT),
    ub_DS = sum(filter(bounds_real_tb_by_sens, family == 'ST_by_sens')$ub_DS),
    inf_t_DT = NA,
    inf_t_DS = NA)


#bounds_real_tb <- bind_rows(TEA_summary, bounds_real_tb_by_sens, ST_bysens_summary) |>
  #filter(family %in% c("TEA_all", "Bantu", "Bantu_subset", "Bantu_subset2", "IE", "ST","ST_bysens")) |> 
  #mutate(family = ifelse(family == "TEA_all", "TEA", family)) |>
  #rename(familyx = family) |>
  #rowwise() |>
  #mutate(
  #  inf_t_DS = compute_inf_t_DS(pi0, pi1, q * mu, depth = filter(tipages_summary, family == familyx)$depth),
  #  inf_t_DT = compute_inf_t_DT(k = k, q = q * mu, depth = filter(tipages_summary, family == familyx)$depth, s = filter(tipages_summary, family == familyx)$s)
  #  ) |>
  #rename(family = familyx) |>
  #relocate(inf_t_DT, .after = ub_DT) |>
  #select(-concept) |> 
  #relocate(ub_DS, .before= ub_DT ) |>
  #relocate(inf_t_DS, .after = ub_DT )


bounds_real_tb <- bounds_real_tb_by_sens |> 
  filter(family %in% c("Bantu", "IE", "ST")) |> 
  rename(familyx = family) |>
  #mutate(inf_t_DT = compute_inf_t_DT(k = k, q = q * mu, depth = filter(tipages_summary, family == familyx)$depth, 
  #                                   s = filter(tipages_summary, family == familyx)$s)) |>
  rename(family=familyx) |>
  select(-concept, -n_cogsets) 

write_csv(bounds_real_tb, here("output/results/bounds_real_tb.csv"))

# Bantu inf_t 
compute_inf_t_DT(k = 3859, q = 34.6 * 1, depth = filter(tipages_summary, family == 'Bantu')$depth, 
                 s = filter(tipages_summary, family == 'Bantu')$s)

compute_inf_t_DT(k = 4990, q = 130 * 1, depth = filter(tipages_summary, family == 'IE')$depth, 
                 s = filter(tipages_summary, family == 'IE')$s)

compute_inf_t_DT(k = 50, q = 9.31 * 1, depth = filter(tipages_summary, family == 'ST')$depth, 
                 s = filter(tipages_summary, family == 'ST')$s)


# Bounds by millenia 
t_values <- seq(0, 20, length.out = 100)
bounds_real_byt_tb <- bounds_real_tb |>
  filter(!family %in% c("ST_bysens", "TEA")) |>
  select(-inf_t_DS, -inf_t_DT, -n_cogsets)|>
  group_by(family) |>
  mutate(count = length(t_values)) |>
  uncount(count) |>
  mutate(t = t_values) |>
  mutate(family = ifelse(family == "ST_bysens", "ST_by_sens", family)) |>
  rename(familyx = family) |>
  rowwise() |>
  mutate(ub_DT = compute_upperbound_DT(k, q*mu, t, filter(tipages_summary, family == familyx)$depth ,filter(tipages_summary, family == familyx)$s),
         ub_DS = compute_upperbound_TauS(pi0, pi1, q*mu, t, filter(tipages_summary, family == familyx)$depth)
         ) |>
  rename(family = familyx)


# Bounds by millennial for rate heterogeneity
t_values <- seq(0, 20, length.out = 100)
bounds_real_byt_tb_by_sens <- bounds_real_tb_by_sens |>
  filter(family %in% c("ST_by_sens", "TEA")) |>
  group_by(concept, family) |>
  select(-inf_t_DS, -inf_t_DT) |> 
  mutate(count = length(t_values)) |>
  uncount(count) |>
  mutate(t = t_values) |>
  rename(familyx = family) |>
  rowwise() |>
  mutate(
    ub_DS = compute_upperbound_TauS(pi0, pi1, q*mu, t, 
                                  filter(tipages_summary, family == familyx)$depth),
    ub_DT = compute_upperbound_DT(n_cogsets, q*mu, t, 
                                  filter(tipages_summary, family == familyx)$depth, 
                                  filter(tipages_summary, family == familyx)$s)
  ) |>
  ungroup() |>
  group_by(familyx, t) |>  
  summarise(across(-ub_DT, mean), ub_DT = sum(ub_DT), .groups = "drop") |>
  select(-concept) |>
  rename(family = familyx)


bounds_real_byt_tb |> 
  bind_rows(bounds_real_byt_tb_by_sens) |>
  write_csv(here("output/results/bounds_real_byt_tb.csv"))

# On actualise bounds_real_tb ici maintenant qu'on à calculé inf_t_DS et inf_t_DT
bounds_real_tb = bounds_real_tb |> 
  mutate(inf_t_DS = if_else(family == 'TEA', 29.6, inf_t_DS)) |>
  mutate(inf_t_DS = if_else(family == 'ST_bysens', 26.6, inf_t_DS)) |>
  mutate(inf_t_DT = if_else(family == 'TEA', 62.6, inf_t_DT)) |>
  mutate(inf_t_DT = if_else(family == 'ST_bysens', 42, inf_t_DT)) 


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
