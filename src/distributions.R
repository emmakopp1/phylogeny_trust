library(here)
library(tidyverse)
library(openxlsx) 
library(tracerer)

burnin <- .2

ntipschars <- read_csv(here("output/results/ntipschars.csv"))


# Tip ages --------------------------------------------------------------------------------------------------------

tipages_bantu <- read_csv(here("output/results/bantu/bantu_ctmc-strict-bd_ages.csv.bz")) |>
  mutate(family = "Bantu") 
tipages_bantu_subsample <- read_csv(here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tipages.csv")) |>
  mutate(family = "Bantu_subset") 
tipages_bantu_subsample2 <- read_csv(here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tipages.csv")) |>
  mutate(family = "Bantu_subset2")
tipages_ie <- read_csv(here("output/results/ie/iecor_ctmc-strict-M1_tipages.csv.bz")) |>
  mutate(family = "IE")
tipages_st <- read_csv(here("output/results/st/st_ctmc-strict-fbd_tipages.csv.bz")) |>
  mutate(family = "ST")
tipages_st_by_sens <- read_csv(here("output/results/st_by_sens/st_ctmc-strict-fbd_by_sens_tipages.csv.bz")) |>
  mutate(family = "ST_by_sens")
tipages_tea <- read_csv(here("output/results/tea/tea_ctmc-strict-fbd-constrained_tipages.csv")) |>
  mutate(family = "TEA") |> 
  mutate(age = 0.1 * age) |> 
  mutate(depth = 0.1 * depth)

tipages_summary <- bind_rows(tipages_bantu, tipages_bantu_subsample, tipages_bantu_subsample2, tipages_ie, tipages_st, tipages_st_by_sens, tipages_tea) |>
  group_by(family) |>
  filter(tree > ceiling(max(tree) * burnin)) |>
  group_by(family, tip) |>
  summarise(age = median(age), depth = median(depth)) |>
  mutate(root_age = round(depth + age,2)) |>
  mutate(depth = round(depth,2)) 

write_csv(tipages_summary, here("output/results/tipages_summary.csv"))


# Trace logs -------------------------------------------------------------------------------------------------------

tracelog_bantu <- read_csv(here("output/results/bantu/bantu_ctmc-strict-bd_tracelog.csv")) |>
  mutate(family = "Bantu")
tracelog_bantu_subsample <- read_csv(here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tracelog.csv")) |>
  mutate(family = "Bantu_subset")
tracelog_bantu_subsample2 <- read_csv(here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tracelog.csv")) |>
  mutate(family = "Bantu_subset2")
tracelog_ie <- read_csv(here("output/results/ie/iecor_ctmc-strict-M1_tracelog.csv")) |>
  mutate(family = "IE")
tracelog_st <- read_csv(here("output/results/st/st_ctmc-strict-fbd_tracelog.csv")) |>
  mutate(family = "ST")
tracelog_st_by_sens <- read_csv(here("output/results/st_by_sens/st_ctmc-strict-fbd_by_sens_tracelog.csv")) |>
  mutate(family = "ST_by_sens")
tracelog_tea <- read_csv(here("output/results/tea/tea_ctmc-strict-fbd-constrained_tracelog.csv")) |>
  mutate(family = "TEA") |>
  rename(clockRate.c.clock = clockrate.c.clock)


n_cogids_tea <- here("data/real/tea_ctmc-strict-fbd-constrained/tea.nex") |>
  read_lines() |>
  str_subset("^charset") |>
  str_remove_all("^charset |;|\\?") |>
  enframe(name = NULL, value = "concept") |>
  separate(concept, into = c("concept", "sets"), sep = " = ") |>
  separate(sets, into = c("start", "end"), sep = "-") |>
  mutate(n_cogsets = as.integer(end) - as.integer(start) + 1) |>
  select(concept, n_cogsets)

n_cogids_st_by_sens <- here("data/real/st_ctmc-strict-fbd-by-sens/st.nex") |>
  read_lines() |>
  str_subset("charset") |>
  str_remove_all("    charset |;|\\?") |>
  enframe(name = NULL, value = "concept") |>
  separate(concept, into = c("concept", "sets"), sep = " = ") |>
  separate(sets, into = c("start", "end"), sep = "-") |>
  mutate(n_cogsets = as.integer(end) - as.integer(start) + 1) |>
  select(concept, n_cogsets)



# Tracelog summaries  
tracelog_tea_by_sens_summary <- tracelog_tea |>
  # Count rows
  add_tally(name = "n_trees") |>
  # Delete burn-in
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  select(Sample, family, n_trees, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
  pivot_longer(cols = matches("freqParameter|mutationRate"), names_to = "name", values_to = "value") |>
  mutate(
    name = str_replace(name, "^freqParameter.s.", "pi@"),
    name = str_replace(name, "^mutationRate.s.", "mu@"),
    name = str_replace(name, "\\.(\\d+)$", "@\\1")
  ) |>
  # Split the name column into three new columns: prefix, main_name, and suffix, using @ as the delimiter
  separate(name, into = c("prefix", "main_name", "suffix"), sep = "@", fill = "right")|>
  mutate(suffix = ifelse(is.na(suffix), "", suffix)) |>
  mutate(variable = case_when(
    prefix == "pi" & suffix != "" ~ paste0(prefix, suffix),
    prefix == "mu" ~ prefix,
    TRUE ~ prefix
  )) |>
  pivot_wider(names_from = variable, values_from = value)|>
  # Agréger les données sans modifier les colonnes non concernées
  group_by(Sample, family, n_trees, clockRate.c.clock, TreeHeight.t.tree,main_name) |>
  summarise(
    pi1 = sum(pi1, na.rm = TRUE),
    pi2 = sum(pi2, na.rm = TRUE),
    mu = sum(mu, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  mutate(concept = main_name) |>
  select(Sample,family, clockRate.c.clock, TreeHeight.t.tree, concept, pi1, pi2, mu) |>
  rename(t_R = TreeHeight.t.tree) |>
  rename(pi0 = pi1, pi1 = pi2) |>
  rename(clock_rate = clockRate.c.clock) |>
  left_join(n_cogids_tea, by = "concept") |>
  relocate(n_cogsets, .after = concept) |> 
  mutate(t_R = 0.1 * t_R)



# Sino-Tibetan family
tracelog_st_by_sens_summary = tracelog_st_by_sens |>
  # Count rows
  add_tally(name = "n_trees") |>
  # Delete burn-in
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  select(Sample, family, n_trees, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
  pivot_longer(cols = matches("freqParameter|mutationRate"), names_to = "name", values_to = "value") |>
  mutate(
    name = str_replace(name, "^freqParameter.s.", "pi@"),
    name = str_replace(name, "^mutationRate.s.", "mu@"),
    name = str_replace(name, "\\.(\\d+)$", "@\\1")
  ) |>
  # Split the name column into three new columns: prefix, main_name, and suffix, using @ as the delimiter
  separate(name, into = c("prefix", "main_name", "suffix"), sep = "@", fill = "right")|>
  mutate(suffix = ifelse(is.na(suffix), "", suffix)) |>
  mutate(variable = case_when(
    prefix == "pi" & suffix != "" ~ paste0(prefix, suffix),
    prefix == "mu" ~ prefix,
    TRUE ~ prefix
  )) |>
  pivot_wider(names_from = variable, values_from = value)|>
  group_by(Sample, family, n_trees, clockRate.c.clock, TreeHeight.t.tree, main_name) |>
  summarise(
    pi1 = sum(pi1, na.rm = TRUE),
    pi2 = sum(pi2, na.rm = TRUE),
    mu = sum(mu, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  mutate(concept = main_name) |>
  select(Sample, family, n_trees, clockRate.c.clock, TreeHeight.t.tree, concept, pi1, pi2, mu) |>
  rename(t_R = TreeHeight.t.tree) |>
  rename(pi0 = pi1, pi1 = pi2) |>
  rename(clock_rate = clockRate.c.clock) |>
  left_join(n_cogids_st_by_sens, by = "concept") |>
  relocate(n_cogsets, .after = concept)


tracelog_summary <- list(tracelog_bantu, tracelog_bantu_subsample, tracelog_bantu_subsample2, tracelog_ie, tracelog_st) |>
  map(~ .x |>
        add_tally(name = "n_trees") |>
        filter(Sample > ceiling(max(Sample) * burnin)) |>
        select(Sample, family, n_trees, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
        rename(t_R = TreeHeight.t.tree) |>
        rename(clock_rate = clockRate.c.clock) |>
        rename_with(~ str_replace(.x, "freqParameter.+(?=\\d$)", "pi")) |>
        rename_with(~ str_replace(.x, "mutationRate\\.s\\.(.*)", "mu")) |>
        rename(pi0 = pi1, pi1 = pi2)) |> 
  bind_rows(tracelog_tea_by_sens_summary) |>
  bind_rows(tracelog_st_by_sens_summary) |> 
  mutate(q = pi0 + pi1) |> 
  relocate(q, .after = pi1) |>
  relocate(mu, .before = q) |> 
  left_join(ntipschars)

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
compute_upperbound_DS <- function(pi0, pi1, q, t, depth) {
  max(as.numeric(pi0), as.numeric(pi1)) + sum(exp(-as.numeric(q) * (t - depth)))
}

# Compute the time threshold beyond which the upper bound of the probability
# of correctly inferring ancestral states falls below 1
compute_inf_t_DS <- function(pi0, pi1, q, depth, interval = c(0, 20), tol = 1e-6, maxiter = 10000) {
  uniroot(function(t) {
    compute_upperbound_DS(pi0, pi1, q, t, depth) - 1
  }, extendInt = "yes", interval = interval, tol = tol, maxiter = maxiter)$root
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


tracelog_summary = tracelog_summary |> 
  filter(!(family == "ST_by_sens" & concept == "the_name")) |>
  filter(!(family == "ST_by_sens" & concept == "four"))


# Bounds
bounds_real_tb_by_sens_distribution <- tracelog_summary |>
  filter(family=="IE") |> 
  mutate(n_cogsets = if_else(is.na(n_cogsets), k, n_cogsets)) |>
  relocate(c(N, k), .after = n_trees) |>
  rowwise() |>
  mutate(
    ub_DS = compute_upperbound_DS(pi0, pi1, q*mu,
                                  t = t_R,
                                  depth = filter(tipages_summary, family == "IE")$depth
    ),
    ub_DT = compute_upperbound_DT(k = n_cogsets, 
                                  q = q * mu, 
                                  t = t_R, 
                                  depth = filter(tipages_summary, family == "IE")$depth, 
                                  s = filter(tipages_summary, family == "IE")$s)
  ) |>
  ungroup() 


d = sapply(bounds_real_tb_by_sens_distribution$ub_DT, as.numeric)
a = sapply(bounds_real_tb_by_sens_distribution$ub_DS, as.numeric)

hist(d, breaks = 70, xlim = c(min(d), max(d)))
hist(a, breaks = 70, xlim = c(min(a), max(a)))


quantile(d, 0.05)
quantile(d, 0.95)


quantile(a, 0.05)
quantile(a, 0.95)









