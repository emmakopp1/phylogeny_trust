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
tipages_kd <- read_csv(here("output/results/kd/kd_ctmc-strict-bd_tipages.csv.bz")) |>
  mutate(family = "KD")

tipages_summary <- bind_rows(tipages_bantu, tipages_bantu_subsample, tipages_bantu_subsample2, tipages_ie, tipages_st, tipages_st_by_sens, tipages_tea, tipages_kd) |>
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
tracelog_kd <- read_csv(here("output/results/kd/kd_ctmc-strict-bd_tracelog.csv")) |>
  mutate(family = "KD")


n_cogids_tea <- here("data/real/tea_ctmc-strict-fbd-constrained/tea.nex") |>
  read_lines() |>
  str_subset("^charset") |>
  str_remove_all("^charset |;|\\?") |>
  enframe(name = NULL, value = "concept") |>
  separate(concept, into = c("concept", "sets"), sep = " = ") |>
  separate(sets, into = c("start", "end"), sep = "-") |>
  mutate(n_cogsets = as.integer(end) - as.integer(start) + 1) |>
  select(concept, n_cogsets)

n_cogids_st_by_sens <- here("data/real/st_ctmc-strict-fbd-ht/st.nex") |>
  read_lines() |>
  str_subset("charset") |>
  str_remove_all("    charset |;|\\?") |>
  enframe(name = NULL, value = "concept") |>
  separate(concept, into = c("concept", "sets"), sep = " = ") |>
  separate(sets, into = c("start", "end"), sep = "-") |>
  mutate(n_cogsets = as.integer(end) - as.integer(start) + 1) |>
  select(concept, n_cogsets)

n_cogids_kd <- here("data/real/kd_ctmc-strict-bd-ht/kd.nex") |>
  read_lines() |>
  str_subset("charset") |>
  str_remove_all("    charset |;|\\?") |>
  enframe(name = NULL, value = "concept") |>
  separate(concept, into = c("concept", "sets"), sep = " = ") |>
  separate(sets, into = c("start", "end"), sep = "-") |>
  mutate(n_cogsets = as.integer(end) - as.integer(start) + 1) |>
  select(concept, n_cogsets) |> 
  mutate(concept = str_remove(concept, "^charset "))

# Tracelog summaries  
# TEA 
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
  pivot_wider(names_from = variable, values_from = value) |>
  # Agréger les données sans modifier les colonnes non concernées
  group_by(Sample, n_trees, family, clockRate.c.clock, TreeHeight.t.tree, main_name) |>
  summarise(
    pi1 = sum(pi1, na.rm = TRUE),
    pi2 = sum(pi2, na.rm = TRUE),
    mu = sum(mu, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  mutate(concept = main_name) |>
  select(family, n_trees, clockRate.c.clock, TreeHeight.t.tree, concept, pi1, pi2, mu) |>
  rename(t_R = TreeHeight.t.tree) |>
  rename(pi0 = pi1, pi1 = pi2) |>
  rename(clock_rate = clockRate.c.clock) |>
  group_by(family, concept) |>
  summarise(across(c(t_R, pi0, pi1, mu, clock_rate), ~ median(.x))) |>
  ungroup() |>
  left_join(n_cogids_tea, by = "concept") |>
  relocate(n_cogsets, .after = concept) |> 
  mutate(t_R = 0.1 * t_R)



# ST_by_sens
tracelog_st_by_sens_summary = tracelog_st_by_sens |>
  #select(-st_ess$parameter)|>
  # Ajouter un compteur de lignes si nécessaire
  add_tally(name = "n_trees") |>
  # Filtrer les échantillons en fonction de burnin
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  # Sélectionner les colonnes d'intérêt
  select(family, n_trees, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
  # Transformer en format long
  pivot_longer(cols = matches("freqParameter|mutationRate"), names_to = "name", values_to = "value") |>
  mutate(
    name = str_replace(name, "^freqParameter.s.", "pi@"),
    name = str_replace(name, "^mutationRate.s.", "mu@"),
    name = str_replace(name, "\\.(\\d+)$", "@\\1")
  ) |>
  # Séparer la colonne 'name' en 'prefix', 'main_name', 'suffix' en utilisant '@' comme séparateur
  separate(name, into = c("prefix", "main_name", "suffix"), sep = "@", fill = "right")|>
  # Remplacer NA dans 'suffix' par une chaîne vide
  mutate(suffix = ifelse(is.na(suffix), "", suffix)) |>
  # Créer les noms de variables pour pivot_wider
  mutate(variable = case_when(
    prefix == "pi" & suffix != "" ~ paste0(prefix, suffix),
    prefix == "mu" ~ prefix,
    TRUE ~ prefix
  )) |>
  # Transformer en format large
  pivot_wider(names_from = variable, values_from = value)|>
  # Agréger les données sans modifier les colonnes non concernées
  group_by(family, n_trees, clockRate.c.clock, TreeHeight.t.tree, main_name) |>
  summarise(
    pi1 = sum(pi1, na.rm = TRUE),
    pi2 = sum(pi2, na.rm = TRUE),
    mu = sum(mu, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  mutate(concept = main_name) |>
  select(family, n_trees, clockRate.c.clock, TreeHeight.t.tree, concept, pi1, pi2, mu) |>
  rename(t_R = TreeHeight.t.tree) |>
  rename(pi0 = pi1, pi1 = pi2) |>
  rename(clock_rate = clockRate.c.clock) |>
  group_by(family, n_trees, concept) |>
  summarise(across(c(t_R, pi0, pi1, mu, clock_rate), ~ median(.x))) |>
  ungroup() |>
  left_join(n_cogids_st_by_sens, by = "concept") |>
  relocate(n_cogsets, .after = concept)


# KD 
tracelog_kd_summary <- tracelog_kd |>
  # Count rows
  add_tally(name = "n_trees") |>
  # Delete burn-in
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  select(Sample, family, n_trees, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
  rename_with(~ str_replace(.x, "^freqParameter.*\\.(\\d+)$", "pi\\1"), starts_with("freqParameter")) |> 
  rename_with(~ str_replace(.x, "^mutationRate.*", "mu"), starts_with("mutationRate")) |> 
  group_by(Sample, family, clockRate.c.clock, TreeHeight.t.tree) |>
  summarise(
    pi1 = sum(pi1, na.rm = TRUE),
    pi2 = sum(pi2, na.rm = TRUE),
    mu = sum(mu, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  select(family, n_trees, clockRate.c.clock, TreeHeight.t.tree, pi1, pi2, mu) |>
  rename(t_R = TreeHeight.t.tree) |>
  rename(pi0 = pi1, pi1 = pi2) |>
  rename(clock_rate = clockRate.c.clock) |>
  group_by(family, concept) |>
  summarise(across(c(t_R, pi0, pi1, mu, clock_rate), ~ median(.x))) |>
  ungroup() |>
  left_join(n_cogids_tea, by = "concept") |>
  relocate(n_cogsets, .after = concept) |> 
  mutate(t_R = 0.1 * t_R)


tracelog_summary <- list(tracelog_bantu, tracelog_bantu_subsample, tracelog_bantu_subsample2, tracelog_ie, tracelog_st, tracelog_kd) |>
  purrr::map(~ .x |>
    add_tally(name = "n_trees")|>
    filter(Sample > ceiling(max(Sample) * burnin)) |>
    select(family, n_trees, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
    summarise(family = unique(family), across(-family, ~ median(.x))) |>
    rename(t_R = TreeHeight.t.tree) |>
    rename(clock_rate = clockRate.c.clock) |>
    rename_with(~ str_replace(.x, "freqParameter.+(?=\\d$)", "pi")) |>
    rename_with(~ str_replace(.x, "mutationRate\\.s\\.(.*)", "mu")) |>
    rename(pi0 = pi1, pi1 = pi2)) |> 
  bind_rows(tracelog_tea_by_sens_summary) |>
  bind_rows(tracelog_st_by_sens_summary) |> 
  mutate(q = (pi0 + pi1)/(2 * pi0 * pi1)) |> 
  relocate(q, .after = pi1) |>
  relocate(mu, .before = q) |> 
  mutate(t_R = ifelse(family == "KD", 1e-3 * t_R, t_R))
  left_join(ntipschars)


write_csv(tracelog_summary, here("output/results/tracelog_summary.csv"))


#  ESS  --------------------------------------------------------------------------- ---------------------------------
ess <- list(tracelog_bantu, tracelog_bantu_subsample, tracelog_bantu_subsample2, tracelog_ie, tracelog_st, tracelog_kd) |>
  map_df(function(x) {
    # Extraire la colonne family
    family_name <- unique(x$family)
    
    # Effectuer les calculs sur le data frame sans family
    x |>
      add_tally(name = "n_trees") |>
      select(Sample, starts_with("freqParameter"), clockRate.c.clock, TreeHeight.t.tree, starts_with("mutationRate")) |>
      filter(Sample > ceiling(max(Sample) * burnin)) |>
      rowid_to_column() |>
      mutate(burnin = rowid <= max(rowid) * burnin) |>
      filter(burnin == FALSE) |>
      select(-rowid, -burnin) |>
      as.data.frame() |>
      calc_esses(sample_interval = max(x$Sample) / (nrow(x) - 1)) |>
      as_tibble() |>
      mutate(family = family_name) |>
      rename(t_R = TreeHeight.t.tree) |>
      rename(clock_rate = clockRate.c.clock) |>
      rename_with(~ str_replace(.x, "freqParameter.+(?=\\d$)", "pi")) |>
      rename_with(~ str_replace(.x, "mutationRate\\.s\\.(.*)", "mu")) |>
      rename(pi0 = pi1, pi1 = pi2)
  }) |> 
  bind_rows() |> 
  relocate(family, .before = pi0)


# Heterogeneity rate 
st_log <- tracelog_st_by_sens |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "st-by-sens") |> 
  select(-family) 


st_ess <- st_log |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  as.data.frame() |>
  calc_esses(sample_interval = max(st_log$Sample) / (nrow(st_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS") |> 
  mutate(family= "ST_by_sens")


tea_log <- tracelog_tea |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "tea-by-sens") |> 
  select(-family) 

tea_ess <- tea_log |> 
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  as.data.frame() |>
  calc_esses(sample_interval = max(tea_log$Sample) / (nrow(tea_log) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS") |>
  mutate(family= "TEA")

list(tea_ess, st_ess) |>
  bind_rows() |>
  relocate(family, .before= parameter) |>
  write_csv(here("output/results/ess_heterogene.csv"))

write_csv(ess,here("output/results/ess.csv")) 
# # Get the number of taxa and traits from a nexus file
# get_nexus_parameters <- function(file) {
#   group_name <- str_extract(basename(file), "^[^.]+")
#   phydt <- ReadAsPhyDat(file)
#   tibble(N = length(attributes(phydt)$names), k = length(attributes(phydt)$index), family = group_name)
# }
#
# # Get the values of pi0, pi1, the number of generated trees, and compute q
# get_tracerlog_parameters <- function(file, burnin = 0.2) {
#   beast_log_full <- read.csv(file)
#   beast_log <- remove_burn_ins(beast_log_full, burn_in_fraction = burnin)
#   beast_log |>
#     select(starts_with("freqParameter"), TreeHeight.t.tree) |>
#     summarise(across(everything(), ~ median(.x))) |>
#     select(matches("\\d$"), TreeHeight.t.tree) %>%
#     pivot_longer(cols = -c(TreeHeight.t.tree), names_to = "Column", values_to = "Value") %>%
#     mutate(Group = str_extract(Column, "\\d$"),TreeHeight = TreeHeight.t.tree) %>%
#     group_by(Group) %>%
#     summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
#     pivot_wider(names_from = Group, values_from = Mean_Value) %>%
#     rename(pi0 = "1", pi1 = "2") %>%
#     mutate(q = 1 / (pi0^2 + pi1^2)) %>%
#     mutate(t_R = median(beast_log$TreeHeight.t.tree)) %>%
#     relocate(q, .after = pi1)
# }
#
#
# # Sino-tibetain
# directories <- list.dirs(here(), full.names = TRUE, recursive = FALSE)
#
# dt_real <- map2(
#   list.files(directories, pattern = "tracelog.*\\.csv", full.names = TRUE, recursive = TRUE),
#   list.files(here(directories, 'real'), pattern = "\\.nex", full.names = TRUE, recursive = TRUE),
#   ~ bind_cols(get_tracerlog_parameters(.x), get_nexus_parameters(.y))
# ) %>%
#   bind_rows() %>%
#   relocate(c(family,N,k), .before = pi0) %>%
#   mutate(family = case_when(
#     str_detect(family, "^bantusubsample$") ~ "Bantu subset",
#     str_detect(family, "^bantusubsample2$") ~ "Bantu subset 2",
#     str_detect(family, "^bantu") ~ "Bantu",
#     str_detect(family, "^ie") ~ "Indo-European",
#     str_detect(family, "^st$") ~ "Sino-Tibetan",
#     str_detect(family, "^st.+tibetanPrior$") ~ "Sino-Tibetan prior",
#     str_detect(family, "^tea") ~ "Trans-Eurasian"
#   ))
#
# tipages_files = list.files(
#   list.dirs(here("output/results"), full.names = TRUE, recursive = FALSE),
#   "tipages", full.names = TRUE)
#
#
#
# write_csv(dt_real, here("output/results/bounds_real_tb.csv"))
