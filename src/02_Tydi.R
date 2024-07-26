library(here)
library(tidyverse)
library(openxlsx) 

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

# Cutting points -------------------------------------------------------------------------------------------------

tipages_names <- tipages_summary |> 
  select(family, tip, root_age) |> 
  mutate(s=NA)

write.xlsx(tipages_names, here("output/results/tipages_time_prior.xlsx"))

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
  mutate(family = "TEA")

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

tracelog_tea_summary <- tracelog_tea |>
  add_tally(name = "n_trees") |>
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  select(family, n_trees, starts_with("freqParameter"), TreeHeight.t.tree, starts_with("mutationRate")) |>
  pivot_longer(cols = matches("freqParameter|mutationRate")) |>
  mutate(name = str_remove(name, "freqParameter\\.s\\.")) |>
  mutate(name = str_replace(name, "mutationRate", "mu.")) |>
  mutate(name = str_remove(name, "\\.s\\.")) |>
  mutate(name = str_replace_all(name, "(mu.)(.*)", "\\2.mu")) |>
  separate(name, into = c("concept", "pi"), sep = "\\.(?=[12]$|mu$)", remove = FALSE) |>
  mutate(pi = if_else(str_detect(name, "mu$"), "mu", pi)) |>
  select(-name) |>
  mutate(pi = if_else(pi == "mu", "mu", paste0("pi", as.integer(str_remove(pi, "pi")) - 1))) |>
  mutate(concept = str_remove(concept, "\\.$")) |>
  mutate(concept = str_replace(concept, "^fly$", "fly_noun")) |>
  rename(t_R = TreeHeight.t.tree) |>
  mutate(t_R = t_R * 0.1) |>
  pivot_wider(names_from = pi, values_from = value) |>
  group_by(family, n_trees, concept) |>
  summarise(across(c(t_R, pi0, pi1, mu), ~ median(.x))) |>
  ungroup() |>
  left_join(n_cogids_tea, by = "concept") |>
  relocate(n_cogsets, .after = concept)



tracelog_st_by_sens_summary <- tracelog_st_by_sens |>
  add_tally(name = "n_trees") |>
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  select(family, n_trees, starts_with("freqParameter"), TreeHeight.t.tree, starts_with("mutationRate")) |>
  pivot_longer(cols = matches("freqParameter|mutationRate")) |>
  mutate(name = str_remove(name, "freqParameter\\.s\\.")) |>
  mutate(name = str_replace(name, "mutationRate", "mu.")) |>
  mutate(name = str_remove(name, "\\.s\\.")) |>
  mutate(name = str_replace_all(name, "(mu.)(.*)", "\\2.mu")) |>
  separate(name, into = c("concept", "pi"), sep = "\\.(?=[12]$|mu$)", remove = FALSE) |>
  mutate(pi = if_else(str_detect(name, "mu$"), "mu", pi)) |>
  select(-name) |>
  mutate(pi = if_else(pi == "mu", "mu", paste0("pi", as.integer(str_remove(pi, "pi")) - 1))) |>
  mutate(concept = str_remove(concept, "\\.$")) |>
  rename(t_R = TreeHeight.t.tree) |>
  pivot_wider(names_from = pi, values_from = value) |>
  group_by(family, n_trees, concept) |>
  summarise(across(c(t_R, pi0, pi1, mu), ~ median(.x))) |>
  ungroup() |>
  left_join(n_cogids_st_by_sens, by = "concept") |>
  filter(!(concept %in% c("the_.1", "the_.2","the_mud"))) |> # To automatize
  relocate(n_cogsets, .after = concept)


# Check ESS < 200 
## Multiple rates
test <- tracelog_st_by_sens |>
  rowid_to_column() |>
  mutate(burnin = rowid <= max(rowid) * burnin) |>
  mutate(data = "st-by-sens") |> 
  select(-family)

test_final <- test |>
  filter(burnin == FALSE) |>
  select(-rowid, -burnin, -data) |>
  calc_esses(sample_interval = max(test$Sample) / (nrow(test) - 1)) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = "parameter", values_to = "ESS")
filter(kd_lgs_bcov_byconcept_ess, ESS < 200)



tracelog_summary <- list(tracelog_bantu, tracelog_bantu_subsample, tracelog_bantu_subsample2, tracelog_ie, tracelog_st) |>
  map(~ .x |>
    add_tally(name = "n_trees") |>
    filter(Sample > ceiling(max(Sample) * burnin)) |>
    select(family, n_trees, starts_with("freqParameter"), TreeHeight.t.tree, starts_with("mutationRate")) |>
    summarise(family = unique(family), across(-family, ~ median(.x))) |>
    rename(t_R = TreeHeight.t.tree) |>
    rename_with(~ str_replace(.x, "freqParameter.+(?=\\d$)", "pi"))|>
    rename_with(~ str_replace(.x, "mutationRate\\.s\\.(.*)", "mu")) |>
    rename(pi0 = pi1, pi1 = pi2)) |>
  bind_rows(tracelog_tea_summary) |>
  bind_rows(tracelog_st_by_sens_summary) |> 
  mutate(q = 1 / (pi0^2 + pi1^2)) |>
  relocate(q, .after = pi1) |>
  relocate(mu, .before = q) |> 
  left_join(ntipschars)


write_csv(tracelog_summary, here("output/results/tracelog_summary.csv"))


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
