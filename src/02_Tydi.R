library(here)
library(tidyverse)

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
tipages_tea <- read_csv(here("output/results/tea/tea_ctmc-strict-fbd-constrained_tipages.csv")) |>
  mutate(family = "TEA")

tipages_summary <- bind_rows(tipages_bantu, tipages_bantu_subsample, tipages_bantu_subsample2, tipages_ie, tipages_st, tipages_tea) |>
  group_by(family) |>
  filter(tree > ceiling(max(tree) * burnin)) |>
  group_by(family, tip) |>
  summarise(age = median(age))

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
tracelog_tea <- read_csv(here("output/results/tea/tea_ctmc-strict-fbd-constrained_tracelog.csv")) |>
  mutate(family = "TEA")

n_cogids_tea <- here("data/real/tea_ctmc-strict-fbd-constrained/tea.nex") |>
  read_lines() |>
  str_subset("^charset") |>
  str_remove_all("^charset |;|\\?") |>
  enframe(name = NULL, value = "concept") |>
  separate(concept, into = c("concept", "sets"), sep = " = ") |>
  separate(sets, into = c("start", "end"), sep = "-") |>
  mutate(n_cogids = as.integer(end) - as.integer(start) + 1) |>
  select(concept, n_cogids)

tracelog_tea_summary <- tracelog_tea |>
  add_tally(name = "n_trees") |>
  filter(Sample > ceiling(max(Sample) * burnin)) |>
  select(family, n_trees, starts_with("freqParameter"), TreeHeight.t.tree) |>
  pivot_longer(starts_with("freqParameter")) |>
  mutate(name = str_remove(name, "freqParameter\\.s\\.")) |>
  separate(name, into = c("concept", "pi"), sep = "\\.(?=[12]$)") |>
  mutate(pi = paste0("pi", as.integer(pi) - 1)) |>
  mutate(concept = str_remove(concept, "\\.$")) |>
  mutate(concept = str_replace(concept, "^fly$", "fly_noun")) |>
  rename(t_R = TreeHeight.t.tree) |>
  pivot_wider(names_from = pi, values_from = value) |>
  group_by(family, n_trees, concept) |>
  summarise(across(c(t_R, pi0, pi1), ~ median(.x))) |>
  ungroup() |>
  left_join(n_cogids_tea) |>
  relocate(n_cogids, .after = concept)

tracelog_summary <- list(tracelog_bantu, tracelog_bantu_subsample, tracelog_bantu_subsample2, tracelog_st) |>
  map(~ .x |>
    add_tally(name = "n_trees") |>
    filter(Sample > ceiling(max(Sample) * burnin)) |>
    select(family, n_trees, starts_with("freqParameter"), TreeHeight.t.tree) |>
    summarise(family = unique(family), across(-family, ~ median(.x))) |>
    rename(t_R = TreeHeight.t.tree) |>
    rename_with(~ str_replace(.x, "freqParameter.+(?=\\d$)", "pi")) |>
    rename(pi0 = pi1, pi1 = pi2)) |>
  bind_rows(tracelog_tea_summary) |>
  mutate(q = 1 / (pi0^2 + pi1^2)) |>
  relocate(q, .after = pi1) |>
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
