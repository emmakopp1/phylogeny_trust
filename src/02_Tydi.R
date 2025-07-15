library(here)
library(tidyverse)
library(openxlsx) 
library(tracerer)

burnin <- .1

ntipschars <- read_csv(here("output/results/ntipschars.csv"))

# Tip ages --------------------------------------------------------------------------------------------------------

tipages_ie <- read_csv(here("output/results/ie/iecor_ctmc-strict-M1_tipages.csv.bz")) |>
  mutate(family = "IE")
tipages_st <- read_csv(here("output/results/st/st_ctmc-strict-fbd_tipages.csv.bz")) |>
  mutate(family = "ST")

tipages_summary <- bind_rows(tipages_ie, tipages_st) |>
  group_by(family) |>
  filter(tree > ceiling(max(tree) * burnin)) |> 
  group_by(family, tip) |>
  summarise(age = median(age), depth = median(depth)) |>
  mutate(root_age = round(depth + age,2)) |>
  mutate(depth = round(depth,2)) 

write_csv(tipages_summary, here("output/results/tipages_summary.csv"))


# Trace logs -------------------------------------------------------------------------------------------------------

tracelog_ie <- read_csv(here("output/results/ie/iecor_ctmc-strict-M1_tracelog.csv")) |>
  mutate(family = "IE")
tracelog_st <- read_csv(here("output/results/st/st_ctmc-strict-fbd_tracelog.csv")) |>
  mutate(family = "ST")


tracelog_summary <- list(tracelog_ie, tracelog_st) |>
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
  bind_rows() |>
  mutate(q = (pi0 + pi1)) |> 
  relocate(q, .after = pi1) |>
  relocate(mu, .before = q) |> 
  left_join(ntipschars)


write_csv(tracelog_summary, here("output/results/tracelog_summary.csv"))
