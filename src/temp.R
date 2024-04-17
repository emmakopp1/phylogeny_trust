library(here)
library(tidyverse)
library(jsonlite)

library(tracerer)

logfile <- here("data/real/sino-tibet-ctmc-strict-bd-fossilsRemoved/sino-tibetan-ctmc-strict-bd.log")
# logfile <- here("data/real/IECoR-ctmc-strict-fbd/IECoR2-chr_1695819208636.log")
beast_log_full <- parse_beast_tracelog_file(logfile)
beast_log <- remove_burn_ins(beast_log_full, burn_in_fraction = 0.2)
beast_log |> 
  as_tibble() |> 
  select(starts_with("freqParameter"), TreeHeight.t.tree) |> 
  summarise(across(everything(), ~ mean(.x))) |> 
  rename_all(str_replace, pattern = "freq.+(\\d)", replacement = "pi\\1") |> 
  rename(pi0 = pi1, pi1 = pi2, t_R = TreeHeight.t.tree) |> 
  mutate(nTrees = max(beast_log$Sample)) |> 
  mutate(family = "Sino-Tibetan") |> 
  relocate(family, .before = pi0)


max(beast_log_full$Sample)

cfg <- here("src/compute_bounds/config.json")
fromJSON(cfg) |> 
  as_tibble() |> 
  mutate(c(""))
  pivot_longer(everything()) |> 
  pivot_wider(-name)
  # t() |> 
  # enframe()
  # rename(family = 1, path = 2, pi0 = 3, pi1 = 4, k = 5, path_cognates = 6, t_conv = 7)
