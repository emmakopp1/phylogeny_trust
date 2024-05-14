library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)
library(TreeTools)


# Get the number of taxa and traits from a nexus file
get_nexus_parameters <- function(file) {
  group_name <- str_extract(basename(file), "^[^.]+")
  phydt <- ReadAsPhyDat(file)
  tibble(N = length(attributes(phydt)$names), k = length(attributes(phydt)$index), Group = group_name)
}

# Get the values of pi0, pi1, the number of generated trees, and compute q
get_tracerlog_parameters <- function(file, burnin = 0.2) {
  beast_log_full <- read.csv(file)
  beast_log <- remove_burn_ins(beast_log_full, burn_in_fraction = burnin)
  beast_log |>
    select(starts_with("freqParameter"), TreeHeight.t.tree) |>
    summarise(across(everything(), ~ median(.x))) |>
    select(matches("\\d$"), TreeHeight.t.tree) %>%
    pivot_longer(cols = -c(TreeHeight.t.tree), names_to = "Column", values_to = "Value") %>%
    mutate(Group = str_extract(Column, "\\d$"),TreeHeight = TreeHeight.t.tree) %>%
    group_by(Group) %>%
    summarise(Mean_Value = mean(Value, na.rm = TRUE)) %>%
    pivot_wider(names_from = Group, values_from = Mean_Value) %>%
    rename(pi0 = "1", pi1 = "2") %>%
    mutate(q = 1 / (pi0^2 + pi1^2)) %>%
    mutate(t_R = median(beast_log$TreeHeight.t.tree)) %>%
    relocate(q, .after = pi1)
}


# Sino-tibetain
directories <- list.dirs(here(), full.names = TRUE, recursive = FALSE)

dt_real <- map2(
  list.files(directories, pattern = "tracelog.*\\.csv", full.names = TRUE, recursive = TRUE),
  list.files(here(directories, 'real'), pattern = "\\.nex", full.names = TRUE, recursive = TRUE),
  ~ bind_cols(get_tracerlog_parameters(.x), get_nexus_parameters(.y))
) %>%
  bind_rows()

dt_real

write_csv(dt_real, here("output/results/bounds_real_tb.csv"))

