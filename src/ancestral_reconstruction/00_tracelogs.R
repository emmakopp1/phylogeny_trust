# ------------------------------------------------------------------------------
# Script Name: 00_tracelogs.R
# Description: This script processes BEAST output files (trees and trace logs) for 
#              ancestral reconstruction analyses of two language families. It performs:
#                - Parsing of BEAST trace log files for parameter estimates
#                - Calculation of summary statistics with burnin removal
#                - Processing of frequency parameters and mutation rates
# ------------------------------------------------------------------------------

library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)
library(TreeTools)
library(readr)


# Get the number of taxa and traits from a nexus file
get_nexus_parameters <- function(file) {
  if (str_detect(file, "tea")) {
    phydt <- ReadAsPhyDat(file)
  } else {
    phydt <- read.nexus.data(file) |>
      PhyDat()
  }
  tibble(N = length(attributes(phydt)$names), k = length(attributes(phydt)$index))
}


# Indo-European ---------------------------------------------------------------------------------------------------

phylo_ie <- read.nexus(here("data/real/iecor_ctmc-strict-M1/IECoR_M1_CTMC_Gamma_1_Rate_For_All_Mgs_combined.trees"))
tracelog_ie <- parse_beast_tracelog_file(here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.log"))
ntipschars_ie <- get_nexus_parameters(here("data/real/iecor_ctmc-strict-M1/iecor.nex")) |>
  mutate(family = "IE")

# Sino-Tibetan ----------------------------------------------------------------------------------------------------

phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees"))
tracelog_st <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.log"))
ntipschars_st <- get_nexus_parameters(here("data/real/st_ctmc-strict-fbd-uni/st.nex")) |>
  mutate(family = "ST")

# Number of tips and characters ------------------------------------------------

ntipschars <- bind_rows(ntipschars_ie, ntipschars_st) |>
  relocate(family, 1)
write_csv(ntipschars, here("output/results/ntipschars.csv"))

# Trace logs  -------------------------------------------------------------------------

burnin <- .1

tracelog_ie <- tracelog_ie |>
  mutate(family = "IE") |> 
  rename()
tracelog_st <- tracelog_st |>
  mutate(family = "ST")


tracelog_summary <- list(tracelog_ie, tracelog_st) |>
  purrr::map(~ .x |>
               add_tally(name = "n_trees")|>
               filter(Sample > ceiling(max(Sample) * burnin)) |>
               select(family, n_trees, starts_with("freqParameter"), clockRate.c.clock , TreeHeight.t.tree, starts_with("mutationRate")) |>
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

