library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)
library(TreeTools)

dir.create(here("output/results/bantu"))
dir.create(here("output/results/bantu_subsample"))
dir.create(here("output/results/bantu_subsample2"))
dir.create(here("output/results/ie"))
dir.create(here("output/results/st"))
dir.create(here("output/results/tea"))

# Get the ages for all tips of each tree in a multiPhylo object
get_tip_ages <- function(phylo) {
  ntips <- Ntip(phylo[[1]])
  map_df(1:(length(phylo)), function(i) {
    ages <- node.depth.edgelength(phylo[[i]])[1:ntips]
    tibble(tree = i, tip = phylo[[i]]$tip.label, age = ages)
  })
  # %>%
  #   group_by(tip) %>%
  #   summarise(age = mean(age, na.rm = TRUE))
}

# Get the number of taxa and traits from a nexus file
get_nexus_parameters <- function(file) {
  phydt <- read.nexus.data(file) |> 
    PhyDat()
  tibble(N = length(attributes(phydt)$names), k = length(attributes(phydt)$index))
}

# Bantu -----------------------------------------------------------------------------------------------------------

phylo_bantu <- read.nexus(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.trees"))

ages_bantu <- get_tip_ages(phylo_bantu)
write_csv(ages_bantu, bzfile(here("output/results/bantu/bantu_ctmc-strict-bd_ages.csv.bz")))

trace_bantu <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.log"))
write_csv(trace_bantu, here("output/results/bantu/bantu_ctmc-strict-bd_tracelog.csv"))

ntipschars_bantu <- get_nexus_parameters(here("data/real/bantu_ctmc-strict-bd/bantu.nex")) |>
  mutate(family = "Bantu")


# Bantu subsample -------------------------------------------------------------------------------------------------

phylo_bantu_subset <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.trees"))

ages_bantu_subset <- get_tip_ages(phylo_bantu_subset)
write_csv(ages_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tipages.csv"))

trace_bantu_subset <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.log"))
write_csv(trace_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tracelog.csv"))

ntipschars_bantu_subset <- get_nexus_parameters(here("data/real/bantu_ctmc-strict-bd-subsample/bantusubsample.nex")) |>
  mutate(family = "Bantu_subset")


# Bantu subsample 2 -----------------------------------------------------------------------------------------------

phylo_bantu_subset2 <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.trees"))

ages_bantu_subset2 <- get_tip_ages(phylo_bantu_subset2)
write_csv(ages_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tipages.csv"))

trace_bantu_subset2 <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.log"))
write_csv(trace_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tracelog.csv"))

ntipschars_bantu_subset2 <- get_nexus_parameters(here("data/real/bantu_ctmc-strict-bd-subsample2/bantusubsample2.nex")) |>
  mutate(family = "Bantu_subset2")


# Indo-European ---------------------------------------------------------------------------------------------------

phylo_ie <- read.nexus(here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.trees"))

ages_ie <- get_tip_ages(phylo_ie)
write_csv(ages_ie, bzfile(here("output/results/ie/iecor_ctmc-strict-M1_tipages.csv.bz")))

trace_ie <- parse_beast_tracelog_file(here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.log"))
write_csv(trace_ie, here("output/results/ie/iecor_ctmc-strict-M1.csv"))

ntipschars_ie <- get_nexus_parameters(here("data/real/iecor_ctmc-strict-M1/iecor.nex")) |>
  mutate(family = "IE")


# Sino-Tibetan ----------------------------------------------------------------------------------------------------

phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees"))

ages_st <- get_tip_ages(phylo_st)
write_csv(ages_st, bzfile(here("output/results/st/st_ctmc-strict-fbd_tipages.csv.bz")))

trace_st <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.log"))
write_csv(trace_st, here("output/results/st/st_ctmc-strict-fbd_tracelog.csv"))

ntipschars_st <- get_nexus_parameters(here("data/real/st_ctmc-strict-fbd/st.nex")) |>
  mutate(family = "ST")


# Transeurasian ---------------------------------------------------------------------------------------------------

phylo_tea <- read.nexus(here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees"))

ages_tea <- get_tip_ages(phylo_tea)
write_csv(ages_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tipages.csv"))

trace_tea <- parse_beast_tracelog_file(here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.log"))
write_csv(trace_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tracelog.csv"))

ntipschars_tea <- get_nexus_parameters(here("data/real/tea_ctmc-strict-fbd-constrained/tea.nex")) |>
  mutate(family = "TEA")


# Number of tips and characters -----------------------------------------------------------------------------------

ntipschars <- bind_rows(ntipschars_bantu, ntipschars_bantu_subset, ntipschars_bantu_subset2, ntipschars_ie, ntipschars_st, ntipschars_tea) |>
  relocate(family, 1)
write_csv(ntipschars, here("output/results/ntipschars.csv"))
