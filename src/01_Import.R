library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)
library(TreeTools)
library(readr)


# Create directories

dir.create(here("output/results/ie"))
dir.create(here("output/results/st"))
dir.create(here("output/trees")) # ?

# Get the ages for all tips of each tree in a multiPhylo object
get_tip_ages <- function(phylo) {
  ntips <- Ntip(phylo[[1]])
  map_df(1:(length(phylo)), function(i) {
    ages <- node.depth.edgelength(phylo[[i]])[1:ntips]
    depths <- round(max(ages) - ages,2)
    tibble(tree = i, tip = phylo[[i]]$tip.label, age = ages, depth = depths)
  })
}

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

tree_ie= phylo_ie[[length(phylo_ie)]]
write.tree(tree_ie, here("output/results/ie/iecor_ctmc-strict-M1_tree.nex"))

ages_ie <- get_tip_ages(phylo_ie)
write_csv(ages_ie, bzfile(here("output/results/ie/iecor_ctmc-strict-M1_tipages.csv.bz")))

trace_ie <- parse_beast_tracelog_file(here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.log"))
write_csv(trace_ie, here("output/results/ie/iecor_ctmc-strict-M1_tracelog.csv"))

ntipschars_ie <- get_nexus_parameters(here("data/real/iecor_ctmc-strict-M1/iecor.nex")) |>
  mutate(family = "IE")


# Sino-Tibetan ----------------------------------------------------------------------------------------------------

phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees"))

tree_st= phylo_st[[length(phylo_st)]]
write.tree(tree_st, here("output/results/st/st_ctmc-strict-fbd_tree.nex"))

ages_st <- get_tip_ages(phylo_st)
write_csv(ages_st, bzfile(here("output/results/st/st_ctmc-strict-fbd_tipages.csv.bz")))

trace_st <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.log"))
write_csv(trace_st, here("output/results/st/st_ctmc-strict-fbd_tracelog.csv"))

ntipschars_st <- get_nexus_parameters(here("data/real/st_ctmc-strict-fbd-uni/st.nex")) |>
  mutate(family = "ST")


# Number of tips and characters -----------------------------------------------------------------------------------

ntipschars <- bind_rows(ntipschars_ie, ntipschars_st) |>
  relocate(family, 1)
write_csv(ntipschars, here("output/results/ntipschars.csv"))
