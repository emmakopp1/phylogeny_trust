library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)
library(TreeTools)
library(readr)


# Add the missing "End;" line at the end of the beast tree files ---------------
write_file("End;",
           here("data/real/st_ctmc-strict-fbd-by-sens/st_ctmc-strict-fbd-by-sens.trees"),
           append = TRUE
)


write_file("End;",
           here("data/real/bantu_ctmc-strict-bd-subsample-filtered/bantu_ctmc-strict-bd-subsample-filtered.trees"),
           append = TRUE
)


write_file("End;",
           here("data/real/bantu_ctmc-strict-bd-subsample2-filtered/bantu_ctmc-strict-bd-subsample2-filtered.trees"),
           append = TRUE
)

write_file("End;",
           here("data/real/kd_ctmc-strict-bd-ht/kd_ctmc-strict-bd-ht.trees"),
           append = TRUE
)



# Create directories
dir.create(here("output/results/bantu"))
dir.create(here("output/results/bantu_subsample"))
dir.create(here("output/results/bantu_subsample2"))
dir.create(here("output/results/ie"))
dir.create(here("output/results/st"))
dir.create(here("output/results/st_by_sens"))
dir.create(here("output/results/tea"))
dir.create(here("output/results/kd"))
dir.create(here("output/trees"))

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


# Bantu -----------------------------------------------------------------------------------------------------------

phylo_bantu <- read.nexus(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.trees"))

tree_bantu = phylo_bantu[[length(phylo_bantu)]]
write.tree(tree_bantu, here("output/results/bantu/bantu_ctmc-strict-bd_tree.nex"))

ages_bantu <- get_tip_ages(phylo_bantu)
write_csv(ages_bantu, bzfile(here("output/results/bantu/bantu_ctmc-strict-bd_ages.csv.bz")))

trace_bantu <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.log"))
write_csv(trace_bantu, here("output/results/bantu/bantu_ctmc-strict-bd_tracelog.csv"))

ntipschars_bantu <- get_nexus_parameters(here("data/real/bantu_ctmc-strict-bd/bantu.nex")) |>
  mutate(family = "Bantu")


# Bantu subsample -------------------------------------------------------------------------------------------------

phylo_bantu_subset <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample-filtered/bantu_ctmc-strict-bd-subsample-filtered.trees"))

tree_bantu_subset = phylo_bantu_subset[[length(phylo_bantu_subset)]]
write.tree(tree_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tree.nex"))

ages_bantu_subset <- get_tip_ages(phylo_bantu_subset)
write_csv(ages_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tipages.csv"))

trace_bantu_subset <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample-filtered/bantu_ctmc-strict-bd-subsample-filtered.log"))
write_csv(trace_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tracelog.csv"))

ntipschars_bantu_subset <- get_nexus_parameters(here("data/real/bantu_ctmc-strict-bd-subsample/bantusubsample-filtered.nex")) |>
  mutate(family = "Bantu_subset")


# Bantu subsample 2 -----------------------------------------------------------------------------------------------

phylo_bantu_subset2 <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample2-filtered/bantu_ctmc-strict-bd-subsample2-filtered.trees"))

tree_bantu_subset2 = phylo_bantu_subset2[[length(phylo_bantu_subset2)]]
write.tree(tree_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tree.nex"))

ages_bantu_subset2 <- get_tip_ages(phylo_bantu_subset2)
write_csv(ages_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tipages.csv"))

trace_bantu_subset2 <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample2-filtered/bantu_ctmc-strict-bd-subsample2-filtered.log"))
write_csv(trace_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tracelog.csv"))

ntipschars_bantu_subset2 <- get_nexus_parameters(here("data/real/bantu_ctmc-strict-bd-subsample2/bantusubsample2-filtered.nex")) |>
  mutate(family = "Bantu_subset2")


# Indo-European ---------------------------------------------------------------------------------------------------

phylo_ie <- read.nexus(here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.trees"))

tree_ie= phylo_ie[[length(phylo_ie)]]
write.tree(tree_ie, here("output/results/ie/iecor_ctmc-strict-M1_tree.nex"))

ages_ie <- get_tip_ages(phylo_ie)
write_csv(ages_ie, bzfile(here("output/results/ie/iecor_ctmc-strict-M1_tipages.csv.bz")))

trace_ie <- parse_beast_tracelog_file(here("data/real/iecor_ctmc-strict-M1/iecor_ctmc-strict-M1.log"))
write_csv(trace_ie, here("output/results/ie/iecor_ctmc-strict-M1_tracelog.csv"))

ntipschars_ie <- get_nexus_parameters(here("data/real/iecor_ctmc-strict-M1/iecor.nex")) |>
  mutate(family = "IE")


# Sino-Tibetan ----------------------------------------------------------------------------------------------------

phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees"))

tree_st= phylo_st[[length(phylo_st)]]
write.tree(tree_st, here("output/results/st/st_ctmc-strict-fbd_tree.nex"))

ages_st <- get_tip_ages(phylo_st)
write_csv(ages_st, bzfile(here("output/results/st/st_ctmc-strict-fbd_tipages.csv.bz")))

trace_st <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.log"))
write_csv(trace_st, here("output/results/st/st_ctmc-strict-fbd_tracelog.csv"))

ntipschars_st <- get_nexus_parameters(here("data/real/st_ctmc-strict-fbd-ht/st.nex")) |>
  mutate(family = "ST")



# Sino-Tibetan by sens -----------------------------------------------------------------------------------------

phylo_st_by_sens <- read.nexus(here("data/real/st_ctmc-strict-fbd-by-sens/st_ctmc-strict-fbd-by-sens.trees"))

tree_st_by_sens = phylo_st_by_sens[[length(phylo_st_by_sens)]]
write.tree(tree_st_by_sens, here("output/results/st_by_sens/st_ctmc-strict-fbd_by_sens_tree.nex"))

ages_st_by_sens <- get_tip_ages(phylo_st_by_sens)
write_csv(ages_st_by_sens, bzfile(here("output/results/st_by_sens/st_ctmc-strict-fbd_by_sens_tipages.csv.bz")))

trace_st_by_sens <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd-by-sens/st_ctmc-strict-fbd-by-sens.log"))
write_csv(trace_st_by_sens, here("output/results/st_by_sens/st_ctmc-strict-fbd_by_sens_tracelog.csv"))

ntipschars_st_by_sens <- get_nexus_parameters(here("data/real/st_ctmc-strict-fbd-ht/st.nex")) |>
  mutate(family = "ST_by_sens")


# Transeurasian ---------------------------------------------------------------------------------------------------

phylo_tea <- read.nexus(here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees"))

tree_tea= phylo_tea[[length(phylo_tea)]]
write.tree(tree_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tree.nex"))

ages_tea <- get_tip_ages(phylo_tea)
write_csv(ages_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tipages.csv"))

trace_tea <- parse_beast_tracelog_file(here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.log"))
write_csv(trace_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tracelog.csv"))

ntipschars_tea <- get_nexus_parameters(here("data/real/tea_ctmc-strict-fbd-constrained/tea.nex")) |>
  mutate(family = "TEA")


# Kra-Dai  -----------------------------------------------------------------------------------------

phylo_kd <- read.nexus(here("data/real/kd_ctmc-strict-bd-ht/kd_ctmc-strict-bd-ht.trees"))

tree_kd = phylo_kd[[length(phylo_kd)]]
write.tree(tree_kd, here("output/results/kd/kd_ctmc-strict-bd_tree.nex"))

ages_kd <- get_tip_ages(phylo_kd)
write_csv(ages_kd, bzfile(here("output/results/kd/kd_ctmc-strict-bd_tipages.csv.bz")))

trace_kd <- parse_beast_tracelog_file(here("data/real/kd_ctmc-strict-bd-ht/kd_ctmc-strict-bd-ht.log"))
write_csv(trace_kd, here("output/results/kd/kd_ctmc-strict-bd_tracelog.csv"))

ntipschars_kd <- get_nexus_parameters(here("data/real/kd_ctmc-strict-bd-ht/kd.nex")) |>
  mutate(family = "KD")

# Number of tips and characters -----------------------------------------------------------------------------------

ntipschars <- bind_rows(ntipschars_bantu, ntipschars_bantu_subset, ntipschars_bantu_subset2, ntipschars_ie, ntipschars_st, ntipschars_st_by_sens, ntipschars_tea, ntipschars_kd) |>
  relocate(family, 1)
write_csv(ntipschars, here("output/results/ntipschars.csv"))
