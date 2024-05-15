library(here)
library(tidyverse)
library(ape)
library(treeio)
library(tracerer)

get_tip_ages <- function(phylo) {
  ntips <- Ntip(phylo[[1]])
  map_df(1:(length(phylo)), function(i) {
    ages <- node.depth.edgelength(phylo[[i]])[1:ntips]
    tibble(tree = i, tip = phylo[[i]]$tip.label, age = ages)
  })%>%
    group_by(tip) %>%
    summarise(age = mean(age, na.rm = TRUE))
}

# Bantu
phylo_bantu <- read.nexus(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.trees"))
ages_bantu <- get_tip_ages(phylo_bantu)
# write_csv(ages_bantu, here("output/results/bantu_ctmc-strict-bd_tipages.csv"))
#write_csv(ages_bantu, xzfile(here("output/results/bantu_ctmc-strict-bd_ages.csv.xz"))) #can't open after
write_rds(ages_bantu, here("output/results/bantu/bantu_ctmc-strict-bd_tipages.rds"), compress = "xz")
trace_bantu <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.log"))
write_csv(trace_bantu, here("output/results/bantu/bantu_ctmc-strict-bd_tracelog.csv"))
write_csv(ages_bantu, bzfile(here("output/results/bantu/bantu_ctmc-strict-bd_ages.csv.bz")))

# Bantu subsample
phylo_bantu_subset <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.trees"))
ages_bantu_subset <- get_tip_ages(phylo_bantu_subset)
write_csv(ages_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tipages.csv"))
trace_bantu_subset <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.log"))
write_csv(trace_bantu_subset, here("output/results/bantu_subsample/bantu_ctmc-strict-bd-subsample_tracelog.csv"))

# Bantu subsample 2
phylo_bantu_subset2 <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.trees"))
ages_bantu_subset2 <- get_tip_ages(phylo_bantu_subset2)
write_csv(ages_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tipages.csv"))
trace_bantu_subset2 <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.log"))
write_csv(trace_bantu_subset2, here("output/results/bantu_subsample2/bantu_ctmc-strict-bd-subsample2_tracelog.csv"))

# Sino-tibetain 
phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees"))
ages_st <- get_tip_ages(phylo_st)
write_csv(ages_st, here("output/results/st/st_ctmc-strict-fbd_tipages.csv"))
trace_st <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.log"))
write_csv(trace_st, here("output/results/st/st_ctmc-strict-fbd_tracelog.csv"))

# Transeurasien
phylo_tea <- read.nexus(here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees"))
ages_tea <- get_tip_ages(phylo_tea)
write_csv(ages_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tipages.csv"))
trace_tea <- parse_beast_tracelog_file(here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.log"))
write_csv(trace_tea, here("output/results/tea/tea_ctmc-strict-fbd-constrained_tracelog.csv"))



