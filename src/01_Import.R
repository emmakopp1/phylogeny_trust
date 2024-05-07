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
  })
}

phylo_bantu <- read.nexus(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.trees"))
ages_bantu <- get_tip_ages(phylo_bantu)
# write_csv(ages_bantu, here("output/results/bantu_ctmc-strict-bd_tipages.csv"))
write_csv(ages_bantu, xzfile(here("output/results/bantu_ctmc-strict-bd_ages.csv.xz")))
write_rds(ages_bantu, here("output/results/bantu_ctmc-strict-bd_tipages.rds"), compress = "xz")
trace_bantu <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd/bantu_ctmc-strict-bd.log"))
write_csv(trace_bantu, here("output/results/bantu_ctmc-strict-bd_tracelog.csv"))
write_csv(ages_bantu, bzfile(here("output/results/bantu_ctmc-strict-bd_ages.csv.bz")))


phylo_bantu_subset <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.trees"))
ages_bantu_subset <- get_tip_ages(phylo_bantu_subset)
write_csv(ages_bantu_subset, here("output/results/bantu_ctmc-strict-bd-subsample_tipages.csv"))
trace_bantu_subset <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample/bantu_ctmc-strict-bd-subsample.log"))
write_csv(trace_bantu_subset, here("output/results/bantu_ctmc-strict-bd-subsample_tracelog.csv"))

phylo_bantu_subset2 <- read.nexus(here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.trees"))
ages_bantu_subset2 <- get_tip_ages(phylo_bantu_subset2)
write_csv(ages_bantu_subset2, here("output/results/bantu_ctmc-strict-bd-subsample2_tipages.csv"))
trace_bantu_subset2 <- parse_beast_tracelog_file(here("data/real/bantu_ctmc-strict-bd-subsample2/bantu_ctmc-strict-bd-subsample2.log"))
write_csv(trace_bantu_subset2, here("output/results/bantu_ctmc-strict-bd-subsample2_tracelog.csv"))
