rm(list=ls())
library(ape)
library(here)
library(phangorn)
library(phytools)
library(TreeTools)
library(tracerer)
library(tidyverse)

write_binary_nexus <- function(x, file) {
  write.phyDat(x, file, format = "nexus")
  read_lines(file) |>
    str_replace('symbols="0123456789"', 'symbols="01"') |>
    str_subset("^(?!\\[Data written by)") |>
    write_lines(file)
}


df <- read.nexus.data(here("data/real/st_ctmc-strict-fbd-by-sens/st.nex")) |>
  purrr::map(unlist) |>
  purrr::map_df(~ as.data.frame(t(.x)), .id = "Taxon") |>
  column_to_rownames("Taxon") |>
  dplyr::select(-c(233, 234, 2023, 2024)) |>  
  as.matrix() |>
  MatrixToPhyDat()


write_binary_nexus(
  df,
  here("data/real/st_ctmc-strict-fbd-by-sens-no-outliers/st-no-outlier.nex")
)








