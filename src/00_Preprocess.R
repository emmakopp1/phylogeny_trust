library(ape)
library(tidyverse)
library(here)

# Registration process : delete nulls columns --------------------------------------------------------------------
# Bantu subsample
read.nexus.data(here("data/real/bantu_ctmc-strict-bd-subsample-filtered/bantusubsample.nex"))|>
  as_tibble() |>
  mutate(across(everything(), ~ na_if(.x, "?"))) |>
  mutate(across(everything(), as.numeric))  |>
  mutate(total = rowSums(across(everything()), na.rm=T)) |>
  filter(total != 0) |> 
  select(-total) |>
  mutate(across(everything(), as.character))  |>
  mutate(across(everything(), ~ replace_na(.x, "?"))) |>
  write.nexus.data(
    here("data/real/bantu_ctmc-strict-bd-subsample-filtered/bantusubsample-filtered.nex"),
    format = "STANDARD")

  

# Bantu subsample 2 
read.nexus.data(here("data/real/bantu_ctmc-strict-bd-subsample2-filtered/bantusubsample2.nex")) |>
  as.tibble() |>
  mutate(across(everything(), ~ na_if(.x, "?"))) |>
  mutate(across(everything(), as.numeric))  |>
  mutate(total = rowSums(across(everything()), na.rm=T)) |>
  filter(total != 0) |> 
  select(-total) |>
  mutate(across(everything(), as.character))  |>
  mutate(across(everything(), ~ replace_na(.x, "?"))) |>
  write.nexus.data(
    here("data/real/bantu_ctmc-strict-bd-subsample2-filtered/bantusubsample2-filtered.nex"),
    format = "STANDARD")






