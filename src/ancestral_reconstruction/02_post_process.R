library(ape)
library(here)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)

# sino-tibetan -----------------------------------------------------------------
# Load dataframe
path_st <- here("output/results/ancestral_reconstruction_st.csv")
data_st <- read.csv2(path_st, header = TRUE, sep = ',')


# Load and clean linguistic data 
Y_chinese <- read.nexus.data(here("data/real/st_ctmc-strict-fbd-uni/st.nex"))
Y_chinese <- lapply(Y_chinese, function(col) replace(col, col == "?", NA))
Y_chinese <- lapply(Y_chinese, as.numeric)
Y_chinese <- as.data.frame(Y_chinese)


# Chinese languages
chinese <- c(
  "trait", 
  "SiniticBeijing",
  "SiniticChaozhou",
  "SiniticGuangzhou",
  "SiniticJieyang",
  "SiniticLonggang",
  "SiniticOldChinese",
  "SiniticXingning")

# Select Chinese languages
Y_chinese <- Y_chinese |>
  mutate(trait = row_number()) |> 
  relocate(trait, .before = everything()) |> 
  select(all_of(chinese))

# For each trait check if there is a Chinese language presence
data_st_main <- data_st |>
  left_join(Y_chinese, by = "trait") |> 
  mutate(any_sinitic = if_any(starts_with("Sinitic"), ~ .x == 1)) |> 
  relocate(any_sinitic, .after = trait)

# For each tree and sens keep the maximum depth reconstruction (value)
max_depth_data_st <- data_st_main |> 
  mutate(value = as.numeric(value)) |> 
  group_by(sens, tree) |>
  slice_max(value, n = 1, with_ties = FALSE) |>
  ungroup()

# For each sens average the maximum depth reconstruction value and the 
# presence of Chinese
summary_data_st <- max_depth_data_st |> 
  group_by(sens) |> 
  summarize(mean_max_depth = mean(value, na.rm = TRUE), 
            mean_sinitic = mean(any_sinitic, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(
    sens = sens |>
      str_replace_all("_", " ") |>
      str_remove_all("\\bthe\\b") |>
      str_remove_all("\\bto\\b") |>
      str_trim() |>
      str_squish()
  )
  
summary_data_st$mean_sinitic[is.nan(summary_data_st$mean_sinitic)] <- 0

write_csv(summary_data_st, here("output/results/ancestral_reconstruction_summary_st.csv"))

# indo-european ---------------------------------------------------------------
# Load dataframe
path_ie <- here("output/results/ancestral_reconstruction_ie.csv")
data_ie <- read.csv2(path_ie, header = TRUE, sep = ',')

# Load and clean linguistic data 
Y_ie <- read.nexus.data(here("data/real/iecor_ctmc-strict-M1/iecor.nex"))
Y_ie <- lapply(Y_ie, function(col) replace(col, col == "?", NA))
Y_ie <- lapply(Y_ie, as.numeric)
Y_ie <- as.data.frame(Y_ie)


# tocharian and anatolian languages
tocharian_anatolian <- c(
  "trait",
  "TocharianA",
  "TocharianB",
  "Hittite",
  "Luvian"
  )

# Select Chinese languages
Y_tocharian_anatolian <- Y_ie |>
  mutate(trait = row_number()) |> 
  relocate(trait, .before = everything()) |> 
  select(all_of(tocharian_anatolian))

# For each trait check if there is a Chinese language presence
data_ie_main <- data_ie |>
  left_join(Y_tocharian_anatolian, by = "trait") |> 
  mutate(any_outgroup = rowSums(across(tocharian_anatolian, ~ .x == 1), na.rm = TRUE) > 0) |> 
  relocate(any_outgroup, .after = trait)

# For each tree and sens keep the maximum depth reconstruction (value)
max_depth_data_ie <- data_ie_main |> 
  mutate(value = as.numeric(value)) |> 
  group_by(sens, tree) |>
  slice_max(value, n = 1, with_ties = FALSE) |>
  ungroup()

# For each sens average the maximum depth reconstruction value and the 
# presence of Chinese
summary_data_ie <- max_depth_data_ie |> 
  group_by(sens) |> 
  summarize(mean_max_depth = mean(value, na.rm = TRUE), 
            mean_outgroup = mean(any_outgroup, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(
    sens = sens |>
      str_replace_all("_", " ") |>
      str_remove_all("\\bthe\\b") |>
      str_remove_all("\\bto\\b") |>
      str_trim() |>
      str_squish()
  )

write_csv(summary_data_ie, here("output/results/ancestral_reconstruction_summary_ie.csv"))
