library(ape)
library(here)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(adephylo)
library(stringr)

# sino-tibetan -----------------------------------------------------------------
# Load dataframe
path_st <- here("output/results/ancestral_reconstruction_st_final.csv")
data_st <- read.csv(path_st, header = TRUE, sep = ',', row.names = NULL)

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
# Number of trait per meaning
trait_per_meaning_ie <- read.csv(here("data/real/meanings_sets_ie.csv")) |>
  mutate(n_traits = end - start + 1) |>
  select(meaning, n_traits) |>
  rename(sens = meaning)

# Load dataframe
path_ie <- here("output/results/ancestral_reconstruction_ie_final.csv")
data_ie <- read.csv2(path_ie, header = TRUE, sep = ',') |> 
  distinct()

# proto meaning/trait value
roots_ie = read.csv(here("data/real/iecor_ctmc-strict-M1/iecor_roots.csv")) |>
  rename(cognate_id = id)

# Nexus files
nexus_text_ie <- readLines(here("data/real/iecor_ctmc-strict-M1/iecor.nex"))
# Extraire les charstatelabels
start <- which(str_detect(nexus_text_ie, "charstatelabels"))
end <- which(str_detect(nexus_text_ie, "^\\s*;"))
end <- end[end > start][1]
labels_lines <- nexus_text_ie[(start+1):(end-1)]

# Parser chaque ligne
trait_map <- labels_lines |>
  str_trim() |>
  str_remove(",$") |>
  str_match("^(\\d+)\\s+(\\S+)$") |>
  as.data.frame() |>
  setNames(c("full", "trait_num", "label")) |>
  filter(!is.na(trait_num)) |>
  mutate(
    trait_num = as.integer(trait_num),
    type = case_when(
      str_detect(label, "_group$") ~ "group",
      str_detect(label, "_cognate_") ~ "cognate",
      TRUE ~ "other"
    ),
    sens = str_extract(label, "^[^_]+"),
    cognate_id = if_else(type == "cognate",
                         as.integer(str_extract(label, "\\d+$")),
                         NA_integer_)
  )

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
  mutate(any_tocharian = if_any(c("TocharianA", "TocharianB", "Hittite", "Luvian"), ~ .x == 1)) |> 
  mutate(any_tocharian = as.integer(any_tocharian)) |> 
  relocate(any_tocharian, .after = trait)

# For each tree and sens keep the maximum depth reconstruction (value)
max_depth_data_ie <- data_ie_main |> 
  mutate(value = as.numeric(value)) |> 
  group_by(sens, tree) |>
  slice_max(value, n = 1, with_ties = FALSE) |>
  ungroup()

# For each sens average the maximum depth reconstruction value and the 
# presence of Chinese
summary_data_ie <- max_depth_data_ie |> 
  left_join(
    trait_map |> 
      select(trait_num, cognate_id) |> 
      rename(trait = trait_num),
    by = "trait"
  ) |> 
  group_by(sens) |> 
  summarize(
    mean_max_depth = mean(value, na.rm = TRUE), 
    mean_tocharian = mean(any_tocharian, na.rm = TRUE),
    cognate_id = cognate_id[which.max(value)],
    .groups = "drop"
  )|> 
  mutate(
    sens = sens |>
      str_replace_all("_", " ") |>
      str_remove_all("\\bthe\\b") |>
      str_remove_all("\\bto\\b") |>
      str_trim() |>
      str_squish()
  ) |> 
  rename(mean_outgroup = mean_tocharian) |> 
  left_join(roots_ie,by="cognate_id") |>
  select(-root_language) 

summary_data_ie$mean_outgroup[is.nan(summary_data_ie$mean_outgroup)] <- 0

write_csv(summary_data_ie, here("output/results/ancestral_reconstruction_summary_ie.csv"))
 
### compute for each trait : la   ---------------------------------------------------

# Calcul de l'âge moyen de la racine sur la posterieur
length_phylo <- 200
phylo_ie <- read.nexus(here("data/real/iecor_ctmc-strict-M1/IECoR_M1_CTMC_Gamma_1_Rate_For_All_Mgs_combined.trees"))
M_ie <- length(phylo_ie)
phylo_ie <- phylo_ie[seq(0.8 * M_ie, M_ie, length = length_phylo)]
M_ie <- length(phylo_ie)  # update M to number of retained trees

mean_root_age <- mean(sapply(1:M_ie, function(i) max(distRoot(phylo_ie[[i]]))))


# Final dataframe 
data_ie_by_trait_and_sens <- data_ie |>
  select(-node) |>
  mutate(value = as.numeric(value)) |>
  group_by(trait, sens) |>
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(trait_num = as.integer(trait)) |>
  left_join(
    trait_map |> select(trait_num, label, type, cognate_id),
    by = "trait_num"
  ) |>
  select(-trait_num, -type, -label)  |>
  group_by(sens) |>
  mutate(max_value = max(mean_value, na.rm = TRUE)) |>
  ungroup() |>
  mutate(ratio_value_max = mean_value / max_value) |>
  left_join(trait_per_meaning_ie, by = "sens") |>
  mutate(mean_value = round(mean_value, 3)) |>
  mutate(max_value = round(max_value, 3)) |>
  mutate(ratio_value_max = round(ratio_value_max, 3)) |>
  right_join(Y_tocharian_anatolian, by = "trait") |>
  mutate(any_outgroup = rowSums(across(all_of(tocharian_anatolian), ~ .x == 1), na.rm = TRUE) > 0) |>
  group_by(sens) |>
  mutate(mean_outgroup = mean(any_outgroup, na.rm = TRUE))|>
  ungroup() |>
  select(-TocharianA, -TocharianB, -Hittite, -Luvian, -any_outgroup) |>
  mutate(mean_outgroup = round(mean_outgroup,3)) |>
  left_join(roots_ie,by="cognate_id") |>
  select(-root_language) |> 
  mutate(mean_root_age = mean_root_age)

head(data_ie_by_trait_and_sens)

#data_ie_water = data_ie_by_trait_and_sens |> filter(sens == "water")

write_csv(data_ie_by_trait_and_sens, here("output/results/ancestral_reconstruction_summary_ie_by_sens.csv"))

head(data_ie_by_trait_and_sens)



