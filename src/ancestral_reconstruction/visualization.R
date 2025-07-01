# ==============================================================================
# SCRIPT: visualization.R
# DESCRIPTION: Visualize the average reconstruction depth per semantic meaning
# ==============================================================================

library(here)
library(tidyverse)
library(patchwork)

# Load files -------------------------------------------------------------------
path_ie <- here("src/ancestral_reconstruction_ie/results.txt")
data_ie <- read.csv2(path_ie, header = TRUE, sep = '\t')

path_st <- here("src/ancestral_reconstruction_st/results.txt")
data_st <- read.csv2(path_st, header = TRUE, sep = '\t')

# Function to process reconstruction data --------------------------------------
process_reconstruction_data <- function(data) {
  data |>
    group_by(sens, tree) |>
    summarize(max_depth = max(value), .groups = "drop") |> 
    ungroup() |> 
    mutate(max_depth = as.numeric(max_depth)) |> 
    group_by(sens) |> 
    summarize(mean_max_depth = mean(max_depth, na.rm = TRUE), .groups = "drop") |> 
    ungroup() |>
    # Add random y-coordinate for text positioning (to avoid overlap)
    mutate(y_position = runif(n(), min = 0, max = 1)) |> 
    mutate(
      sens = sens |>
        str_replace_all("_", " ") |>
        str_remove_all("\\bthe\\b") |>
        str_remove_all("\\bto\\b") |>
        str_trim() |>
        str_squish()
    )
}

# Process data for both language families -------------------------------------
result_ie <- process_reconstruction_data(data_ie)
result_st <- process_reconstruction_data(data_st)

# Create plotting function ----------------------------------------------------
create_depth_plot <- function(data, title) {
  ggplot(data, aes(x = mean_max_depth, y = y_position, label = sens)) +
    geom_text(nudge_y = 0.03, size = 2) +
    labs(
      title = title,
      x = "Average maximum reconstruction depth",
      y = ""
    ) +
    theme_minimal() + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    xlim(0, 10)
}

# Create plots ----------------------------------------------------------------
plot_ie <- create_depth_plot(result_ie, "Indo-European")
plot_st <- create_depth_plot(result_st, "Sino-Tibetan")

# Combine and save plots ------------------------------------------------------
combined_plot <- plot_ie + plot_st

# Save the combined plot
ggsave(
  filename = here("output/figs/ancestral_reconstruction_by_semantic_meaning.pdf"),
  plot = combined_plot,
  width = 20,
  height = 10,
  units = "in"
)
