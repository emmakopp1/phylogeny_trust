# ==============================================================================
# SCRIPT: visualization.R
# DESCRIPTION: Visualize the average reconstruction depth per semantic meaning
# ==============================================================================

library(here)
library(tidyverse)
library(patchwork)

# load files -------------------------------------------------------------------
summary_data_st <- read_csv(here("output/results/ancestral_reconstruction_summary_st.csv"))
summary_data_ie <- read_csv(here("output/results/ancestral_reconstruction_summary_ie.csv"))

# sino-tibetan
plot_st <- ggplot(summary_data_st, aes(x = mean_max_depth, y = mean_sinitic, label = sens)) +
  geom_text(size = 2, 
            position = position_jitter(width = 0.1, height = 0.05, seed = 123)) +
  labs(
    title = "Sino-Tibetan",
    x = "depth",
    y = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlim(0, 10)

# indo-european
plot_ie <- ggplot(summary_data_ie, aes(x = mean_max_depth, y = mean_outgroup, label = sens)) +
  geom_text(size = 2, 
            position = position_jitter(width = 0.1, height = 0.05, seed = 123)) +
  labs(
    title = "Indo-European",
    x = "depth",
    y = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlim(0, 10)


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
