# ------------------------------------------------------------------------------
# Script Name: visualization.R
# Run : local
# Description: Creates heatmap visualizations for shared cognate analysis results
# -----------------------------------------------------------------------------------------

library(here)
library(pheatmap)
library(tidyverse)

# load data --------------------------------------------------------------------
results <- readRDS(here("output/results/shared_cognates.rds"))

# dataset for age 1
pdf(here("output/figs/shared_cognate_heatmap_1.pdf"), width = 10, height = 10)
pheatmap(as.matrix(results$analysis[[1]]$similarity_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = str_glue(""))
dev.off()


# dataset for age 9
pdf(here("output/figs/shared_cognate_heatmap_9.pdf"), width = 10, height = 10)
pheatmap(as.matrix(results$analysis[[9]]$similarity_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = str_glue(""))
dev.off()

# dataset for age 17
pdf(here("output/figs/shared_cognate_heatmap_17.pdf"), width = 10, height = 10)
pheatmap(as.matrix(results$analysis[[17]]$similarity_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = str_glue(""))
dev.off()

# summary table 
summary_table <- read_csv(here("output/results/shared_cognate_summary_table.csv"))

# number of shared cognates between the outgroup and the ingroup
plot_summary_table <- ggplot(summary_table, aes(x = age, y = n_shared_mean)) +
  geom_line(color = "darkblue", size = 1) +
  geom_point(color = "darkblue", size = 2) +
  labs(
    x = "age (millenia)",
    y = "",
    title = ""
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = 1:17)

# save results
ggsave(here("output/figs/shared_cognate_outgroup.pdf"), 
       plot = plot_summary_table, 
       width = 8, 
       height = 6, 
       units = "in")



