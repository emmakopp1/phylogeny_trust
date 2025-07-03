# ------------------------------------------------------------------------------
# Script Name: 04_visualization.R
# Description: This script generates all the visualizations for the simulation 
#              study. It uses pre-processed results to create publication-ready 
#              plots comparing phylogenetic reconstruction performance across:
#                - Tree summarization methods (MCC vs. consensus)
#                - Reconstruction accuracy of the first split (outgroup)
#                - Total number of internal nodes per method and age
#                - Proportion of true, false, and uncertain clades
#              The figures are exported as PDF files for inclusion in reports.
# -----------------------------------------------------------------------------------------
library(here)
library(patchwork)
library(tidyverse)
library(ape)

# load data --------------------------------------------------------------------
df_number_of_nodes_avg =  read_csv(here("output/results/number_of_nodes_summary.csv"))
prob_first_split_mcc = read_csv(here("output/results/prob_first_split_mcc.csv"))
marginal_probability_first_split_ic = read_csv(here("output/results/marginal_probability_first_split_ic.csv"))

# count the number of true, false and uncertain node from the true to consensus tree
count_true_to_cs = readRDS(here("output/results/count_true_to_cs.rds")) |>
  group_by(age) |>
  mutate(
    mean_present = mean_n[value == "1"],
    mean_uncertain = mean_n[value == "2"],
    mean_absent = mean_n[value == "0"],
    y_label = case_when(
      value == "1" ~ 0,
      value == "2" ~ mean_present + mean_uncertain / 2,
      value == "0" ~ mean_present + mean_uncertain + mean_absent # total = haut
    )
  ) |>
  ungroup() 

# number of well reconstructed node in the mcc tree
count_true_to_mcc = readRDS(here("output/results/count_true_to_mcc.csv"))

# count the numer of true node from de consensus tree to the true tree
count_cs_to_true = read_csv(here("output/results/count_cs_to_true.csv"))

# proportion of true,false and uncertain nodes in the consensus tree (true -> consensus)
prop_true_to_cs = read_csv(here("output/results/prop_true_to_cs.csv"))

# pour le plot 4.bis mcc
prop_mcc_to_true = read_csv(here("output/results/prop_mcc_to_true.csv"))

# proportion of true, false node from the consensus to the true tree
prop_cs_to_true = read_csv(here("output/results/prop_cs_to_true.csv"))  

# proportion of true, false node from the mcc to the true tree
prop_mcc_to_true = read_csv(here("output/results/prop_mcc_to_true.csv"))

# plots ------------------------------------------------------------------------
# 1. Marginal probability of the first split in the posterior ------------------
plot_marginal_probability_first_split_ic <- ggplot(
  marginal_probability_first_split_ic,
  aes(x = age, y = prob)
) +
  geom_ribbon(aes(ymin = inf, ymax = sup), fill = "skyblue", alpha = 0.5) +
  geom_line(color = "darkblue", size = 1) +
  geom_point(color = "darkblue", size = 2) +
  labs(
    title = "Marginal probability of the first split",
    x = "age",
    y = "probability"
  ) +
  theme_minimal(base_size = 12)

plot_marginal_probability_first_split_ic
ggsave(here("output/figs/marginal_probability_first_split_ic.pdf"))

# 2. plot of the marginal probability of first_split in the mcc & consensus --------
plot_mcc_posterior_prob <- ggplot(prob_first_split_mcc, aes(x = age, y = mean_mcc_prob, color = "MCC")) +
  geom_line(size = 1, show.legend = FALSE) +
  geom_point(size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("MCC" = "darkblue")) +
  labs(
    title = "MCC",
    x = "age",
    y = "proportion of good outgroup",
    color = "" # Pour ne pas afficher "color" dans la légende
  ) +
  theme_minimal(base_size = 10)

plot_mcc_posterior_prob
ggsave(here("output/figs/marginal_probability_first_split_mcc.pdf"))



# 3. number of nodes in the summary tree ---------------------------------------
my_colors <- c("n_mcc" = "darkred", "n_consensus" = "darkblue")

plot_number_of_nodes = ggplot(df_number_of_nodes_avg, aes(x = age)) +
  geom_line(aes(y = n_mcc, color = "n_mcc"), size = 1) +
  geom_point(aes(y = n_mcc, color = "n_mcc"), size = 2) +
  geom_line(aes(y = n_consensus, color = "n_consensus"), size = 1) +
  geom_point(aes(y = n_consensus, color = "n_consensus"), size = 2) +
  scale_color_manual(values = my_colors, labels = c("MCC", "Consensus")) +
  labs(
    title = "",
    x = "age",
    y = "number of nodes",
    color = "method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10)
  )

plot_number_of_nodes
ggsave(here("output/figs/number_of_nodes_consensus_mcc.pdf"))

# 4. true, false and uncertain nodes in the summary tree (summary to true) -----
# consensus
plot_cs_incertain <- ggplot(count_true_to_cs, aes(x = factor(age), y = mean_n, fill = value)) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = y_label, label = round(mean_n, 0)),
    color = "white", size = 3,
    vjust = case_when(
      count_true_to_cs$value == "1" ~ -0.3,
      count_true_to_cs$value == "0" ~ 1.3,
      TRUE ~ 0
    )
  ) +
  labs(
    title = "Consensus",
    x = "Age",
    y = "Average number of nodes",
    fill = "Value"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue", "2" = "darkorange"),
    breaks = c("0", "2", "1"),
    labels = c("Absent", "Uncertain", "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()


# mcc
plot_mcc_incertain <- ggplot(count_true_to_mcc, aes(x = factor(age), y = n_mean, fill = exist)) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = y_label, label = round(n_mean, 0)),
    vjust = ifelse(count_true_to_mcc$exist == "1", -0.5, 1.2), # vers le bas ou vers le haut
    color = "white", size = 3
  ) +
  labs(
    title = "MCC",
    x = "age",
    y = "number of nodes",
    fill = "Existence"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()

plot_cs_incertain + plot_mcc_incertain

ggsave(here("output/figs/barplot_resume_to_true.pdf"))

# 5. Plot of number of true nodes in the summary tree (true to summary) ---------
# mcc
plot_count_true_to_mcc <- ggplot(count_true_to_mcc, aes(x = factor(age), y = n_mean, fill = exist)) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = y_label, label = round(n_mean, 0)),
    vjust = ifelse(count_true_to_mcc$exist == "1", -0.5, 1.2), # vers le bas ou vers le haut
    color = "white", size = 3
  ) +
  labs(
    title = "MCC",
    x = "age",
    y = "number of nodes",
    fill = "Existence"
  ) +
  ylim(0,50) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()

# consensus
plot_count_cs_to_true <- ggplot(count_cs_to_true, aes(x = factor(age), y = n_mean, fill = factor(exist))) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = y_label, label = round(n_mean, 0)),
    vjust = ifelse(count_cs_to_true$exist == "1", -0.5, 1.2), # vers le bas ou vers le haut
    color = "white", size = 3
  ) +
  labs(
    title = "Consensus",
    x = "age",
    y = "number of nodes",
    fill = "Existence"
  ) +
  ylim(0,50) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()

plot_count_cs_to_true + plot_count_true_to_mcc

ggsave(here("output/figs/barplot_true_to_resume.pdf"))

# 4bis. true, false and uncertain nodes in the summary tree (summary to true) ----
# consensus data

plot_prop_true_to_cs <- ggplot(prop_true_to_cs, aes(x = factor(age), y = mean_n, fill = factor(value))) +
  geom_col(position = "fill") +
  labs(
    title = "Consensus",
    x = "Age",
    y = "Average number of nodes",
    fill = "Value"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred","2" = "darkorange", "1" = "darkblue"),
    breaks = c("0", "2", "1"),
    labels = c("Absent", "Uncertain", "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()


# plot
plot_prop_mcc_to_true <- ggplot(prop_mcc_to_true, aes(x = factor(age), y = n_mean, fill = factor(exist))) +
  geom_col(position = "fill") +
  labs(
    title = "MCC",
    x = "age",
    y = "number of nodes",
    fill = "Existence"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()


plot_prop_true_to_cs + plot_prop_mcc_to_true

ggsave(here("output/figs/barplot_prop_true_to_resume.pdf"))

# 5bis. Plot of number of true nodes in the summary tree (true to resume) ------
# consensus
# plots
plot_prop_mcc_to_true <- ggplot(prop_mcc_to_true, aes(x = factor(age), y = n_mean, fill = factor(exist))) +
  geom_col(position = "fill") +
  labs(
    title = "MCC",
    x = "age",
    y = "number of nodes",
    fill = "Existence"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()


plot_prop_cs_to_true <- ggplot(prop_cs_to_true, aes(x = factor(age), y = n_mean, fill = factor(exist))) +
  geom_col(position = "fill") +
  labs(
    title = "Consensus",
    x = "age",
    y = "number of nodes",
    fill = "Existence"
  ) +
  scale_fill_manual(
    values = c("0" = "darkred", "1" = "darkblue"),
    labels = c("0" = "Absent", "1" = "Present")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal()

plot_prop_cs_to_true + plot_prop_mcc_to_true

ggsave(here("output/figs/barplot_prop_resume_to_true.pdf"))

# 6. plot of one plausible rake node and one not plausible rake node in the consensus
# import one consensus and one true tree of age 10 simulation 1 
tree_cs <- read.tree(here('data/simulated-2025-05-13/beast-data-sim-1/beast-data-sim-1-10/consensus-10.tree'))
tree_true <- read.tree(here('data/simulated-2025-05-13/beast-data-sim-1/beast-data-sim-1-10/tree-sim-1-10.tree'))

# analyse a plausible rake node
node <- 62
descendant <- Descendants(tree_true, node)[[1]]
tip <- tree_true$tip.label[descendant]
mrca <- getMRCA(tree_cs, tip)

# analyse a not plausible rake node
node_67 <- 67
descendant_67 <- Descendants(tree_true, node_67)[[1]]
tip_67 <- tree_true$tip.label[descendant_67]
mrca_67 <- getMRCA(tree_cs, tip_67)

# tip colors
tip_colors_tt <- rep("black", length(tree_true$tip.label))
# node 62 in blue
tip_colors_tt[descendant] <- "blue"        
# node 67 in red
tip_colors_tt[descendant_67] <- "red"       
tip_colors_cs <- rep("black", length(tree_cs$tip.label))

# identify tip position in the consensus tree
tip_positions_cs <- match(tip, tree_cs$tip.label)
tip_colors_cs[tip_positions_cs] <- "blue"   

tip_positions_cs_67 <- match(tip_67, tree_cs$tip.label)
tip_colors_cs[tip_positions_cs_67] <- "red"

par(mfrow=c(1,2))
# true tree
plot(tree_true, tip.color = tip_colors_tt, cex = 0.7)
nodelabels(node = node, frame = 'circle', cex = 0.5)
nodelabels(node = node_67, frame = 'circle', cex = 0.5)

# consensus tree
plot(tree_cs, direction = "leftwards", tip.color = tip_colors_cs, cex=0.7)
nodelabels(node = mrca, frame = 'circle', cex = 0.5)
nodelabels(node = mrca_67, frame = 'circle', cex = 0.5)

ggsave(here("output/figs/plausible_node.pdf"))
