library(tidyverse)
library(here)
library(ggh4x)
library(ggrepel)
library(khroma)
library(knitr)
library(patchwork)
library(treeio)
library(TreeTools)
library(phangorn)
library(ggtree)

width <- 13.5
height <- 19
base_font <- "Noto Sans Condensed"
base_font2 <- "Noto Sans ExtraCondensed"
plt <- color("vibrant")(3)
plt2 <- color("highcontrast")(3)
theme_set(
  theme_minimal(base_family = base_font, base_size = 9) +
    theme(
      aspect.ratio = .618,
      strip.text = element_text(size = 9),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = .35)
    )
)

N_traits <- 3000
pi1 <- 0.94305
pi0 <- 0.05695
clock_rate <- 0.018
Q <- matrix(
  c(
    -clock_rate / (2 * pi0),
    clock_rate / (2 * pi0),
    clock_rate / (2 * pi1),
    -clock_rate / (2 * pi1)
  ),
  nrow = 2,
  byrow = TRUE
)

shared_cognates <- readRDS(here("output/results/shared_cognates.rds")) |>
  rename(age = tree_age, prop = value)
shared_cognates |>
  ggplot() +
  geom_pointpath(
    aes(x = age, y = prop),
    color = plt[2],
    linewidth = 1,
    stroke = .1
  ) +
  geom_errorbar(
    aes(x = age, ymin = inf, ymax = sup),
    color = plt[2],
    width = .2
  ) +
  ylab("Proportion of shared cognates") +
  xlab("Age of the most recent common ancestor (ka)") +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_x_continuous(breaks = seq(0, 17, 1)) +
  coord_cartesian(clip = "off")
ggsave(
  here("output/figs/shared_cognate_no_homoplasie.pdf"),
  width = width,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/shared_cognate_no_homoplasie.pdf"))

shared_cognates_thq <- read_csv(here(
  "output/figs/shared_cognate_thq_no_homoplasie.csv"
)) |>
  rename(age = tree_age, prop = S_root) 

shared_cognates_thq |>
  ggplot() +
  geom_path(
    aes(x = age, y = prop),
    color = plt[2],
    linewidth = 1,
    stroke = .1
  ) +
  ylab("Proportion of shared cognates without homoplasie") +
  xlab("Age of the most recent common ancestor (ka)") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 17, 5)) +
  geom_segment(
    aes(
      x = as.numeric(8.60856),
      xend = as.numeric(8.60856),
      y = as.numeric(0.3113052 + pi1),
      yend = -Inf
    ),
    linetype = "dashed",
    linewidth = .35,
  ) +
  geom_segment(
    aes(
      x = -Inf,
      xend = as.numeric(8.60856),
      y = as.numeric(0.3113052 ),
      yend = as.numeric(0.3113052)
    ),
    linetype = "dashed",
    linewidth = .35,
  ) +
  coord_cartesian(clip = "off") +
  scale_color_vibrant(
    name = NULL,
    guide = guide_legend(
      override.aes = list(linewidth = 1.5)
    )
  ) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 1)
  )
ggsave(
  here("output/figs/shared_cognate_thq_no_homoplasie.pdf"),
  width = width * .8,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/shared_cognate_thq_no_homoplasie.pdf"))

prop_shared_tip_pair <- readRDS(here(
  "output/results/shared_cognate_tip_pair.rds"
)) |>
  rowid_to_column("age")

prop_shared_tip_pair |>
  ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0, color = "Computed")) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 0,
    color = "Glottochronology with\nestimated rate"
  )) +
  geom_segment(
    aes(
      x = -Inf,
      xend = min(filter(prop_shared_tip_pair, prop <= .1)$age),
      y = .1,
      yend = .1
    ),
    linetype = "dashed",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = min(filter(prop_shared_tip_pair, prop <= .1)$age),
      xend = min(filter(prop_shared_tip_pair, prop <= .1)$age),
      y = .1,
      yend = -Inf
    ),
    linetype = "dashed",
    linewidth = .35
  ) +
  geom_line(
    data = tibble(age = seq(0, 17, .1)),
    aes(x = age, y = exp(-2 * Q[1, 2] * age)),
    linewidth = .75,
    col = plt[2]
  ) +
  # geom_line(data = tibble(age = seq(0, 17, .1)), aes(x = age, y = exp(-2 * .14 * age)), col = plt[1]) +
  geom_pointpath(
    aes(x = age, y = prop),
    color = plt[1],
    linewidth = 1,
    stroke = .1
  ) +
  ylab("Proportion of shared cognates") +
  xlab("Age of the most recent common ancestor (ka)") +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_x_continuous(breaks = seq(0, 17, 1)) +
  scale_color_vibrant(
    name = NULL,
    guide = guide_legend(
      override.aes = list(linewidth = 1.5)
    )
  ) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 1)
  )

ggsave(
  here("output/figs/shared_cognate_tip_pair.pdf"),
  width = width,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/shared_cognate_tip_pair.pdf"))

shared_cognates |>
  mutate(type = "any pair") |>
  bind_rows(
    prop_shared_tip_pair |> mutate(type = "sample pair")
  ) |>
  ggplot() +
  geom_segment(
    aes(
      x = -Inf,
      xend = min(filter(prop_shared_tip_pair, prop <= .1)$age),
      y = .1,
      yend = .1
    ),
    linetype = "dashed",
    color = "grey40",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = min(filter(prop_shared_tip_pair, prop <= .1)$age),
      xend = min(filter(prop_shared_tip_pair, prop <= .1)$age),
      y = .1,
      yend = -Inf
    ),
    linetype = "dashed",
    color = "grey40",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = -Inf,
      xend = min(filter(prop_shared_tip_pair, prop <= .5)$age),
      y = .5,
      yend = .5
    ),
    linetype = "dotted",
    color = "grey40",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = min(filter(prop_shared_tip_pair, prop <= .5)$age),
      xend = min(filter(prop_shared_tip_pair, prop <= .5)$age),
      y = .5,
      yend = -Inf
    ),
    linetype = "dotted",
    color = "grey40",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = -Inf,
      xend = min(filter(shared_cognates, prop <= .5)$age),
      y = .5,
      yend = .5
    ),
    linetype = "twodash",
    color = "grey40",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = min(filter(shared_cognates, prop <= .5)$age),
      xend = min(filter(shared_cognates, prop <= .5)$age),
      y = .5,
      yend = -Inf
    ),
    linetype = "twodash",
    color = "grey40",
    linewidth = .35
  ) +
  geom_pointpath(
    aes(x = age, y = prop, color = type),
    linewidth = 1,
    stroke = .1
  ) +
  geom_errorbar(
    aes(x = age, ymin = inf, ymax = sup, color = type),
    width = .2
  ) +
  scale_color_vibrant(name = NULL) +
  ylab(
    "Proportion of cognates shared\nbetween languages across the 1st split"
  ) +
  xlab("Age of the most recent common ancestor (ka BP)") +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_x_continuous(breaks = seq(0, 17, 1)) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 1)
  )
ggsave(
  here("output/figs/shared_cognates.pdf"),
  width = width,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/shared_cognates.pdf"))

prop_true_to_cs <- read_csv(here("output/results/prop_true_to_cs.csv")) |>
  mutate(type = "Consensus", value = as.character(value))
prop_cs_to_true <- read_csv(here("output/results/prop_cs_to_true.csv")) |>
  rename(value = exist, mean_n = n_mean) |>
  mutate(type = "Consensus", value = as.character(value))
prop_mcc_to_true <- read_csv(here("output/results/prop_mcc_to_true.csv")) |>
  rename(value = exist, mean_n = n_mean) |>
  mutate(type = "MCC", value = as.character(value))
prop_hipstr_to_true <- read_csv(here(
  "output/results/prop_hipstr_to_true.csv"
)) |>
  rename(value = exist, mean_n = n_mean) |>
  mutate(type = "HIPSTR", value = as.character(value))

bind_rows(prop_true_to_cs, prop_mcc_to_true, prop_hipstr_to_true) |>
  mutate(type = fct_inorder(type)) |>
  mutate(
    value = case_when(
      value == "1" ~ "Concordant",
      value == "2" ~ "Reconcilable",
      value == "0" ~ "Discordant"
    )
  ) |>
  mutate(
    value = factor(
      value,
      levels = c("Discordant", "Reconcilable", "Concordant")
    )
  ) |>
  ggplot(aes(x = (age), y = mean_n, fill = value)) +
  geom_col(position = "fill", linewidth = .15) +
  geom_hline(
    yintercept = .5,
    linetype = "dashed",
    linewidth = .5,
    color = "white"
  ) +
  labs(
    x = "Age (ka BP)",
    y = "Average proportion\nof nodes",
    fill = ""
  ) +
  scale_fill_highcontrast(reverse = TRUE) +
  scale_color_highcontrast(reverse = TRUE) +
  # scale_x_continuous(breaks = c(5,10,15)) +
  # scale_x_discrete(
  #   guide = guide_axis(check.overlap = TRUE)
  # ) +
  coord_cartesian(clip = "off", expand = FALSE) +
  facet_wrap(~type) +
  theme(
    # axis.text = element_text(family = base_font2, size = 8),
    # axis.text.y.left = element_text(margin = margin(r = 0.175, unit = "lines")),
    # axis.text.x.bottom = element_text(margin = margin(t = 0.25, unit = "lines")),
    axis.ticks = element_line(size = .25, color = "grey40"),
    axis.ticks.length = unit(0.15, "lines"),
    # axis.title = element_text(size = 9),
    # axis.title.y.left = element_text(margin = margin(b = 0.25, unit = "lines")),
    # legend.text = element_text(size = 10),
    legend.key.size = unit(.75, "line"),
    # panel.spacing.x = unit(0.15, "lines"),
    legend.position = "bottom",
    legend.margin = margin(t = -.5, r = 0, b = 0, l = 0, unit = "lines"),
  )
ggsave(
  here("output/figs/barplot_prop_true_to_resume.pdf"),
  width = width,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/barplot_prop_true_to_resume.pdf"))

bind_rows(prop_cs_to_true, prop_mcc_to_true, prop_hipstr_to_true) |>
  mutate(type = fct_inorder(type)) |>
  mutate(
    value = case_when(
      value == "1" ~ "Concordant",
      value == "2" ~ "Reconcialable",
      value == "0" ~ "Discordant"
    )
  ) |>
  mutate(
    value = factor(
      value,
      levels = c("Discordant", "Reconcilable", "Concordant")
    )
  ) |>
  ggplot(aes(x = age, y = mean_n, fill = value)) +
  geom_col(position = "fill", linewidth = .15) +
  geom_hline(
    yintercept = .5,
    linetype = "dashed",
    linewidth = .5,
    color = "white"
  ) +
  labs(
    x = "Age (ka BP)",
    y = "Average proportion\nof nodes",
    fill = ""
  ) +
  scale_fill_manual(values = rev(color("high contrast")(3)[-2])) +
  coord_cartesian(clip = "off", expand = FALSE) +
  facet_wrap(~type) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -.5, r = 0, b = 0, l = 0, unit = "lines"),
  )
ggsave(
  here("output/figs/barplot_prop_resume_to_true.pdf"),
  width = width,
  height = height * 1.25,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/barplot_prop_resume_to_true.pdf"))

summary_data_st <- read_csv(here(
  "output/results/ancestral_reconstruction_summary_st.csv"
)) |>
  mutate(family = "Sino-Tibetan") |>
  rename(mean_outgroup = mean_sinitic)
summary_data_ie <- read_csv(here(
  "output/results/ancestral_reconstruction_summary_ie.csv"
)) |>
  mutate(family = "Indo-European")

concepts <- bind_rows(summary_data_st, summary_data_ie) |>
  filter(!is.na(mean_outgroup)) |>
  mutate(x = mean_outgroup >= .5) |>
  group_by(family, x) |>
  filter(
    mean_max_depth == max(mean_max_depth, na.rm = TRUE) |
      mean_max_depth == min(mean_max_depth, na.rm = TRUE)
  ) |>
  mutate(
    max = mean_max_depth == max(mean_max_depth, na.rm = TRUE),
    min = mean_max_depth == min(mean_max_depth, na.rm = TRUE)
  ) |>
  ungroup() |>
  group_by(family, x, max, min) |>
  slice(1) |>
  ungroup() |>
  mutate(sens = str_remove_all(sens, "hide ")) |>
  mutate(sens = str_replace_all(sens, " of weight", "\n(of weight)")) |>
  mutate(sens = str_replace_all(sens, "I first person singular", "1SG"))

bind_rows(summary_data_st, summary_data_ie) |>
  ggplot() +
  geom_point(
    aes(x = mean_max_depth, y = mean_outgroup),
    color = plt[2],
    alpha = .5,
    size = 1
  ) +
  geom_point(
    data = concepts,
    aes(x = mean_max_depth, y = mean_outgroup),
    color = plt[1],
    size = 1
  ) +
  geom_text_repel(
    data = concepts,
    aes(x = mean_max_depth, y = mean_outgroup, label = sens),
    seed = 123,
    min.segment.length = 0,
    segment.size = .35,
    # point.padding = .5,
    family = base_font,
    lineheight = .8,
    color = plt[1],
    bg.color = "white",
    bg.r = 0.05,
    size = 9 / .pt
  ) +
  facet_wrap(~family, scales = "free") +
  coord_cartesian(clip = "off") +
  xlab("Mean maximum age (ka BP)") +
  ylab("Probability of presence in\nthe early-diverging lineage")
ggsave(
  here("output/figs/ancestral_reconstruction_by_semantic_meaning.pdf"),
  width = width,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/ancestral_reconstruction_by_semantic_meaning.pdf"))

# summary_data_st |>
#   mutate(age = floor(mean_max_depth * 2) / 2) |>
#   group_by(age) |>
#   arrange(sens) |>
#   mutate(y = row_number() - 1) |>
#   ggplot() +
#   geom_text(aes(x = age, y = y, label = sens), size = 3, family = base_font)

bind_rows(summary_data_st, summary_data_ie) |>
  ggplot() +
  geom_histogram(
    aes(x = mean_max_depth),
    fill = plt[2],
    color = "white",
    binwidth = .5
  ) +
  xlab("Mean maximum age (ka BP)") +
  ylab("Number of concepts with\na reconstructible root") +
  facet_wrap(~family, scales = "free")
ggsave(
  here("output/figs/hist_age_concepts.pdf"),
  width = width,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/hist_age_concepts.pdf"))

# summary_data_st |>
#   ggplot() +
#   geom_point(aes(x = mean_max_depth, y = n_traits))

prob_first_split_mcc <- read_csv(here(
  "output/results/prob_first_split_mcc.csv"
)) |>
  rename(p = mean_mcc_prob) |>
  mutate(type = "MCC")
prob_first_split_hipstr <- read_csv(here(
  "output/results/prob_first_split_hipstr.csv"
)) |>
  rename(p = mean_hipstr_prob) |>
  mutate(type = "HIPSTR")

marginal_probability_first_split_hdi <- read_csv(here(
  "output/results/marginal_prob_first_split_hdi.csv"
))

marginal_probability_first_split_hdi |>
  ggplot() +
  geom_ribbon(
    aes(x = tree_age, ymin = prob_inf, ymax = prob_sup),
    fill = plt[3],
    alpha = .25
  ) +
  geom_segment(
    aes(
      x = -Inf,
      xend = min(
        filter(marginal_probability_first_split_hdi, prob_mean <= .5)$tree_age
      ),
      y = .5,
      yend = .5
    ),
    linetype = "dashed",
    color = "grey40",
    linewidth = .35
  ) +
  geom_segment(
    aes(
      x = min(
        filter(marginal_probability_first_split_hdi, prob_mean <= .5)$tree_age
      ),
      xend = min(
        filter(marginal_probability_first_split_hdi, prob_mean <= .5)$tree_age
      ),
      y = .5,
      yend = -Inf
    ),
    linetype = "dashed",
    color = "grey40",
    linewidth = .35
  ) +
  geom_pointpath(
    aes(x = tree_age, y = prob_mean),
    color = plt[2],
    linewidth = .85,
    stroke = .1
  ) +
  scale_x_continuous(breaks = seq(0, 17, 1)) +
  xlab("Age (ka BP)") +
  ylab("Mean probability of correctly\ninferring the first split")
ggsave(
  here("output/figs/marginal_probability_first_split_hdi.pdf"),
  width = width * .8,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/marginal_probability_first_split_hdi.pdf"))

count_true_to_cs_data_long <- read_csv(here(
  "output/results/count_true_to_cs_data_long.csv"
)) |>
  mutate(summary_type = "CS") |>
  rename(type = value) |>
  mutate(
    type = case_when(
      type == "0" ~ "Discordant",
      type == "1" ~ "Concordant",
      type == "2" ~ "Reconcilable"
    )
  )
count_true_to_mcc_data_long <- read_csv(here(
  "output/results/count_true_to_mcc_data_long.csv"
)) |>
  mutate(summary_type = "MCC") |>
  rename(type = exist) |>
  mutate(
    type = case_when(
      type == "0" ~ "Discordant",
      type == "1" ~ "Concordant"
    )
  )

bind_rows(count_true_to_cs_data_long, count_true_to_mcc_data_long) |>
  mutate(
    type = factor(type, levels = c("Discordant", "Reconcilable", "Concordant"))
  ) |>
  mutate(
    p = mean_n / sum(mean_n, na.rm = TRUE),
    .by = c(summary_type, n_trait)
  ) |>
  ggplot(aes(x = as.factor(n_trait), y = p, fill = as.factor(type))) +
  geom_bar(
    stat = "identity",
    position = "fill",
    linewidth = .15
    # width = .75
  ) +
  geom_hline(
    yintercept = .5,
    linetype = "dashed",
    linewidth = .5,
    color = "white"
  ) +
  scale_fill_highcontrast(reverse = TRUE) +
  labs(
    x = "Number of traits",
    y = "Average proportion\nof nodes",
    fill = ""
  ) +
  facet_wrap(~summary_type) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

ggsave(
  here("output/figs/number_of_traits_influence.pdf"),
  width = width * .8,
  height = height,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/number_of_traits_influence.pdf"))

# 6. plot of one plausible rake node and one not plausible rake node in the consensus
# import one consensus and one true tree of age 10 simulation 1
tree_cs <- read.tree(here(
  'data/simulated-2025-07-28/beast-data-sim-7/beast-data-sim-7-15/consensus-15.tree'
))
labs_lex <- sort(tree_cs$tip.label)
map_lex <- setNames(seq_along(labs_lex), labs_lex)
tree_cs$tip.label <- paste0("T", map_lex[tree_cs$tip.label])
tree_true <- read.tree(here(
  'data/simulated-2025-07-28/beast-data-sim-7/beast-data-sim-7-15/tree-sim-7-15.tree'
))
labs_lex <- sort(tree_true$tip.label)
map_lex <- setNames(seq_along(labs_lex), labs_lex)
tree_true$tip.label <- paste0("T", map_lex[tree_true$tip.label])

# analyse a concordant rake node in the consensus tree
node_plausible <- 53 # in the true tree
descendant_plausible <- Descendants(tree_true, node_plausible)[[1]]
tip_plausible <- tree_true$tip.label[descendant_plausible]
mrca_plausible <- getMRCA(tree_cs, tip_plausible) # mrca in the cs tree

# analyse a discordant rake node in the consensus tree
node_not_plausible <- 71 # node in the true tree
descendant_not_plausible <- Descendants(tree_true, node_not_plausible)[[1]]
tip_not_plausible <- tree_true$tip.label[descendant_not_plausible]
mrca_not_plausible <- getMRCA(tree_cs, tip_not_plausible) # mrca in the cs tree

# analyse a reconciliable rake node in the consensus tree
node_rec <- 99 # in the true tree
descendant_rec <- Descendants(tree_true, node_rec)[[1]]
tip_rec <- tree_true$tip.label[descendant_rec]
mrca_rec <- getMRCA(tree_cs, tip_rec) # mrca in the cs tree

# colors in the true tree
# tip colors
tip_colors_tt <- tirep("black", length(tree_true$tip.label))
# concordent in blue
tip_colors_tt[descendant_plausible] <- plt2[1]
# discordant node in red
tip_colors_tt[descendant_not_plausible] <- plt2[3]
# reconciliable in yellow
tip_colors_tt[descendant_rec] <- plt2[2]

tb <- tibble(tip.label = tree_true$tip.label) |>
  mutate(
    type = case_when(
      tip.label %in% tip_plausible ~ "Concordant",
      tip.label %in% tip_not_plausible ~ "Discordant",
      tip.label %in% tip_rec ~ "Reconcilable",
      .default = NA
    ) |>
      factor(
        levels = c("Discordant", "Reconcilable", "Concordant")
      )
  )

# colors in the consensus tree
tip_colors_cs <- rep("black", length(tree_cs$tip.label))

# identify tip position in the consensus tree
# concordant
tip_positions_cs <- match(tip_plausible, tree_cs$tip.label)
tip_colors_cs[tip_positions_cs] <- plt2[1]
# discordant
tip_positions_cs_not_plausible <- match(tip_not_plausible, tree_cs$tip.label)
tip_colors_cs[tip_positions_cs_not_plausible] <- plt2[3]
# reconciliable
tip_positions_cs_rec <- match(tip_rec, tree_cs$tip.label)
tip_colors_cs[tip_positions_cs_rec] <- plt2[2]

tree_true_plot <- ggtree(tree_true, linewidth = .25, ladderize = FALSE) %<+%
  tb +
  geom_tiplab(
    aes(label = label, color = type),
    size = 8 / .pt,
    # color = tip_colors_tt,
    family = base_font
  ) +
  coord_cartesian(clip = "off")
tree_cs_plot <- ggtree(tree_cs, linewidth = .25, ladderize = FALSE) %<+%
  tb +
  geom_tiplab(
    aes(label = label, color = type, fill = type),
    size = 8 / .pt,
    hjust = 1,
    # key_glyph = draw_key_rect,
    family = base_font
  ) +
  geom_tippoint(
    aes(color = type, fill = type),
    color = NA,
    shape = 22,
    alpha = 0
  ) +
  scale_x_reverse() +
  coord_cartesian(clip = "off")
(tree_true_plot +
  hexpand(.05) +
  guides(color = "none") +
  ggtitle("True tree") +
  tree_cs_plot +
  guides(
    fill = guide_legend(override.aes = list(size = 5, alpha = 1))
  ) +
  ggtitle("Consensus tree") +
  hexpand(.05, direction = 1) &
  theme_void(base_family = base_font, base_size = 9) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 9, hjust = 0.5)
  ) &
  scale_color_highcontrast(reverse = TRUE, na.value = "black", guide = "none") &
  scale_fill_highcontrast(reverse = TRUE, na.translate = FALSE) &
  labs(color = "", fill = "")) +
  # guide_area() +
  plot_layout(guides = 'collect')
ggsave(
  here("output/figs/plausible_node.pdf"),
  width = width,
  height = height / 1.5,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/plausible_node.pdf"))
