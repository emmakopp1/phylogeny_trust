library(here)
library(tidyverse)
library(ggthemes)
library(knitr)
library(kableExtra)

font <- "Noto Sans SemiCondensed"
theme_set(
  theme_minimal(base_family = font, base_size = 10) +
    theme(
      legend.position = "bottom",
      aspect.ratio = .618
    )
)
wdt <- 14
hgt <- wdt * .7


bounds_tb <- read_csv(here("output/results/bounds_tb.csv"))
bounds_tb |>
  mutate(across(where(is.numeric), ~ as.character(round(.x, 2)))) |>
  mutate(DT = ifelse(DT >= 1, "\\geq 1", DT)) |>
  arrange(family) |>
  pivot_longer(-family, names_to = "parameter") |>
  mutate(value = str_replace(value, "^(\\D{2,})$", "{\\1}")) |>
  mutate(family = paste0("{", family, "}")) |>
  mutate(parameter = str_replace(parameter, "pi", "\\\\pi_")) |>
  mutate(parameter = str_replace(parameter, "^(\\S+)$", "$\\1$")) |>
  mutate(parameter = str_replace(parameter, "\\$D([TS])\\$", "upper bound of $\\\\Delta^\\1(t_R)$")) |>
  mutate(parameter = str_replace(parameter, "inf_", "\\\\inf_{\\\\text{t}}\\\\{\\\\Delta^")) |>
  mutate(parameter = str_replace(parameter, "topo", "T(t) = 1\\\\}")) |>
  mutate(parameter = str_replace(parameter, "root", "R(t) = 1\\\\}")) |>
  pivot_wider(names_from = family) |>
  mutate(parameter = str_replace(parameter, "t\\$", "t_R$ (root age, ka BP)")) |>
  mutate(parameter = str_replace(parameter, "k\\$", "k$ (number of traits)")) |>
  mutate(parameter = str_replace(parameter, "N\\$", "N$ (number of languages)")) |>
  mutate(parameter = ifelse(str_detect(parameter, "inf"), paste0(parameter, " (ka BP)"), parameter)) |>
  rename(" " = parameter) |>
  kbl(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE, align = c("l", "S", "S", "S", "S")) |>
  write_lines(here("output/tabs/tab_upperbound.tex"))

bounds_byt_tb <- read_csv(here("output/results/bounds_byt_tb.csv"))
fig_bounds <- bounds_byt_tb |>
  dplyr::filter(!str_detect(family, "subset")) |>
  mutate(Delta = factor(Delta, levels = c("T", "R"))) |>
  ggplot(aes(x = t, y = value, linetype = family, color = family)) +
  geom_line() +
  xlab("age (ka BP)") +
  ylab("upper bound") +
  scale_color_few("Dark") +
  facet_wrap(~Delta, scales = "free", labeller = label_bquote(Delta^italic(.(as.character(Delta)))))
ggsave(here("output/figs/fig_bounds.pdf"), fig_bounds, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_bounds.pdf"))

fig_bounds_bantu <- bounds_byt_tb |>
  dplyr::filter(str_detect(family, "Bantu")) |>
  mutate(`number of languages` = case_when(
    family == "Bantu" ~ 424,
    family == "Bantu subset" ~ 107,
    family == "Bantu subset 2" ~ 52,
  )) |>
  mutate(`number of languages` = fct_rev(factor(`number of languages`))) |>
  mutate(Delta = factor(Delta, levels = c("T", "R"))) |>
  ggplot(aes(x = t, y = value, linetype = `number of languages`, color = `number of languages`)) +
  geom_line() +
  xlab("age (ka BP)") +
  ylab("upper bound") +
  scale_color_few("Dark") +
  facet_wrap(~Delta, scales = "free", labeller = label_bquote(Delta^italic(.(as.character(Delta)))))
ggsave(here("output/figs/fig_bounds_bantu.pdf"), fig_bounds_bantu, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_bounds_bantu.pdf"))

qs_tb <- read_csv(here("output/results/qs_tb.csv"))
fig_qs <- qs_tb |>
  ggplot(aes(x = age, y = q_theo, group = family, color = family, linetype = family)) +
  geom_line() +
  xlab("age (ka BP)") +
  ylab(expression(italic(Q[s]))) +
  scale_y_continuous(n.breaks = 10) +
  scale_color_few("Dark")
ggsave(here("output/figs/fig_qs.pdf"), fig_qs, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_qs.pdf"))

node_probs_tb <- read_csv(here("output/results/node_probs_tb.csv"))
fig_nodeprobs <- node_probs_tb |>
  ggplot(aes(x = factor(age), y = p, group = factor(age))) +
  geom_boxplot(fill = "gray90", outliers = FALSE) +
  geom_point(position = position_jitter(seed = 0, width = .3), size = 1.5, alpha = 1, color = few_pal("Dark")(2)[1], shape = 1) +
  xlab("age (ka BP)") +
  ylab("proportion")
ggsave(here("output/figs/fig_nodeprobs.pdf"), fig_nodeprobs, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_nodeprobs.pdf"))

font <- "Noto Sans Condensed"

library(phangorn)
library(treeio)
library(ggtree)
plt_dark <- few_pal("Dark")(8)

t5_true <- read.newick(here("data/simulated/beast-data-sim-5/tree-sim-5.tree"))
clds <- c(71, 95, 52)
clds_color <- tibble(clade = factor(0:2), ca = clds, color = plt_dark[1:3])
clds_nds <- map_df(clds, ~ tibble(ca = .x, node = unlist(Descendants(t5_true, .x, type = "tips")))) |>
  arrange(node) |>
  mutate(label = t5_true$tip.label) |>
  left_join(clds_color)
rn <- getRoot(t5_true)

fig_t5true <- t5_true |>
  groupClade(.node = c(95, 52)) |>
  full_join(clds_nds) |>
  ggtree() +
  geom_rootedge(.25) +
  geom_tree(aes(color = group)) +
  geom_tree(data = ~ filter(.x, group == 0 & x < .5)) +
  geom_tiplab(aes(color = clade), family = font, size = 9 / .pt) +
  scale_color_few("Dark") +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(.5, 1, .5, 0, unit = "line"), legend.position = "none")
fig_t5true

t5_consensus <- read.newick(here("data/simulated/beast-data-sim-5/consensus-5.tree"))
t5_consensus$root.edge.length <- 0

fig_t5consensus <- t5_consensus |>
  groupClade(.node = c(69, 53)) |>
  full_join(select(clds_nds, -node)) |>
  ggtree(ladderize = FALSE) +
  geom_rootedge(.25) +
  geom_tree(aes(color = group)) +
  geom_tree(data = ~ slice(.x, 1:52)) +
  geom_tiplab(aes(color = clade), family = font, size = 9 / .pt) +
  # geom_nodelab(aes(label = node)) +
  scale_color_few("Dark") +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(.5, 0.5, .5, 1, unit = "line"), legend.position = "none")
fig_t5consensus

ggsave(here("output/figs/fig_t5trueconsensus.pdf"), fig_t5true + fig_t5consensus, device = cairo_pdf, width = wdt, height = hgt * 1.7, units = "cm")
plot_crop(here("output/figs/fig_t5trueconsensus.pdf"))

t9_consensus <- read.newick(here("data/simulated/beast-data-sim-9/consensus-9.tree"))
t9_consensus$root.edge.length <- 0

fig_t9consensus <- t9_consensus |>
  groupClade(.node = c(c(52, 80, 78, 87), c(82, 86), c(67, 72, 75, 84))) |>
  full_join(select(clds_nds, -node)) |>
  ggtree(ladderize = FALSE, color = "gray") +
  geom_tree(data = ~ filter(.x, group %in% 5:6 | x < 1), color = plt_dark[2]) +
  geom_tree(data = ~ filter(.x, group %in% 7:10 | x == 0), color = plt_dark[3]) +
  geom_tree(data = ~ filter(.x, group %in% 1:4 | x < 1.9), color = plt_dark[1]) +
  geom_tree(aes(color = clade)) +
  geom_tiplab(aes(color = clade), family = font, size = 9 / .pt) +
  geom_rootedge(.25) +
  # geom_nodelab(aes(label = paste(node, round(x,2)))) +
  scale_color_few("Dark") +
  coord_cartesian(clip = "off", expand = FALSE) +
  # theme_minimal() +
  theme(plot.margin = margin(.5, 1.5, .5, 1, unit = "line"), legend.position = "none")
fig_t9consensus

ggsave(here("output/figs/fig_t9trueconsensus.pdf"), fig_t5true + fig_t9consensus, device = cairo_pdf, width = wdt, height = hgt * 1.7, units = "cm")
plot_crop(here("output/figs/fig_t9trueconsensus.pdf"))

st_consensus <- read.newick(here("data/real/sino-tibet-ctmc-strict-bd-fossilsRemoved/sino-tibetan-ctmc-strict-bd-consensus.tree"))
st_consensus$root.edge.length <- 0
fig_stconsensus <- ggtree(st_consensus) +
  geom_tiplab(family = font, size = 10 / .pt) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE, xlim = c(-.25, 11)) +
  theme(plot.margin = margin(.5, 0, .5, 0, unit = "line"))
ggsave(here("output/figs/fig_stconsensus.pdf"), fig_stconsensus, device = cairo_pdf, width = wdt, height = hgt * 2, units = "cm")
plot_crop(here("output/figs/fig_stconsensus.pdf"))

ie_consensus <- read.newick(here("data/real/IECoR-ctmc-strict-fbd/IECoR2-chr_consensus.tree"))
ie_consensus$root.edge.length <- 0
fig_ieconsensus <- ggtree(ie_consensus) +
  geom_tiplab(family = font, size = 9 / .pt) +
  geom_rootedge(.25) +
  coord_cartesian(clip = "off", expand = FALSE, xlim = c(-.25, 9.35)) +
  theme(plot.margin = margin(.5, 0, .5, 0, unit = "line"))
ggsave(here("output/figs/fig_ieconsensus.pdf"), fig_ieconsensus, device = cairo_pdf, width = wdt, height = hgt * 2, units = "cm")
plot_crop(here("output/figs/fig_ieconsensus.pdf"))
