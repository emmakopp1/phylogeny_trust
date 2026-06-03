library(tidyverse)
library(here)
library(ggh4x)
library(legendry)
library(khroma)
library(knitr)

base_font <- "Noto Sans Condensed"
plt <- color("vibrant")(3)

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

shared_cognates <- readRDS(here("output/results/shared_cognates.rds"))


shared_cognates |>
  ggplot() +
  geom_pointpath(
    aes(x = tree_age, y = value),
    color = plt[2],
    linewidth = 1,
    stroke = .1
  ) +
  geom_errorbar(
    aes(x = tree_age, ymin = inf, ymax = sup),
    color = plt[2],
    width = .2
  ) +
  ylab("Proportion of shared cognates") +
  xlab("Age of the most recent common ancestor (ka)") +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_x_continuous(breaks = seq(0, 17, 1)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = base_font) +
  theme(
    aspect.ratio = .618,
    panel.grid.minor = element_blank()
  )
ggsave(
  here("output/figs/shared_cognate_no_homoplasie.pdf"),
  width = 12,
  height = 12,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/shared_cognate_no_homoplasie.pdf"))

prop_shared_tip_pair <- readRDS(here(
  "output/results/shared_cognate_tip_pair.pdf"
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
  theme_minimal(base_family = base_font) +
  theme(
    aspect.ratio = .618,
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.justification = c(1, 1),
    legend.text = element_text(size = 12),
  )

ggsave(
  here("output/figs/shared_cognate_tip_pair.pdf"),
  width = 12,
  height = 12,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/shared_cognate_tip_pair.pdf"))

count_true_to_mcc <- read_rds(here("output/results/count_true_to_mcc.csv")) |>
  select(age, value = exist, mean_n = n_mean) |>
  mutate(type = "MCC", value = as.character(value))
count_true_to_cs <- readRDS(here("output/results/count_true_to_cs.rds")) |>
  mutate(type = "Majority-rule consensus", value = as.character(value))

bind_rows(count_true_to_mcc, count_true_to_cs) |>
  mutate(
    value = case_when(
      value == "1" ~ "Present",
      value == "2" ~ "Uncertain",
      value == "0" ~ "Absent"
    )
  ) |>
  mutate(value = factor(value, levels = c("Absent", "Uncertain", "Present"))) |>
  # group_by(type, value) |>
  # mutate(p = mean_n)
  # ungroup() |>
  ggplot(aes(x = factor(age), y = mean_n, fill = value)) +
  geom_col(position = "fill") +
  labs(
    x = "Age (ka BP)",
    y = "Average proportion of nodes",
    fill = ""
  ) +
  scale_fill_highcontrast(reverse = TRUE) +
  coord_cartesian(clip = "off", expand = FALSE) +
  facet_wrap(~type) +
  theme_minimal() +
  theme_minimal(base_family = base_font) +
  theme(
    aspect.ratio = .618,
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12),
  )
ggsave(
  here("output/figs/barplot_prop_true_to_resume.pdf"),
  width = 12 * 1.35,
  height = 12,
  units = "cm",
  device = cairo_pdf
)
plot_crop(here("output/figs/barplot_prop_true_to_resume.pdf"))
