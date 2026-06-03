library(tidyverse)
library(here)
library(ggh4x)
library(legendry)
library(khroma)
library(knitr)

base_font <- "Noto Sans Condensed"

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
prop_shared_tip_pair <- readRDS(here(
  "output/results/shared_cognate_tip_pair.pdf"
)) |>
  rowid_to_column("age")

plt <- color("vibrant")(3)
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
