library(tidyverse)
library(ggthemes)
library(knitr)

theme_set(
  theme_minimal(base_family = "Noto Sans", base_size = 10) +
    theme(legend.position = "bottom", aspect.ratio = .618)
)
wdt <- 18
hgt <- wdt * .7


bounds_byt_tb <- write_csv(here("output/results/bounds_byt_tb.csv"))
fig_bounds <- bounds_byt_tb |>
  mutate(Delta = factor(Delta, levels = c("T", "R"))) |>
  ggplot(aes(x = t, y = value, linetype = family, color = family)) +
  geom_line() +
  xlab("age (ka BP)") +
  ylab(expression(Delta)) +
  scale_color_few("Dark") +
  facet_wrap(~Delta, scales = "free", labeller = label_bquote(Delta^.(as.character(Delta))))
ggsave(here("output/figs/fig_bounds.pdf"), fig_bounds, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_bounds.pdf"))


qs_tb <- read_csv(here("output/results/qs_tb.csv"))
fig_qs <- qs_tb |> 
  ggplot(aes(x = age, y = q_theo, group = language, color = language, linetype = language)) +
  geom_line() +
  xlab("age (ka BP)") +
  ylab("Q") +
  scale_y_continuous(n.breaks = 10) +
  scale_color_few("Dark") +
  theme(axis.title.y = element_text(face = "italic"))
ggsave(here("output/figs/fig_qs.pdf"), fig_qs, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_qs.pdf"))
