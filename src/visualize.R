library(here)
library(tidyverse)
library(ggthemes)
library(knitr)
library(kableExtra)

theme_set(
  theme_minimal(base_family = "Noto Sans SemiCondensed", base_size = 10) +
    theme(axis.title.y = element_text(face = "italic"),
          legend.position = "bottom", 
          aspect.ratio = .618)
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
  mutate(parameter = str_replace(parameter, "^(\\S+)$", "$\\1$")) |>
  mutate(parameter = str_replace(parameter, "D(?=[TR])", "\\\\Delta^")) |>
  mutate(parameter = str_replace(parameter, "inf_", "\\\\inf_{\\\\text{x}}\\\\{\\\\Delta^")) |>
  mutate(parameter = str_replace(parameter, "topo", "T(x) = 1\\\\}")) |>
  mutate(parameter = str_replace(parameter, "root", "R(x) = 1\\\\}")) |>
  pivot_wider(names_from = family) |>
  mutate(parameter = str_replace(parameter, "t\\$", "t$ (root age, ka BP)")) |> 
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
  ylab(expression(Delta)) +
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
  ylab(expression(Delta)) +
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
  ggplot(aes(x = factor(age), y = p, group = factor(age)))+ 
  geom_boxplot(fill = "gray90", outliers = FALSE) +
  geom_point(position = position_jitter(seed = 0, width = .3), size = 1.5, alpha = 1, color = few_pal("Dark")(2)[1], shape = 1) +
  xlab("age (ka BP)")
ggsave(here("output/figs/fig_nodeprobs.pdf"), fig_nodeprobs, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_nodeprobs.pdf"))
