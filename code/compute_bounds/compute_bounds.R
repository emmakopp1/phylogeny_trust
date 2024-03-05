# rm(list=ls())
library(jsonlite)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)
library(here)

here("compute_bounds/functions.R")

# Config
path_compute_bounds_function <- here("code/compute_bounds/functions.R")
path_to_config <- here("code/compute_bounds/config.json")
suppressWarnings({
  source(path_compute_bounds_function)
})

# Data
data_st <- get_tree_par_fun(path_to_config, "sino-tibetan", n_tree = 1)
data_iecor <- get_tree_par_fun(path_to_config, "iecor", n_tree = 1)
data_bantu <- get_tree_par_fun(path_to_config, "bantu", n_tree = 1)

# Topology function
f_topology_st <- data_st$f_topology
f_topology_bantu <- data_bantu$f_topology
f_topology_iecor <- data_iecor$f_topology

# Root function
f_root_st <- data_st$f_root
f_root_bantu <- data_bantu$f_root
f_root_iecor <- data_iecor$f_root

# Bantu subset
data_bantu_sub <- get_tree_par_fun(path_to_config, "bantu_subsample")
data_bantu_sub2 <- get_tree_par_fun(path_to_config, "bantu_subsample_2")

f_topology_bantu_sub <- data_bantu_sub$f_topology
f_topology_bantu_sub2 <- data_bantu_sub2$f_topology

f_root_bantu_sub <- data_bantu_sub$f_root
f_root_bantu_sub2 <- data_bantu_sub2$f_root

# Valeur issue de Tracer
data_st$param$t
data_iecor$param$t
data_bantu$param$t

library(tidyverse)
library(kableExtra)
dt_params <- list("Sino-Tibetan" = data_st$param, "Bantu" = data_bantu$param, "Indo-European" = data_iecor$param)
bounds_tb <- dt_params |>
  map_df(
    ~ tibble(
      t = .x$t, k = .x$k, N = .x$n,
      DT = compute_upper_bound_topology(t, k, .x$Q, N),
      DR = compute_upper_bound_root(t, .x$Q, N),
      inf_topo = find_t_value(k, .x$Q, N),
      inf_root = find_t_value_root(.x$Q, N)
    )
  ) |>
  mutate(
    family = names(dt_params),
    "Substitution model" = "CTMC",
    "Clock model" = "strict",
    "Tree model" = "BD",
    .before = t
  )

bounds_tb |>
  mutate(across(where(is.numeric), ~ as.character(round(.x, 2)))) |>
  pivot_longer(-family, names_to = " ") |>
  mutate(value = str_replace(value, "(\\D{2,})", "{\\1}")) |>
  mutate(family = paste0("{", family, "}")) |>
  mutate(` ` = str_replace(` `, "^(\\S+)$", "$\\1$")) |>
  mutate(` ` = str_replace(` `, "D(?=[TR])", "\\\\Delta^")) |>
  mutate(` ` = str_replace(` `, "inf_", "\\\\inf_{\\\\text{x}}\\\\{\\\\Delta^")) |>
  mutate(` ` = str_replace(` `, "topo", "T(x) = 1\\\\}")) |>
  mutate(` ` = str_replace(` `, "root", "R(x) = 1\\\\}")) |>
  pivot_wider(names_from = family) |>
  kbl(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE, align = c("l", "S", "S", "S")) |>
  write_lines("tab_upperbound.tex")

# Bounds values
# Sino tibetan
compute_upper_bound_topology(
  data_st$param$t,
  data_st$param$k,
  data_st$param$Q,
  data_st$param$n
)

# Bantu
compute_upper_bound_topology(
  data_bantu$param$t,
  data_bantu$param$k,
  data_bantu$param$Q,
  data_bantu$param$n
)

compute_upper_bound_topology(
  data_iecor$param$t,
  data_iecor$param$k,
  data_iecor$param$Q,
  data_iecor$param$n
)

# You can do the same by using the function compute_upper_bound_root


# ----------------- Graph : dataset comparison ---------------------
# Grayscale color palette
gray_palette <- c("#000000", "#000000", "#000000")

# Data generation
t_values <- seq(0, 20, length.out = 100)

data <- tibble(
  t = rep(t_values, 3),
  "Delta_T" = c(f_topology_st(t_values), f_topology_bantu(t_values), f_topology_iecor(t_values)),
  "Delta_R" = c(f_root_st(t_values), f_root_bantu(t_values), f_root_iecor(t_values)),
  family = rep(c("Sino-tibetan", "Bantu", "Indo-european"), each = 100)
)


# Plotting

library(ggthemes)
library(extrafont)
library(knitr)
theme_set(theme_minimal(base_family = "Noto Sans"))
wdt <- 18
hgt <- wdt * .6

fig_bounds <- data |>
  pivot_longer(-c(t, family)) |>
  mutate(name = str_remove(name, "Delta_")) |>
  mutate(name = factor(name, levels = c("T", "R"))) |> 
  ggplot(aes(x = t, y = value, linetype = family, color = family)) +
  geom_line() +
  xlab(expression(italic(t))) +
  ylab(expression(Delta)) +
  scale_color_few("Dark") +
  # scale_color_ptol("Vibrant") +
  facet_wrap(~name, scales = "free", labeller = label_bquote(Delta^.(as.character(name)))) +
  theme(legend.position = "bottom")

# fig_bounds <- (fig_topo_bounds + fig_root_bounds) & plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(here("figs/fig_bounds.pdf"), fig_bounds, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("figs/fig_bounds.pdf"))

# ----------------- Graph : Bantu datasets comparison ---------------------
# Data generation
data_topology <- data.frame(
  t = t_values,
  delta_T_bantu = f_topology_bantu(t_values),
  delta_T_bantu_sub = f_topology_bantu_sub(t_values),
  delta_T_bantu_sub2 = f_topology_bantu_sub2(t_values)
)

data_root <- data.frame(
  t = t_values,
  delta_R_bantu = f_root_bantu(t_values),
  delta_R_bantu_sub = f_root_bantu_sub(t_values),
  delta_R_bantu_sub2 = f_root_bantu_sub2(t_values)
)


# Topology Plot
p_topology <- ggplot(data_topology, aes(x = t)) +
  geom_line(aes(y = delta_T_bantu, linetype = "Bantu N=423"), color = gray_palette[1]) +
  geom_line(aes(y = delta_T_bantu_sub, linetype = "Bantu N=106"), color = gray_palette[2]) +
  geom_line(aes(y = delta_T_bantu_sub2, linetype = "Bantu N=51"), color = gray_palette[3]) +
  ylab(expression(Delta^T)) +
  xlab("t") +
  scale_linetype_manual(values = c("Bantu N=423" = 3, "Bantu N=106" = 2, "Bantu N=51" = 1)) +
  theme_minimal()


# Root Plot
p_root <- ggplot(data_root, aes(x = t)) +
  geom_line(aes(y = delta_R_bantu, linetype = "Bantu N=423"), color = gray_palette[1]) +
  geom_line(aes(y = delta_R_bantu_sub, linetype = "Bantu N=106"), color = gray_palette[2]) +
  geom_line(aes(y = delta_R_bantu_sub2, linetype = "Bantu N=51"), color = gray_palette[3]) +
  ylab(expression(Delta^R)) +
  xlab("t") +
  scale_linetype_manual(values = c("Bantu N=423" = 3, "Bantu N=106" = 2, "Bantu N=51" = 1)) +
  theme_minimal()



# ----------------- Graph : Choose your figure  ---------------------

# For upper bounds depending on the dataset : active p1 and p2
# For upper bounds depending on the size sample of the bantu dataset : active p_root and p_topology

par(mfrow = c(1, 2))
plot_combined <- plot_grid(
  # p1+theme(legend.position = "none"),
  # p2+theme(legend.position = "none"),
  p_root + theme(legend.position = "none"),
  p_topology + theme(legend.position = "none"),
  rel_heights = c(1, 1),
  align = "v"
)

# Print the combined plots
print(plot_combined)


# ---------------- Infima --------------------------------

# For topology
# Sino-tibetain
inf_topology_st <- find_t_value(data_st$param$k, data_st$param$Q, data_st$param$n)
inf_topology_st

# Bantu
inf_topology_bantu <- find_t_value(data_bantu$param$k, data_bantu$param$Q, data_bantu$param$n)
inf_topology_bantu

# Iecor
inf_topology_iecor <- find_t_value(data_iecor$param$k, data_iecor$param$Q, data_iecor$param$n)
inf_topology_iecor


# For root
# Sino-tibetain
inf_root_st <- find_t_value_root(data_st$param$Q, data_st$param$n)
inf_root_st

# Bantu
inf_root_bantu <- find_t_value_root(data_bantu$param$Q, data_bantu$param$n)
inf_root_bantu

# Iecor
inf_root_iecor <- find_t_value_root(data_iecor$param$Q, data_iecor$param$n)
inf_root_iecor
