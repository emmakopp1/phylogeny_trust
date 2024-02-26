rm(list=ls())
library(jsonlite)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)


# Config 
path_compute_bounds_function = "~/Documents/phylogeny_trust/code/compute_bounds/functions.R"
path_to_config = "/Users/kopp/Documents/phylogeny_trust/code/compute_bounds/config.json"
suppressWarnings({source(path_compute_bounds_function)})

# Data 
data_st = get_tree_par_fun(path_to_config,'sino-tibetan',n_tree = 1)
data_iecor = get_tree_par_fun(path_to_config,'iecor',n_tree = 1)
data_bantu = get_tree_par_fun(path_to_config,'bantu',n_tree = 1)


# Topology function 
f_topology_st = data_st$f_topology
f_topology_bantu = data_bantu$f_topology
f_topology_iecor = data_iecor$f_topology

# Root function 
f_root_st = data_st$f_root
f_root_bantu = data_bantu$f_root
f_root_iecor = data_iecor$f_root

# Bantu subset
data_bantu_sub = get_tree_par_fun(path_to_config,'bantu_subsample')
data_bantu_sub2 = get_tree_par_fun(path_to_config,'bantu_subsample_2')

f_topology_bantu_sub = data_bantu_sub$f_topology
f_topology_bantu_sub2 = data_bantu_sub2$f_topology

f_root_bantu_sub = data_bantu_sub$f_root
f_root_bantu_sub2 = data_bantu_sub2$f_root

# Valeur issue de Tracer
data_st$param$t
data_iecor$param$t 
data_bantu$param$t 



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

data <- data.frame(
  t = rep(t_values, 3),
  delta_T = c(f_topology_st(t_values), f_topology_bantu(t_values), f_topology_iecor(t_values)),
  delta_R = c(f_root_st(t_values), f_root_bantu(t_values), f_root_iecor(t_values)),
  Group = rep(c('Sino-tibetan', 'Bantu', 'Indo-european'), each = 100)
)


# Plotting

# Topology upper bound
p1 = ggplot(data, aes(x = t, y = delta_T, color = Group, linetype = Group)) +
  geom_line() +
  scale_color_manual(values = gray_palette) +
  ylab(expression(Delta^T)) +
  xlab("t") +
  theme_minimal()

# Root upper bound
p2 = ggplot(data, aes(x = t, y = delta_R, color = Group, linetype = Group)) +
  geom_line() +
  scale_color_manual(values = gray_palette) +
  ylab(expression(Delta^R)) +
  xlab("t") +
  theme_minimal()



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

par(mfrow=c(1,2))
plot_combined <- plot_grid(
  #p1+theme(legend.position = "none"),
  #p2+theme(legend.position = "none"),
  p_root + theme(legend.position =  "none"),
  p_topology + theme(legend.position =  "none"),
  rel_heights = c(1, 1),
  align = "v"
)

# Print the combined plots
print(plot_combined)


# ---------------- Infima --------------------------------

# For topology
# Sino-tibetain 
inf_topology_st = find_t_value(data_st$param$k, data_st$param$Q, data_st$param$n)
inf_topology_st

# Bantu
inf_topology_bantu = find_t_value(data_bantu$param$k, data_bantu$param$Q, data_bantu$param$n)
inf_topology_bantu

# Iecor 
inf_topology_iecor = find_t_value(data_iecor$param$k, data_iecor$param$Q, data_iecor$param$n)
inf_topology_iecor


# For root 
# Sino-tibetain 
inf_root_st = find_t_value_root(data_st$param$Q, data_st$param$n)
inf_root_st

# Bantu
inf_root_bantu = find_t_value_root(data_bantu$param$Q, data_bantu$param$n)
inf_root_bantu

# Iecor 
inf_root_iecor = find_t_value_root(data_iecor$param$Q, data_iecor$param$n)
inf_root_iecor
