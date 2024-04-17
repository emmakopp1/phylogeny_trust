library(here)
library(tidyverse)
library(jsonlite)

library(tracerer)

list.dirs(here("data/real/"), full.names = TRUE, recursive = FALSE) |> 
  map(~
      # list.files(.x, "\\.log")
      list.files(.x, "\\.trees")
  )

logfile <- here("data/real/sino-tibet-ctmc-strict-bd-fossilsRemoved/sino-tibetan-ctmc-strict-bd.log")
# logfile <- here("data/real/IECoR-ctmc-strict-fbd/IECoR2-chr_1695819208636.log")
beast_log_full <- parse_beast_tracelog_file(logfile)
beast_log <- remove_burn_ins(beast_log_full, burn_in_fraction = 0.2)
beast_log |> 
  as_tibble() |> 
  select(starts_with("freqParameter"), TreeHeight.t.tree) |> 
  summarise(across(everything(), ~ mean(.x))) |> 
  rename_all(str_replace, pattern = "freq.+(\\d)", replacement = "pi\\1") |> 
  rename(pi0 = pi1, pi1 = pi2, t_R = TreeHeight.t.tree) |> 
  mutate(nTrees = max(beast_log$Sample)) |> 
  mutate(family = "Sino-Tibetan") |> 
  relocate(family, .before = pi0)

cfg <- here("src/compute_bounds/config.json")
dt <- fromJSON(cfg) |> 
  as_tibble() |> 
  mutate(parameter = c("path", "pi0", "pi1", "k", "path_cognates", "t_conv")) |> 
  pivot_longer(-parameter, names_to = "family") |>
  pivot_wider(names_from = parameter) |> 
  relocate(path_cognates, .after = path) |> 
  mutate(across(4:7, as.numeric))


# This function computes the substitution model transition matrix
compute_chain <- function(pi0, pi1) {
  return(1 / (pi0^2 + pi1^2) * as.matrix(rbind(c(-pi0, pi0), c(pi1, -pi1))))
}

Q <- compute_chain(dt$pi0[1], dt$pi1[1])

get_tree_param <- function(tree, pi0, pi1, k, t) {
  Q <- compute_chain(pi0, pi1)
  n <- mean(sapply(tree, function(arbre) length(arbre$tip.label)))
  return(list(n = n, t = t, k = k, Q = Q, pi0 = pi0, pi1 = pi1))
}


compute_q <- function(Q) {
  d <- as.numeric(dim(Q)[1])
  return(sum(apply(Q + 100 * diag(d), 2, min)))
}

compute_q(Q)
# This function computes the upper bound of the probability of the exact topology
# reconstruction
compute_upper_bound_topology <- function(t, k, Q, n) {
  d <- as.numeric(dim(Q)[1])
  q <- sum(apply(Q + 100 * diag(d), 2, min))
  return(k * n * exp(-q * t))
}

pi0 <- dt$pi0[1]
pi1 <- dt$pi1[1]
k <- 3785
N <- 46
t <- 8.867

Q <- 1 / (pi0^2 + pi1^2) * as.matrix(rbind(c(-pi0, pi0), c(pi1, -pi1)))
d <- as.numeric(dim(Q)[1])
q <- sum(apply(Q + 100 * diag(d), 2, min))
k * N * exp(-q * t)

t_values <- seq(0, 20, length.out = 101)

(k * N * exp(-q * t_values))
# This function computes the upper bound of the probability of the exact root
# reconstruction
compute_upper_bound_root <- function(t, Q, n, pi0, pi1) {
  m <- max(pi0,pi1)
  q <- compute_q(Q)
  return(m + n * exp(-q * t))
}

m <- max(pi0,pi1)
q <- sum(apply(Q + 100 * diag(d), 2, min))
m + N * exp(-q * t)


# -------------- Infima of both bounds -----------------


find_t_value <- function(k, Q, n, tolerance = 1e-6, max_iter = 1000) {
  objective_function <- function(t) {
    return(compute_upper_bound_topology(t, k, Q, n) - 1)
  }
  
  result <- uniroot(objective_function, interval = c(0, 20), tol = tolerance, maxiter = max_iter)
  
  return(result$root)
}

f <- function(x) {(k * N * exp(-q * x)) - 1}
uniroot(function(x) {(k * N * exp(-q * x)) - 1}, interval = c(0, 20), tol = 1e-6, maxiter = 1000)$root


find_t_value_root <- function(Q, n, pi0, pi1, tolerance = 1e-6, max_iter = 1000) {
  objective_function <- function(t) {
    return(compute_upper_bound_root(t, Q, n, pi0, pi1) - 1)
  }
  
  result <- uniroot(objective_function, interval = c(0, 20), tol = tolerance, maxiter = max_iter)
  return(result$root)
}


