library(here)
library(tidyverse)
library(TreeTools)
library(kableExtra)
library(tracerer)
library(ape)
library(adephylo)
library(phytools)
library(castor)


# Compute the upper bound of the probability of inferring the true tree topology
compute_upperbound_DT <- function(k, N, q, t, s=0) {
  k * N * exp(-q * (t - s))
}

# Compute the upper bound of the probability of correctly inferring ancestral states version2
compute_upperbound_DS <- function(pi0, pi1, q, ages) {
  max(pi0, pi1) + sum(exp(-q*ages))
}

# Compute the time threshold beyond which the upper bound of the probability
# of inferring the true tree topology falls below 1
compute_inf_t_DT <- function(k, N, q, t, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DT(k, N, q, t) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}

# Compute the time threshold beyond which the upper bound of the probability
# of correctly inferring ancestral states falls below 1
compute_inf_t_DS <- function(pi0, pi1, N, q, t, interval = c(0, 20), tol = 1e-6, maxiter = 1000) {
  uniroot(function(t) {
    compute_upperbound_DS(pi0, pi1, N, q, t) - 1
  }, interval = interval, tol = tol, maxiter = maxiter)$root
}
