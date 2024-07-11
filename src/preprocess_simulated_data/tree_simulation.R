### DON'T RUN THIS PROGRAMM
library(here)
library(TreeSim)

path <- here("data/simulated/beast-data-sim-")

# Parameter initialization
N = 50
div = 0.245
r = 0.293
s = 0.226

lambda <- div / (1 - r)
mu <- (r * div) / (1 - r)
psi <- (s / (1 - s)) * (r * div / (1 - r))
t = 1

# Tree simulation
tree <- sim.bd.taxa.age(
  n = N,
  numbsim = 1,
  lambda = lambda,
  mu = mu,
  frac = 1,
  age = t,
  mrca = TRUE
)

tree <- tree[[1]]


path <- paste(path, t, ".tree", sep = "")
#write.nexus(tree, file = pathF)

