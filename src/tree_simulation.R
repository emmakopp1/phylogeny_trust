library(here)
library(TreeSim)
library(ape)
library(purrr)

# 1. Simulation of the initial tree
path <- here("data/simulated_temp/beast-data-sim-")
dir.create( here("data/simulated_temp"))

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

dir.create(sprintf(here("data/simulated_temp/beast-data-sim-1")))
write.nexus(tree, file = sprintf(here("data/simulated_temp/beast-data-sim-%d/tree-sim-%d.tree"), t, t))

# 2. Construction of the scaled trees
l <- seq(1, 17, 1)

for (k in l) {
  new_tree <- tree
  new_tree$edge.length <- k * as.numeric(new_tree$edge.length)
  dir.create(sprintf(here("data/simulated_temp/beast-data-sim-%d"), k))
  write.tree(new_tree,
    file = sprintf(here("data/simulated_temp/beast-data-sim-%d/tree-sim-%d.tree"), k, k)
  )
}

# Files manipulation

files = list.files(here('data/simulated_temp'),full.names = TRUE, recursive = TRUE) 
target_text = read_lines(here("data/beast-data-sim.xml")) %>% paste(collapse = "\n")

process_file <- function(file) {
  file_text <- read_lines(file) |> paste(collapse = "\n")
  str_replace(target_text, "xyz", file_text) |>
    write_lines(
      file %>%
        str_replace("/tree-sim-\\d+\\.tree$", "") %>%  # Remove the tree file part
        str_replace("beast-data-sim-\\d+", "\\0/\\0.xml")
    )
}

updated_texts <- files |>
  map_chr(~ process_file(.x)) 





