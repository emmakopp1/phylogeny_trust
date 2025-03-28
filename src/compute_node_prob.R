library(here)
library(ape)
library(phangorn)
library(Matrix)
library(matrixStats)
library(castor)
library(here)
library(purrr)
library(tidyverse)
library(phytools)
library(TreeTools)
library(adephylo)
library(paleotree)

# functions --------------------------------------------------------------------
# given a tree, compute the nodes from the deepest to youngest
get_deepest_node <- function(tree) {
  result <- list()
  # Recursive function
  preorder_recursive <- function(node_index) {
    result <<- cbind(result, node_index)

    children <- tree$edge[tree$edge[, 1] == node_index, 2]
    for (child in children) {
      preorder_recursive(child)
    }
  }
  root <- castor::find_root(tree)
  preorder_recursive(root)
  return(unlist(result))
}

# remove burnin of a phylogeny
remove_burnin <- function(phylo, burnin) {
  M <- length(phylo)
  phylo[ceiling(burnin * M):M]
}

# compute consensus tree
myconsensus <- function(trees) {
  consensus.edges(trees, consensus.tree = consensus(trees, p = .5), rooted = T)
}

# initialization of the node and the matrix of results
path_true <- here("data/simulated-2025-03-17/beast-data-sim-1/tree-sim-1.tree")
tree_true <- read.tree(path_true)
nodes <- getNodesByDepth(tree_true)[-1]
res <- matrix(NA, nrow = 17, ncol = 10)
colnames(res) <- as.character(nodes[1:10])

# For all age
common_path <- here("data/simulated-2025-03-17/")

for (i in 1:17) {
  # read paths
  path_true <- sprintf("%sbeast-data-sim-%i/tree-sim-%i.tree", common_path, i, i)
  path_phylo <- sprintf("%sbeast-data-sim-%i/ctmc-strict-bd-%i.trees", common_path, i, i)

  # load trees
  tree_phylo <- read.nexus(path_phylo)
  tree_true <- read.tree(path_true)

  # remove burnin of the sample
  tree_phylo <- remove_burnin(tree_phylo, 0.2)

  # for all nodes compute the monophylecy of a group from the simulation with truth
  for (s in 1:10) {
    node <- nodes[s]
    childrens <- Descendants(tree_true, node, type = c("tips"))[[1]]
    res[i, s] <- mean(sapply(tree_phylo, function(t) is.monophyletic(t, tree_true$tip.label[childrens])))
  }
}

View(res)


res.t <- as.data.frame(t(res))
colnames(res.t) <- as.character(seq(1, 17, 1))

boxplot(res.t, col = 'darkblue')
node_probs_tb <- as_tibble(res.t) |>
  pivot_longer(everything(), names_to = "age", values_to = "p") |>
  mutate(age = as.integer(age)) |>
  arrange(age, p)
# change path
# write_csv(node_probs_tb, here("output/results/node_probs_tb.csv"))

# from visualize
# node_probs_tb <- read_csv(here("output/results/node_probs_tb.csv"), show_col_types = FALSE)
fig_nodeprobs <- node_probs_tb |>
  ggplot(aes(x = factor(age), y = p, group = factor(age))) +
  geom_boxplot(fill = "gray90", outliers = FALSE) +
  geom_point(position = position_jitter(seed = 0, width = .3), size = 1.5, alpha = 1, color = few_pal("Dark")(2)[1], shape = 1) +
  xlab("age (ka BP)") +
  ylab("proportion")
ggsave(here("output/figs/fig_nodeprobs2025-03-17.pdf"), fig_nodeprobs, device = cairo_pdf, width = wdt, height = hgt, units = "cm")
plot_crop(here("output/figs/fig_nodeprobs2025-03-17.pdf"))