library(here)
source(here("src/utils.R"))


# Functions
get_deepest_node = function(tree,N){
  ages = as.data.frame(round(dateNodes(tree),2))
  colnames(ages)=c('ages')
  # We don't select the root
  return(order(ages$ages,decreasing = T)[2:(N+1)])
}

remove_burnin = function(trees,burnin_rate){
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}

myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}


# Initialisation of the node and the matrix of results
path_true = here("data/simulated/beast-data-sim-5/tree-sim-5.tree")
tree_true = read.tree(path_true)
nodes = get_deepest_node(tree_true,10)
res = matrix(NA,nrow = 17,ncol=10)
colnames(res)=as.character(nodes)

# For all age 
common_path = here("data/simulated/")

for (i in 1:17){
  # Read path
  path_true = sprintf("%sbeast-data-sim-%i/tree-sim-%i.tree", common_path, i, i)
  path_phylo = sprintf("%sbeast-data-sim-%i/ctmc-strict-bd-%i.trees", common_path, i, i)
  
  # Load trees
  tree_phylo = read.nexus(path_phylo)
  tree_true = read.tree(path_true)
  
  # Remove burnin of the sample
  tree_phylo = remove_burnin(tree_phylo,0.2)
  
  # For all nodes compute the monophylecy of a group from the simulation with truth
  for (s in 1:10){
    node = nodes[s]
    childrens = Descendants(tree_true, node, type = c("tips"))[[1]]
    res[i,s] = mean(sapply(tree_phylo, function(t) is.monophyletic(t, tree_true$tip.label[childrens])))
  }
}

library(ape)
trees = read.nexus("/Users/kopp/Documents/phylogeny_trust/data/simulated/beast-data-sim-7/ctmc-strict-bd-7.trees")
trees = remove_burn_in(trees,0.9)
consensus_tree = myconsensus(trees)
write.tree(consensus_tree,"/Users/kopp/Documents/phylogeny_trust/data/simulated/beast-data-sim-7/consensus-7.tree")



res.t = as.data.frame(t(res))
colnames(res.t) = as.character(seq(1,17,1))
node_probs_tb <- as_tibble(res.t) |> 
  pivot_longer(everything(), names_to = "age", values_to = "p") |> 
  mutate(age = as.integer(age)) |> 
  arrange(age, p)
write_csv(node_probs_tb, here("output/results/node_probs_tb.csv"))
