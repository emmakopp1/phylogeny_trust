rm(list=ls())
library(ape)
path_tree_5k  = "/Users/kopp/Documents/chr_paper/reconstruction/simulated_bd_trees/beast-data-sim-5/tree-sim-5.tree"
tree = read.tree(path_tree_5k)

l = seq(0.2,0.8,0.2)
for (k in l){ 
  new_tree = tree
  new_tree$edge.length=k*as.numeric(new_tree$edge.length)
  t = k*5
  write.tree(new_tree,
             file = sprintf("/Users/kopp/Documents/chr_paper/reconstruction/simulated_bd_trees/beast-data-sim-%g/tree-sim-%g.tree",t,t)
  )}
