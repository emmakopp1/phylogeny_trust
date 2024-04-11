library(ape)

path_tree_5k  = "data/simulated/beast-data-sim-5/tree-sim-5.tree"
tree = read.tree(path_tree_5k)

l = seq(0.2,0.8,0.2)
for (k in l){ 
  new_tree = tree
  new_tree$edge.length=k*as.numeric(new_tree$edge.length)
  t = k*5
  write.tree(new_tree,
             file = sprintf("data/simulated/beast-data-sim-%g/tree-sim-%g.tree",t,t)
  )}
