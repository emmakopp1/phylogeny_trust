library(here)
source(here("src/utils.R"))
source(here("src/init.R"))

# Functions
get_deepest_node = function(tree,N){
  ages = as.data.frame(round(dateNodes(tree),2))
  colnames(ages)=c('ages')
  # We don't select the root
  return(order(ages$ages,decreasing = T)[2:(N+1)])
}


# Initialisation of the node and the matrix of results
path_true = "~/Documents/phylogeny_trust_local/reconstruction/simulated_bd_trees/beast-data-sim-5/tree-sim-5.tree"
tree_true = read.tree(path_true)
nodes = get_deepest_node(tree_true,10)
res = matrix(NA,nrow = 17,ncol=10)
colnames(res)=as.character(nodes)

# For all age 
common_path = "/Users/kopp/Documents/phylogeny_trust_local/reconstruction/simulated_bd_trees/"

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

# PLot 
c.pal = viridis(17) 

res.t = as.data.frame(t(res))
colnames(res.t) = as.character(seq(1,17,1))

pdf(here("src/figs/boxplot.pdf"))
boxplot(res.t,col=c.pal[1:17],xlab=c('age of the phylogeny'))
dev.off()



