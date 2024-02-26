rm(list=ls())

library(paleotree)
library(ape)
library(phytools)
library(TreeTools)
library(dispRity)
library(ips)
library(phangorn)
library(viridisLite)
library(geiger)

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


# Tous les arbres ont la même profondeur relatives
path_true = "/Users/kopp/Documents/chr_paper/reconstruction/simulated_bd_trees/beast-data-sim-5/tree-sim-5.tree"
tree_true = read.tree(path_true)
nodes = get_deepest_node(tree_true,10)
res = matrix(NA,nrow = 17,ncol=10)
colnames(res)=as.character(nodes)

# For all age 
for (i in 1:17){
  #print(i)
  path_true = sprintf("/Users/kopp/Documents/chr_paper/reconstruction/simulated_bd_trees/beast-data-sim-%i/tree-sim-%i.tree",i,i)
  path_phylo = sprintf("/Users/kopp/Documents/chr_paper/reconstruction/simulated_bd_trees/beast-data-sim-%i/ctmc-strict-bd-%i.trees",i,i)
  
  tree_phylo = read.nexus(path_phylo)
  tree_true = read.tree(path_true)
  
  tree_phylo = remove_burnin(tree_phylo,0.2)
  
  # For all nodes 
  for (s in 1:10){
    node = nodes[s]
    #print(node)
    childrens = Descendants(tree_true, node, type = c("tips"))[[1]]
    res[i,s] = mean(sapply(tree_phylo, function(t) is.monophyletic(t, tree_true$tip.label[childrens])))
  }
}

c.pal = viridis(17) # On selectionne deux couleurs


# Save file
#png("/Users/kopp/Documents/chr_paper/report/plots/nodes_prob.png")
#matplot(seq(1,17,1),res,ylab=c(""), type = "l",pch=8,ylim=c(0,1),lty=1,lwd=1,col= c.pal,xlab=c("age of the tree"))
#legend("topright", legend = colnames(res), col=c.pal, lty = 1, cex=0.6)
#dev.off()

res.t= t(res)
res.t = as.data.frame(res.t)
colnames(res.t) = as.character(seq(1,17,1))
boxplot(res.t,col=c.pal[1:17],xlab=c('age of the phylogeny'))

# Convert data to long format
plot(tree_true)
nodelabels(frame="none",cex=0.8)
edgelabels(frame="none",cex=0.8,col="blue", adj = c(1, 2))

node_to_length = sort(colMeans(res))
node_to_length = rbind(node_to_length,tree_true$edge.length[c(38,25,1,10,72,23,71,2,39,88)])





