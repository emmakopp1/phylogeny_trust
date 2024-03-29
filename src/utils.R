library(here)
source(here("src/init.R"))

remove_burnin = function(trees,burnin_rate){
  #'remove_burnin
  #'
  #'Remove the burn-in of a sample of trees
  #'@param trees
  #'@param burnin_rate the rate of burn-in we apply
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}
myconsensus = function(trees){
  #' myconcensus
  #' 
  #' Compute the consensus tree of multiples trees
  #' @param trees The trees.
  #' @example 
  #' path = "/Users/kopp/Documents/chr_paper/beast/bantu-ctmc-strict-bd/ctmc-strict-bd.trees"
  #' tree = read.nexus(path)
  #' tree = remove_burnin(tree,0.9)
  #' consensus_tree = myconsensus(tree)
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}


# For each trees
'for (i in 8:14){
  path = sprintf("/Users/kopp/Documents/chr_paper/reconstruction/one-tree/beast-data-sim-%i/ctmc-strict-bd-%i.trees", i,i)
  tree = read.nexus(path)
  tree = remove_burnin(tree,0.1)
  consensus_tree = myconsensus(tree)
  
  write.tree(consensus_tree, 
             file = sprintf("/Users/kopp/Documents/chr_paper/reconstruction/one-tree/beast-data-sim-%i/consensus-%i.tree",i,i))
}


tree = read.nexus()
'
