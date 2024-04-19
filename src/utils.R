library(here)
source(here("src/init.R"))


remove_burnin = function(trees,burnin_rate){
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}
myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}


library(rwty)
library(treeio)
path = here("data/real/st_ctmc-strict-bd-fossilsRemoved/st_ctmc-strict-bd-fossilRemoved.trees")
trees = rwty::load.trees(path,trim=100)
trees = trees$trees
trees = remove_burnin(trees,0.1)


tips = c(
  "TibetanAlike",
  "TibetanBatang",
  "TibetanLhasa",
  "TibetanXiahe")

mrca_age = function(tree,tips){
  
  mrca = treeio::MRCA(tree,tips)
  age_mrca = distRoot(tree,11)[[1]] - distRoot(tree,mrca)[[1]]
  return(age_mrca)
  
}

suppressWarnings(
  mean(sapply(trees, function(tree) mrca_age(tree,tips)))
  )




