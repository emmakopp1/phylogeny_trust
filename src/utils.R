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




path = "/Users/kopp/Documents/transeurasienne/beast/trans-ctmc-strict-fbd/data2.trees"
#tree = read.nexus(path)
#trees = rwty::load.trees(
#  path,
#  trim=1000,
#  log="/Users/kopp/Documents/transeurasienne/beast/trans-ctmc-strict-fbd/data2.log")


#consensus_tree = myconsensus(trees$trees)
#plot(consensus_tree,cex=0.5)




# Liste des langues à supprimer
path_out = "/Users/kopp/Documents/phylogeny_trust/data/real/sino-tibet-ctmc-strict-bd-fossilsRemoved/sino-tibetanfossilRemoved.nex"
path_in = "/Users/kopp/Documents/sino-tibetan/sino-tibetan.nex"

library(ape)
df = read.nexus.data(
  "/Users/kopp/Desktop/tea254.nex")

ancestors <- c("SiniticOldChinese", "TibetanOldTibetan", "Tangut", "BurmishOldBurmese")

# Supprimer les langues de la liste tt
df <- subset(df, !names(df) %in% ancestors)

write.nexus.data(df,file= path_out,format='standard',missing='?')

