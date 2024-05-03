library(here)
library(dplyr)
library(tibble)
source(here("src/init.R"))

remove_burnin = function(trees,burnin_rate){
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}
myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}


path= here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees")
trees = read.nexus(path)
M = length(trees)


fossils = c("BurmishOldBurmese","Tangut","SiniticOldChinese","TibetanOldTibetan")
tibble(
  age=min((rep(distRoot(trees[[M]],1),length(fossils)) - distRoot(trees[[M]],fossils))*1000)
  )
