library(here)
library(dplyr)
library(tibble)
library(ape)
library(TreeTools)
library(phytools)
library(adephylo)
library(castor)

remove_burnin = function(trees,burnin_rate){
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}

myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}


path= here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees")
trees = read.nexus(path)
trees = remove_burnin(trees,0.2)
M = length(trees)
tree = trees[[M]]






