rm(list=ls())
library(stringr)
library(ape)
library(phytools)
library(rwty)
library(castor)
library(TreeSim)
path = "/Users/kopp/Documents/chr_paper/reconstruction/sim/beast/tree-sim-"

simulate_tree = function(t,N=50,div = 0.245,r = 0.293,s = 0.226){
  # Parameters
  lambda = div/(1-r)
  mu = (r*div)/(1-r)
  psi = (s/(1-s))*(r*div/(1-r))
  
  # Tree simulation
  tree = sim.bd.taxa.age(
    n=N,
    numbsim=1, 
    lambda=lambda, 
    mu=mu, 
    frac = 1, 
    age=t, 
    mrca = TRUE
  )
  
  tree= tree[[1]]
  
  # Write tree 
  pathF = paste(path,t,".tree",sep="")
  print(pathF)
  write.nexus(tree,file=pathF)
}


