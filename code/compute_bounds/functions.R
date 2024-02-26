library(rwty)
library(expm)
library(popdemo)
library(adephylo)
library(phytools)
library(ape)
library(castor)
library(latex2exp)
library(dplyr)
library(jsonlite)


# ----------- To compute bounds  ----------------

#Function which retrieve the data we need depending on the database
get_tree_par_fun = function(path,label,n_tree){
  
  # Extract the data from the json file 
  data = extract_data(path,label)
  
  # Extract tree(s) from the sample
  tree = path_to_trees(data$path,n_tree = 1)
  
  # Extract convergence parameters
  param = get_tree_param(tree,data$pi0,data$pi1,data$k,data$t_conv)
  
  # Functions to compute the graph
  
  # Topology function 
  f_topology = vectorize_function(compute_upper_bound_topology,param,type="topology")
  # Root function 
  f_root = vectorize_function(compute_upper_bound_root,param,type="root")
  
  return(
    list(tree=tree,param=param,f_topology=f_topology,f_root=f_root,path_cognates = data$path_cognates)
  )
}
  
# This function read a json file containing some parameters of the model :
# In the json file for each data base, you have acces to  
# - the path to the tree sample
# - pi1 : convergence parameters of the substitution model
# - pi0 : convergence parameters of the substitution model
# - The number of cognat
extract_data = function(path,label){
  dc_label_data = fromJSON(path)
  data = dc_label_data[[label]]
  
  return(list(
    path = data[[1]], 
    pi0 = as.numeric(data[[2]]), 
    pi1 = as.numeric(data[[3]]), 
    k=as.numeric(data[[4]]),
    path_cognates= data[[5]],
    t_conv = as.numeric(data[[6]])
    ))
}

# This function outputs the trees. 
# You can decide to output multiple trees with the parameter n_tree
# This function takes into account a burn-in of 20% 
path_to_trees = function(path,n_tree){
  phy = read.nexus(path)
  N = length(phy)
  
  # Burn-in 20%
  start = trunc(0.2*N)
  
  # When i want only one tree
  if (n_tree==1){return(phy[length(phy)])}
  
  # If I want more than one
  if (n_tree>1){
    n = sample(start:N,n_tree, replace=F)
    return(phy[n])
  }}

# This function compute the bounds
compute_bounds = function(path,pi0,pi1){
  
  # Import tree 
  tree = path_to_trees(path)
  # Parameters
  par = get_tree_param(tree,pi0,pi1)
  return(c(
    compute_upper_bound_root(par$t,par$Q,par$n),
    compute_upper_bound_topology(par$t,par$k,par$Q,par$n)
    ))
}

# This function compute the substitution model transition matrix
compute_chain = function(pi0,pi1){
  return(
    1/(pi0^2 + pi1^2) * as.matrix(rbind(c(-pi0,pi0),c(pi1,-pi1)))
  )}


get_tree_param = function(tree,pi0,pi1,k,t){
  
  Q = compute_chain(pi0,pi1)
  n = mean(sapply(tree, function(arbre) length(arbre$tip.label)))
  
  return(list(n=n,t=t,k=k,Q=Q))
}


compute_q = function(Q){
  d = as.numeric(dim(Q)[1])
  return(sum(apply(Q+100*diag(d),2,min)))
}

# This function compute the upper bound of the probability of the exact root 
# reconstruction
compute_upper_bound_root = function(t,Q,n){
  q = compute_q(Q)
  return(0.5 + n*exp(-q*t))
}

# This function compute the upper bound of the probability of the exact topology 
# reconstruction
compute_upper_bound_topology = function(t,k,Q,n){
  d = as.numeric(dim(Q)[1])
  q = sum(apply(Q+100*diag(d),2,min))
  return(k * n * exp(-q*t))
}

# ------------------- To compute graphics --------------------

vectorize_function = function(f,param,type){
  if (type=='root'){
    return (
      Vectorize(function(t) min(1,f(t,param$Q,param$n)))
    )}
  
  if (type == 'topology'){
    return(
      Vectorize(function(t) min(1,f(t,param$k,param$Q,param$n)))
    )}
  }


# -------------- Infimum of both bounds -----------------


find_t_value = function(k, Q, n, tolerance = 1e-6, max_iter = 1000) {
  
  objective_function = function(t) {
    return(compute_upper_bound_topology(t, k, Q, n) - 1)
  }
  
  result = uniroot(objective_function, interval = c(0, 20), tol = tolerance, maxiter = max_iter)
  
  return(result$root)
}

find_t_value_root = function(Q, n, tolerance = 1e-6, max_iter = 1000) {

  objective_function = function(t) {
    return(compute_upper_bound_root(t, Q, n) - 1)
  }

  result = uniroot(objective_function, interval = c(0, 20), tol = tolerance, maxiter = max_iter)
  return(result$root)
}




