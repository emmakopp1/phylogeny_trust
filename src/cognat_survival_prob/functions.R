library(ape)
library(castor)
library(phytools)
library(phangorn)
library(stringr)
library(tidyr)



# ----------- Probabilty of survival in both subtrees induces by the root -----------

# This function output the branch length in a tree between two nodes
find_branch_length <- function(tree, father, children) {
  indice <- which(tree$edge[, 1] == father & tree$edge[, 2] == children)
  return(tree$edge.length[indice])
}

# This function outputs the probability of a cognat to survive between node i to node j 
delta <- function(tree, i, j, mu= 0.2052587) {
  return(exp(-mu * find_branch_length(tree, i, j)))
}


# See notations in article Gray&Nicholls 2018
# Sequence u
u <- function(tree, i, mu= 0.2052587) {
  if (i %in% 1:(tree$Nnode + 1)) {
    return(0)
  } 
  else {
    childrens <- Descendants(tree, i, type = c("children"))
    return(
      (1 - delta(tree, i, childrens[1], mu) + delta(tree, i, childrens[1],mu) * u(tree, childrens[1], mu)) *
        (1 - delta(tree, i, childrens[2], mu) + delta(tree, i, childrens[2], mu) * u(tree, childrens[2], mu))
    )}}

# Probability of a trait born in i to survive in 1 or more leave
P <- function(tree, i, mu= 0.2052587) {
  return(1 - u(tree, i, mu))
}


# Prob of a trait born in i to survive in bot subtree
Q <- function(tree, i, mu= 0.2052587) {
  childrens <- Descendants(tree, i, type = c("children"))
  return(
    delta(tree, i, childrens[1], mu) * P(tree, childrens[1], mu) *
      delta(tree, i, childrens[2], mu) * P(tree, childrens[2], mu)
  )
}


# Separate the tip of tree into two set regarding to the subtrees induced by the root
tips_of_subtree = function(tree){
  root = find_root(tree)
  child1 = Descendants(tree, root, type = c("children"))[[1]]
  
  tipsA = tree$tip.label[Descendants(tree, child1, type = c("tips"))[[1]]]
  tipsB = setdiff(tree$tip.label,tipsA)
  
  return(list(A=tipsA, B=tipsB))
} 


# From .txt file to dataframe
data_to_df = function(path){
  # Checking
  
  data = readLines(path)
  
  # Initialisation
  x = str_squish(data[1])
  x_name = unlist(strsplit(x," "))[1]
  x_vec = unlist(strsplit(x," "))[2]
  x_vec = unlist(strsplit(x_vec, ""))
  
  df = data.frame(
    Name = x_name, 
    setNames(as.data.frame(matrix(x_vec, ncol = length(x_vec),byrow = TRUE)), paste0("Bit", 1:length(x_vec)))
  )
  
  # Recurence
  for (k in 2:length(data)){
    x = str_squish(data[k])
    
    # Tip name
    x_name = unlist(strsplit(x," "))[1]
    
    # Vecteur
    x_vec = unlist(strsplit(x," "))[2]
    x_vec = unlist(strsplit(x_vec, ""))
    
    df=rbind(df, c(x_name,x_vec))
  }
  
  # Index
  names = df$Name
  df = df[,2:length(df)]
  df = as.data.frame(sapply(df, as.integer))
  rownames(df) = names
  
  return(df)
}

# ------------------- Survival probability  --------------------

# Theoretical value
compute_survival_prob_by_ages = function(tree,n_sens,mu= 0.2052587){
  
  t = as.numeric(distRoot(tree,1))
  root = find_root(tree)
  ks = (1:20)/t
  
  # create 20 trees
  for (k in ks) {
    tree_k = tree
    tree_k$edge.length = tree_k$edge.length * k 
    assign(paste0("tree_", k*t), tree_k)
  }
  
  q_theo =c()
  for (k in 1:20){
    tree_k = get(paste0("tree_",k))
    q_theo=c(q_theo,Q(tree,root,mu*k/t)*n_sens)
  }
  
  
  df_theo = data.frame(
    age = 1:length(q_theo),
    q_theo = q_theo
  )
  
  
  p = ggplot(df_theo, aes(x = age, y = q_theo)) +
    geom_line(color = "darkblue") +
    labs(
      x = "age of the tree (in millenial)",
      y = "",
      title = ""
    )+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list(
    data=df_theo, plot=p
  ))
  
}

# Empirical value
survival_frequency = function(path_cognates,tree){
  
  root = find_root(tree)
  
  # Subtrees
  df = data_to_df(path_cognates)
  rownames(df) = gsub("'",'', rownames(df))
  tip_set = tips_of_subtree(tree)
  dfA = df[tip_set$A ,]
  dfB = df[tip_set$B,]
  
  return(
    # Nombre de personne dans groupA qui ont chaque traits 
    as.numeric(sum(colSums(dfA, na.rm=T) >=1 & colSums(dfB, na.rm=T) >=1))
  )
  
}









