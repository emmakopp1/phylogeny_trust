library(ape)
library(castor)
library(phytools)
library(phangorn)
library(stringr)
library(tidyr)



# ----------- Probabilty of survival in both subtrees induces by the root -----------

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
  
  
  #p = ggplot(df_theo, aes(x = age, y = q_theo)) +
  #  geom_line(color = "darkblue") +
  #  labs(
  #    x = "age of the tree (in millenial)",
  #    y = "",
  #    title = ""
  #  )+
  #  theme(plot.title = element_text(hjust = 0.5))
  
  return(df_theo)
  
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









