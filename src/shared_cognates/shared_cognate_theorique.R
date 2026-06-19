# Libraries 
library(here) 
library(ape)
library(phytools)
library(castor)
library(phangorn)
library(adephylo)
length_phylo <- 200
# Load BEAST posterior trees and thin sample 
phylo_st <- read.nexus(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees"))
M_st <- length(phylo_st)
phylo_st <- phylo_st[seq(0.8 * M_st, M_st, length = length_phylo)]
M_st <- length(phylo_st)  # update M to number of retained trees


# Load meaning boundaries and filter parameters
meanings_sets_st <- read.csv(here("output/results/meanings_sets_st.csv"))
bounds_real_tb_by_sens <- read.csv(here("output/results/tracelog_summary.csv"))
bounds_real_tb_st <- bounds_real_tb_by_sens[bounds_real_tb_by_sens$family == "ST", ]

# Define constants and substitution model 
#K_st <- length(meanings_sets_st$meaning)
#I_k_st <- meanings_sets_st$end - meanings_sets_st$start + 1

pi_st <- bounds_real_tb_st[, c("pi0", "pi1")]
clock_rate_st = 0.018
lambda_st <- clock_rate_st / (2 * pi_st$pi0)
mu_st <- clock_rate_st / (2 * pi_st$pi1)
Q_st <- cbind(c(-lambda_st, mu_st), c(lambda_st, -mu_st))

#exp(-2*mu_st)

# take the last tree
tree = phylo_st[[length_phylo]]
root = find_root(tree)
root_age = max(distRoot(tree))

# S(x) : probabilité qu'il existe un descendant d de x tels que le trait -------
# est présent apartout sur <x,d> | présent en x 
# on veut calculer S(root)

S <- function(node){
  
  if (node %in% 1:Ntip(tree)){return(1)}
  else{
    # children of node
    children = Descendants(tree, node, type = "children")
    
    # index of edges from node to children in the tree edges
    branch_length_chidren_index = sapply(children,function(i) which(tree$edge[,1] == node & tree$edge[,2] == i))
    
    # (exponential )branch length from <node,y> and <node,z> 
    branch_length_children = tree$edge.length[branch_length_chidren_index]
    exp_branch_length_children = exp(-mu_st * branch_length_children)
    
    return(
      1 - (1 - exp_branch_length_children[1] + exp_branch_length_children[1]*(1 - S(children[1])))
      * (1 - exp_branch_length_children[2] + exp_branch_length_children[2]*(1 - S(children[2])))
    )
  }
  
}

first_split = Descendants(tree, root, type = "children")

# index of edges from node to children in the tree edges
branch_length_chidren_index = sapply(first_split,function(i) which(tree$edge[,1] == root & tree$edge[,2] == i))

# (exponential )branch length from <node,y> and <node,z> 
branch_length_children = tree$edge.length[branch_length_chidren_index]
exp_branch_length_children = exp(-mu_st * branch_length_children)

exp_branch_length_children[1]*S(first_split[1]) * exp_branch_length_children[2]*S(first_split[2])
# je veux proba qu'il est survecu des coté 
# sur 200 mots a la racine 

# Calcul de T_1(x) : il existe des cognat non homoplasique
T_1 <- function(node){
  if (node %in% 1:Ntip(tree)){return(1)}
  else{
    children = Descendants(tree, node, type = "children")
    
    # index of edges from node to children in the tree edges
    branch_length_chidren_index = sapply(children,function(i) which(tree$edge[,1] == node & tree$edge[,2] == i))
    
    # (exponential )branch length from <node,y> and <node,z> 
    branch_length_children = tree$edge.length[branch_length_chidren_index]
    exp_branch_length_children = exp(-mu_st * branch_length_children)
    
    return(
      1 - (exp(-(1-lambda_st) * branch_length_children[1]) * (1 - T_1(children[1])) +
             exp(-lambda_st * branch_length_children[1]) * (1 - T_0(children[1])) )
      * (exp(-(1-lambda_st) * branch_length_children[2]) * (1 - T_1(children[2])) +
           exp(-lambda_st * branch_length_children[2]) * (1 - T_0(children[2])) )
    )
  }
}

T_0 <- function(node){
  if (node %in% 1:Ntip(tree)){return(1)}
  else{
    children = Descendants(tree, node, type = "children")
    
    # index of edges from node to children in the tree edges
    branch_length_chidren_index = sapply(children,function(i) which(tree$edge[,1] == node & tree$edge[,2] == i))
    
    # (exponential )branch length from <node,y> and <node,z> 
    branch_length_children = tree$edge.length[branch_length_chidren_index]
    exp_branch_length_children = exp(-mu_st * branch_length_children)
    
    return(
      1 - (exp(-(1-mu_st) * branch_length_children[1]) * (1 - T_0(children[1])) +
             exp(-mu_st * branch_length_children[1]) * (1 - T_1(children[1])) )
      * (exp(-(1-mu_st) * branch_length_children[2]) * (1 - T_0(children[2])) +
           exp(-mu_st * branch_length_children[2]) * (1 - T_1(children[2])) )
    )
  }
}
Descendants(tree,51, type= "children")

T_1(52) * (1 - exp(-mu_st * branch_length_children[1])) + 
  T_1(58) * (1 - exp(-mu_st * branch_length_children[2])) - 
  T_1(52) * (1 - exp(-mu_st * branch_length_children[1])) *
  T_1(58) * (1 - exp(-mu_st * branch_length_children[2]))

# Sensitivity of T(root) to tree age -------------------------------------------
first_split = Descendants(tree, root, type = "children")

target_ages <- seq(1, 17, by = 0.1)
coefs       <- target_ages / root_age   # multiplicative scaling factors

results_T <- data.frame(
  tree_age = numeric(),
  coef     = numeric(),
  T_root   = numeric()
)

for (i in seq_along(target_ages)) {
  
  # Scale branch lengths
  tree_scaled            <- tree
  tree_scaled$edge.length <- tree$edge.length * coefs[i]
  
  # index of edges from node to children in the tree edges
  branch_length_chidren_index = sapply(first_split,function(i) which(tree_scaled$edge[,1] == root & tree_scaled$edge[,2] == i))
    
  # (exponential )branch length from <node,y> and <node,z> 
  branch_length_children = tree_scaled$edge.length[branch_length_chidren_index]
    
  # peut etre ca 
  #s_val <- T_1(first_split[1]) * ( 1 - exp(-mu_st * branch_length_children[1])) + 
  #  T_1(first_split[2]) * (1 - exp(-mu_st * branch_length_children[2])) - 
  #  T_1(first_split[1]) * ( 1 - exp(-mu_st * branch_length_children[1])) *
  #  T_1(first_split[2]) * (1 - exp(-mu_st * branch_length_children[2]))
  
  # peut etre ca
  s_val <-  1 - ((1 - T_1(first_split[1]) * ( 1 - exp(-mu_st * branch_length_children[1])))*
    (1- T_1(first_split[2]) * (1 - exp(-mu_st * branch_length_children[2]))))
  
  results_T <- rbind(results_T, data.frame(
    tree_age = target_ages[i],
    coef     = coefs[i],
    T_root   = s_val + pi1
  ))
}

shared_cognate_comp <- readRDS(here("output/results/shared_cognates.rds"))

# ── Plot ───────────────────────────────────────────────────────────────────────

plot(
  results_T$tree_age, results_T$T_root,
  type = "l", lwd = 2, col = "darkblue",
  xlab = "Age (ka BP)",
  ylab = "",
  main = ""
)
abline(v = root_age, lty = 2, col = "tomato", lwd = 1.5)
legend("topright",
       legend = c("T(root)", sprintf("Original age (%.2f)", root_age)),
       col    = c("darkblue", "tomato"),
       lty    = c(1, 2), lwd = 2)



# Définir les limites communes
xlim <- range(c(results_T$tree_age, shared_cognate_comp$tree_age))
ylim <- range(c(results_T$T_root, shared_cognate_comp$value))

# Premier plot (tous les 1)
plot(
  results_T$tree_age, log(results_T$T_root),
  type = "l", lwd = 2, col = "darkblue",
  xlab = "Tree age (root depth)",
  ylab = "",
  main = "",
  xlim = xlim
)

curve(exp(-2*x*mu_st), add=T) # avec 2 langues au bout de 10k ans, 5% des cognats sont gardés
# mais (ligne bleu + de diversit -> + de chance de survivre -> pente plus douce)


# T_1(Root) devrait matcher cette fracion
plot(results_T$T_root[seq(1,161,by=10)]/shared_cognate_comp$value) 
# proba qu'un trait a la racine est présent dans au moins une feuille (evenement conditionelle qui nous ennuit)
lines(results_T$T_root[seq(1,161,by=10)]/shared_cognate_comp$value)


