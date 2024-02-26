rm(list=ls())
library(ape,castor)
library(phytools)
library(phangorn,stringr)
library(tidyr)
library(TreeDist)
library(adephylo)
library(ggplot2)
library(stringr)
library(cowplot)
library(patchwork)

##### Import file
path_survival_functions = "/Users/kopp/Documents/chr_paper/code/cognat_survival_prob/function.R"
path_compute_bounds_function = "/Users/kopp/Documents/chr_paper/code/compute_bounds/functions.R"
path_to_config = "/Users/kopp/Documents/chr_paper/code/compute_bounds/config.json"

# Import functions
suppressWarnings({source(path_survival_functions)})
suppressWarnings({source(path_compute_bounds_function)})

# Data 
data_st = get_tree_par_fun(path_to_config,'sino-tibetan',1)
data_iecor = get_tree_par_fun(path_to_config,'iecor',1)
data_bantu = get_tree_par_fun(path_to_config,'bantu',1)


#---------------- 1. Survival probabilities -----------
# n_sens correspond to the number of wornd meaning in each datasets
res_st = compute_survival_prob_by_ages(data_st$tree[[1]],n_sens=200)
res_bantu = compute_survival_prob_by_ages(data_bantu$tree[[1]],n_sens=150)
res_iecor = compute_survival_prob_by_ages( data_iecor$tree[[1]],n_sens=170)

# Bind all res 
df = cbind(res_st$data,res_bantu$data,res_iecor$data)
df = df[, c(1,2,4,6)]
colnames(df) = c("age","q_st","q_bantu","q_iecor")

# Plots for each datasets
ggplot(df, aes(x = age)) +
  geom_line(aes(y = q_st,  linetype = "q_st"), linewidth = 0.5) +
  geom_line(aes(y = q_bantu, linetype = "q_bantu"), linewidth = 0.5) +
  geom_line(aes(y = q_iecor, linetype = "q_iecor"), linewidth = 0.5) +
  scale_color_manual(values = c("q_st" = "black", "q_bantu" = "black", "q_iecor" = "black")) +
  scale_linetype_manual(values = c("q_st" = "solid", "q_bantu" = "dashed", "q_iecor" = "dotted")) +
  labs(title = "", x = "Age", y = "Values") +
  theme_minimal()



#---------------- 2. Empirical value vs theoretical value ----------------
# Here I compute the expected number of meaning surviving from the root two subgroups of the root
# and the observed number of meaning which survived in the two subgroups 
# This part is still in progress

l_theo_empirical_st = t(as.data.frame(
  unname(sapply(data_st$tree, 
                function(arbre) 
                  c(
                    distRoot(arbre,1),
                    as.numeric(Q(arbre,find_root(arbre))*200),
                    survival_frequency(data_st$path_cognates,arbre)
                    ))),
  row.names = c("age","theorique","empirical")
  ))


l_theo_empirical_st

# ------------------- 3. Number of common meanings between a leaf and a groups --------

# In this part we consider that the subgroups of the root are the Sinitic vs Other

# This part is about checking the number of common meaning between a sinitic language
# and the other group 

# Separate the tip of tree into two set regarding to the subtrees induced by the root
tip_set = tips_of_subtree(data_st$tree[[1]])
tip_set$A # Others
tip_set$B # Chinese 


# Convert nexus data into dataframe
df_cognates = data_to_df(data_st$path_cognates)

# Dataframe of Sinitc group & other groups
v1 = df_cognates[tip_set$B,]
v2 = df_cognates[tip_set$A,]


resultats = c()

for (sinitic_row in tip_set$B){
  #x = apply(v2, 1, function(row) sum(v1[sinitic_row,] == row & v1[sinitic_row,] == 1, na.rm = TRUE))
  x = sum(apply(v2, 2, max,na.rm=1)&v1[sinitic_row,],na.rm=1)
  resultats = cbind(resultats,x)
}

colnames(resultats) = tip_set$B
resultats

# Calcul valeur théorique du nombre de cognat commun attendu

# Pour arbre sino-tibetain 
#t = distRoot(data_st$tree[[1]],1)[[1]]
#mu = 0.2052587
#val_theo = mean(rbinom(1000,200, 1/(t*mu)*(1-exp(-mu*t))^2))

#write.csv(resultats,
#          file= "/Users/kopp/Documents/chr_paper/code/cognat_survival_prob/cognat_commun2.csv",
#          sep=" ",
#          row.names=1,
#          col.names=1)





