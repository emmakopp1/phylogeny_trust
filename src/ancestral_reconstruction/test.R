library(here)
library(ape)
library(phytools)
library(castor)
library(phangorn)
library(adephylo)
library(dplyr)

# Charger les données IE
phylo_ie <- read.nexus(here("data/real/iecor_ctmc-strict-M1/IECoR_M1_CTMC_Gamma_1_Rate_For_All_Mgs_combined.trees"))
M_ie <- length(phylo_ie)
phylo_ie <- phylo_ie[seq(0.8 * M_ie, M_ie, length = 200)]
M_ie <- length(phylo_ie)

Y_ie <- read.nexus.data(here("data/real/iecor_ctmc-strict-M1/iecor.nex"))
Y_ie <- lapply(Y_ie, function(col) replace(col, col == "?", NA))
Y_ie <- lapply(Y_ie, as.numeric)
Y_ie <- as.data.frame(Y_ie)

meanings_sets_ie <- read.csv(here("output/results/meanings_sets_ie.csv"))
bounds_real_tb_ie <- read.csv(here("output/results/tracelog_summary.csv"))
bounds_real_tb_ie <- bounds_real_tb_ie[bounds_real_tb_ie$family == "IE", ]

pi_ie <- bounds_real_tb_ie[, c("pi0", "pi1")]
lambda_ie <- 1 / (2 * pi_ie$pi0)
mu_ie <- 1 / (2 * pi_ie$pi1)
Q_ie <- cbind(c(-lambda_ie, mu_ie), c(lambda_ie, -mu_ie))

# Transposer Y pour avoir langues en lignes
Y_full <- t(Y_ie)

# Tester uniquement sur les traits 4699 et 4700
traits_test <- c(4699, 4700)
results <- list()

for (trait_global in traits_test) {
  y_full <- Y_full[, trait_global]
  
  for (t in 1:M_ie) {
    tree_pruned <- phylo_ie[[t]]
    y <- y_full
    ii <- which(is.na(y))
    
    if (length(ii) > 0) {
      y <- y[-ii]
      tree_pruned <- drop.tip(tree_pruned, names(ii))
    }
    
    if (sum(y, na.rm = TRUE) > 1 & length(unique(y)) > 1) {
      rec <- ancr(fitMk(tree_pruned, y, "ARD", fittedQ = Q_ie, pi = as.numeric(pi_ie)))
      indice <- which(rec[["ace"]][, 2] > 0.5)
      nodes <- as.integer(names(indice))
      
      if (length(nodes) > 0) {
        depths <- as.numeric(max(distRoot(tree_pruned)) - distRoot(tree_pruned, nodes))
        best_idx <- which.max(depths)
        value <- depths[best_idx]
        best_node <- nodes[best_idx]
        
        results[[length(results) + 1]] <- data.frame(
          value = value,
          tree = t,
          trait = trait_global,
          node = best_node
        )
      }
    }
  }
}

res_df <- do.call(rbind, results)
saveRDS(res_df, here("output/results/ancr_t4699_t4700.rds"))
res_df = readRDS(here("output/results/ancr_t4699_t4700.rds"))

head(res_df)
t4699 = res_df |> filter(trait == 4699, tree==95)

mean(t4699$value, na.rm=T)

t4700 = res_df |> filter(trait == 4700)

mean(t4700$value, na.rm=T)


tree4699 = readRDS(here("output/trees/ancr/tree_pruned_ie/tree_water_t95_trait4699.rds"))
tree4700 = readRDS(here("output/trees/ancr/tree_pruned_ie/tree_water_t132_trait4700.rds"))

# trait 4699
plot(tree4699, cex=0.4)
nodelabels(node = 294, cex=0.4)

node_4699 = Descendants(tree4699, node=294)[[1]]
# descendant du noeud reconstruit 
tree4699$tip.label[node_4699]
# age du noeud reconstruit 
max(distRoot(tree4699))- distRoot(tree4699, 294) 

# trait 4700
plot(tree4700,cex=0.4)
nodelabels(node = 268, cex=0.4)

node_4700 = Descendants(tree4700, node=268)[[1]]
# descendant du noeud reconstruit 
tree4700$tip.label[node_4700]
# age du noeud reconstruit 
max(distRoot(tree4700))- distRoot(tree4700, 268) 


# extraction des traits présents dans le fichier nexus
Y_test <- Y_full[rownames(Y_full) != "OldBreton", traits_test]

# trait 4699 - cognate set 335 (european)
rownames(Y_test)[Y_test[, 1] == TRUE]
# trait 4700 - cognate set 157 (iranien)
rownames(Y_test)[Y_test[, 2] == TRUE]

# la racine remonte a des noeuds et donc la racine aurait été innové dans 10 branche 
# voir les endroit où le trait a été innové et combien de fois il a été innové. 
# et est ce que c'est réaliste ? 
# la ca serait une homoplasie massive 
# il aurait été innové indépendemment au moins 6 fois sur l'autre groupe 
# on s'attendrait a ce que ce trait puisse etre reconstruit à la racine 

# on peut vérifier les test en voyant les homoplasie
# si deux mots on un trait présent mais structure morphologique différentes 

# il faut un moyen de pénaliser les homoplasie 

# glissement sémentique cas particulier
df_semantiquet_t4699_t4700 <- read.csv(here("output/results/ancestral_reconstruction_ie_glissement_sementiques_t4699_t4700.csv"), row.names = NULL) 
colnames(df_semantiquet_t4699_t4700) = c('node', 'p_node','sens', 'node_parent', 'p_parent', 'tree', 'trait')

# traits 4700
df_semantiquet_t4700 = read.csv(here("output/results/ancestral_reconstruction_ie_glissement_sementiques_t4700.csv"), row.names = NULL) 
colnames(df_semantiquet_t4700) = c('node', 'p_node','sens', 'node_parent', 'p_parent', 'tree', 'trait')

head(df_semantiquet_t4700)

tt = df_semantiquet_t4700 |> 
  select(node, p_node) |> 
  filter(p_node > 0.5)|> 
  arrange(desc(p_node)) |> 
  pull(node) |> 
  as.numeric() |> 
  na.omit() |> 
  as.vector()

plot(tree4700, cex = 0.4)

tiplabels(
  tip = which(Y_full[tree4700$tip.label, 4700] == 1),
  pch = 16,
  col = "red",
  cex = 0.6
)
nodelabels(node=tt, cex=0.4)

# traits 4699
df_semantiquet_t4699= read.csv(here("output/results/ancestral_reconstruction_ie_glissement_sementiques_t4699.csv"), row.names = NULL) 
colnames(df_semantiquet_t4699) = c('node', 'p_node','sens', 'node_parent', 'p_parent', 'tree', 'trait')

head(df_semantiquet_t4699)

tt = df_semantiquet_t4699 |> 
  select(node, p_node) |> 
  filter(p_node > 0.5)|> 
  arrange(desc(p_node)) |> 
  pull(node) |> 
  as.numeric() |> 
  na.omit() |> 
  as.vector()

plot(tree4699, cex = 0.4)

tiplabels(
  tip = which(Y_full[tree4699$tip.label, 4699] == 1),
  pch = 16,
  col = "red",
  cex = 0.6
)
nodelabels(node=tt, cex=0.4)









