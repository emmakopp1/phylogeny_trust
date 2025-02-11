#sink("src/ancestral_reconstruction/test.log")
#print('test')
#sink()
#t1<-Sys.time()

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phytools", quietly = TRUE)) install.packages("phytools")
if (!requireNamespace("castor", quietly = TRUE)) install.packages("castor")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")
if (!requireNamespace("rwty", quietly = TRUE)) install.packages("rwty")
if (!requireNamespace("adephylo", quietly = TRUE)) install.packages("adephylo")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

library(parallel)
library(here) 
library(ape)
library(phytools)
library(castor)
library(phangorn)
library(rwty)
library(adephylo)
library(reshape2)


# output parameters, tree and burnin
bounds_real_tb_by_sens <- read.csv("/home/users/kopp/work/ancestral_reconstruction/bounds_real_tb_by_sens.csv")
tree <- read.nexus("/home/users/kopp/work/ancestral_reconstruction/st_ctmc-strict-fbd-uniform.trees")
M <- length(tree)
tree<- tree[seq(0.2*M,M, length = 1000)] # thin-in
#tree = tree[M]
M <- length(tree)

# check treeHeight
ntips = Ntip(tree[1])
mean(sapply(tree, function(single_tree) max(node.depth.edgelength(single_tree)[1:ntips])))

# data
Y <- read.nexus.data('/home/users/kopp/work/ancestral_reconstruction/st.nex')
Y <- lapply(Y, function(col) replace(col, col == "?", NA))
Y <- lapply(Y, as.numeric)
Y <- as.data.frame(Y)


# indexes for each indice in the data
meanings_sets <- read.csv('/home/users/kopp/work/ancestral_reconstruction/meanings_sets.csv')

# sino tibetan study for homogeneous traits
bounds_real_tb <- bounds_real_tb_by_sens[bounds_real_tb_by_sens$family == "ST", ]

# number of meanings
K <- length(meanings_sets$meaning)
I_k <- meanings_sets$end - meanings_sets$start + 1

# substituion parameters
pi <- bounds_real_tb[, c("pi0", "pi1")]
lambda <-1/(2*pi$pi0)
mu <-1/(2*pi$pi1)
Q <-  cbind(c(-lambda,mu),c(lambda,-mu))

# Nombre de cœurs disponibles

ncl = 30
cl = makeCluster(ncl, type="FORK")
clusterSetRNGStream(cl)

process_k <- function(k) {
  
  meaning_k <- meanings_sets$meaning[k]
  I <- I_k[k]
  
  # data indices
  start <- as.integer( meanings_sets[meanings_sets$meaning == meaning_k, c("start") ])
  end <- as.integer( meanings_sets[meanings_sets$meaning == meaning_k, c("end") ])

  # data
  Y_pruned <- t(Y[start:end, ])
  
  # results
  results <- list()
  
  for (t in 1:M) {
    tree_pruned <- tree[[t]]
    for (trait in 1:ncol(Y_pruned)) {
      y <- Y_pruned[,trait] # vecteur avec des noms
      ii <- which(is.na(y))
      if (length(ii)>0){
        y <- y[-ii]
        tip_to_drop <- names(ii)
        tree_pruned <- drop.tip(tree_pruned, tip_to_drop)
        }
      
      if (sum(y, na.rm=T) > 1) {
        
        rec <- ancr(fitMk(tree_pruned, y, "ARD", fittedQ = Q, pi = as.numeric(pi)))

        indice <- which(rec[["ace"]][,2] > 0.5)
        nodes <- as.integer(names(indice))
        if(length(nodes) >0 ){
          value <- max(as.numeric(distRoot(tree_pruned,1) - distRoot(tree_pruned,nodes)))
          results[[length(results) + 1]] <- c(value = value, sens = meaning_k, tree = t)
        }
      }
    }
  }
  
  # Vérification avant la conversion
  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    # Retourner une matrice vide avec des noms de colonnes corrects
    matrix(NA, nrow = 0, ncol = 3, dimnames = list(NULL, c("value", "sens", "tree")))
  }
}



clusterExport(cl, varlist = c("process_k"))
res_list <- parLapply(cl, 1:K, function(k) {
  sprintf("Starting process_k for k = %d", k) # Ajouter un message de débogage
  tryCatch(
    {
      result <- process_k(k)
      message(sprintf("Result for k = %d: %s", k, result))  # Message de débogage après l'exécution de process_k
      return(result)
    },
    error = function(e) {
      print(sprintf("Error in process_k for k = %d: %s", k, e$message))
      return(matrix(NA, nrow = 0, ncol = 3, dimnames = list(NULL, c("value", "sens", "tree"))))
    }
  )
})

# save results
res_combined <- do.call(rbind, res_list) 
res_combined <- as.data.frame(res_combined)


write.table(res_combined, file ='/home/users/kopp/work/ancestral_reconstruction/results.txt', sep = "\t", row.names = FALSE, quote = FALSE)
stopCluster(cl)



# Combinaison des résultats
#res_df <- do.call(rbind, lapply(res_list, as.data.frame))

#res_df <- bind_rows(lapply(res_list, as.data.frame)) |>
#  mutate(across(everything(), ~ ifelse(is.na(.), NA, .))) |>
#  drop_na()








