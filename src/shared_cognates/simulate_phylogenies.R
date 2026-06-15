library(xml2); library(ape); library(phytools); library(glue)
library(here); library(castor); library(phangorn); library(tidyverse)
library(tracerer)

dir.create(here('data/shared_cognates'), showWarnings = FALSE)

# Load Sino Tibetan phylogenetic analysis outputs ------------------------------

tracelog_st <- parse_beast_tracelog_file(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.log"))
tracelog_st[round(nrow(tracelog_st) * 0.2) : nrow(tracelog_st), ]
pi1_density <- tracelog_st$freqParameter.s.sino.tibetan.2
pi0_density <- tracelog_st$freqParameter.s.sino.tibetan.1
clock_rate_density <- tracelog_st$clockRate.c.clock

q_density <- clock_rate_density/(2*pi0_density) + clock_rate_density/(2*pi1_density)

library(posterior)
trace_tail <- ess_tail(q_density)
trace_bulk <- ess_bulk(q_density)

# Define constants and substitution model 
#K_ie <- length(meanings_sets_ie$meaning)
#I_k_ie <- meanings_sets_ie$end - meanings_sets_ie$start + 1



# Paramètres -------------------------------------------------------------------
N_sim      <- 1
N_traits   <- 3000
pi1        <- 0.94305
pi0        <- 0.05695
clock_rate <- 0.018

Q <- matrix(
  c(-clock_rate/(2*pi0),  clock_rate/(2*pi0),
    clock_rate/(2*pi1), -clock_rate/(2*pi1)),
  nrow = 2, byrow = TRUE
)

q <- clock_rate/(2*pi0) + clock_rate/(2*pi1)

res        <- vector("list", 17)
root_state <- vector("list", 17)

# Boucle de simulation sur les âges d'arbre ------------------------------------
for (tree_age in 1:17) {
  
  tree_path <- glue(here(
    "data/simulated-2025-07-28/beast-data-sim-{N_sim}",
    "beast-data-sim-{N_sim}-{tree_age}",
    "tree-sim-{N_sim}-{tree_age}.tree"
  ))
  
  tree   <- read.tree(tree_path)
  n_tips <- length(tree$tip.label)
  mat    <- matrix(NA_real_, nrow = n_tips, ncol = N_traits)
  rownames(mat) <- tree$tip.label
  
  # Simulation des histoires de traits le long de l'arbre
  tt <- sim.history(tree, Q, nsim = N_traits, direction = 'row_to_column')
  
  # Identification de la racine et de ses deux sous-arbres descendants
  root          <- find_root(tree)
  root_children <- tree$edge[tree$edge[, 1] == root, 2]
  tips_A <- tree$tip.label[Descendants(tree, root_children[1], type = "tips")[[1]]]
  tips_B <- tree$tip.label[Descendants(tree, root_children[2], type = "tips")[[1]]]
  
  # Sauvegarde de l'état à la racine pour chaque trait (1 = présent, 0 = absent)
  root_state[[tree_age]] <- sapply(1:3000, function(i) {
    val <- tt[[i]]$node.states[root, 1]
    ifelse(val == "1", 1L, 0L)
  })
  
  # Toutes les paires de feuilles entre le groupe interne (A) et externe (B)
  pairs   <- expand.grid(A = tips_A, B = tips_B, stringsAsFactors = FALSE)
  N_pairs <- nrow(pairs)
  
  # Précalcul des chemins d'arêtes pour chaque paire de feuilles
  pair_edges <- vector("list", N_pairs)
  for (pair in seq_len(N_pairs)) {
    tipA_num   <- which(tree$tip.label == pairs[pair, 1])
    tipB_num   <- which(tree$tip.label == pairs[pair, 2])
    path_nodes <- nodepath(tree, tipA_num, tipB_num)
    
    edges <- integer(length(path_nodes) - 1)
    for (i in seq_len(length(path_nodes) - 1)) {
      nf  <- path_nodes[i]
      nt  <- path_nodes[i + 1]
      idx <- which(tree$edge[, 1] == nf & tree$edge[, 2] == nt)
      if (!length(idx)) idx <- which(tree$edge[, 1] == nt & tree$edge[, 2] == nf)
      edges[i] <- if (length(idx)) idx else NA_integer_
    }
    pair_edges[[pair]] <- edges[!is.na(edges)]
  }
  
  edge_root <- which(tree$edge[, 1] == root)[1]
  
  # Pour chaque trait : vérifier s'il était présent à la racine et
  # reconstructible (i.e. partagé par au moins une paire A/B sans homoplasie)
  for (trait in seq_len(N_traits)) {
    sim    <- tt[[trait]]
    states <- sim$states
    ns     <- sim$node.states
    
    # On ignore les traits absents à la racine
    if (ns[edge_root, 1] != "1") next
    
    tip_vals       <- states[tree$tip.label]
    trait_conserve <- FALSE
    
    for (pair in seq_len(N_pairs)) {
      tipA <- pairs[pair, 1]
      tipB <- pairs[pair, 2]
      
      # On ignore les paires où le trait est absent dans l'une des deux feuilles
      if (tip_vals[tipA] != "1" || tip_vals[tipB] != "1") next
      
      edges <- pair_edges[[pair]]
      if (!length(edges)) {
        trait_conserve <- TRUE; break
      }
      
      # Détection vectorisée des naissances homoplasiques (0 -> 1) sur le chemin
      naissances <- any(ns[edges, 1] == "2" & ns[edges, 2] == "1")
      if (!naissances) {
        trait_conserve <- TRUE; break
      }
    }
    
    if (trait_conserve)
      mat[, trait] <- ifelse(tip_vals[rownames(mat)] == "1", 1L, 0L)
  }
  
  res[[tree_age]] <- mat
}

# Nombre de cognats partagés entre deux sous arbres-----------------------------
shared_cognates <- data.frame(
  tree_age = 1:17,
  value    = numeric(17),
  inf      = numeric(17),
  sup      = numeric(17)
)

for (tree_age in 1:17) {
  tree_path <- glue(here(
    "data/simulated-2025-07-28/beast-data-sim-{N_sim}",
    "beast-data-sim-{N_sim}-{tree_age}",
    "tree-sim-{N_sim}-{tree_age}.tree"
  ))
  tree <- read.tree(tree_path)
  root <- find_root(tree)
  
  root_children <- tree$edge[tree$edge[, 1] == Ntip(tree) + 1, 2]
  tips_A <- tree$tip.label[Descendants(tree, root_children[1], type = "tips")[[1]]]
  tips_B <- tree$tip.label[Descendants(tree, root_children[2], type = "tips")[[1]]]
  
  mat_A <- res[[tree_age]][tips_A, ]
  mat_B <- res[[tree_age]][tips_B, ]
  
  # Un cognat est "partagé" s'il est présent dans au moins une feuille de chaque sous-arbre
  present_A    <- colSums(mat_A == 1, na.rm = TRUE) > 0
  present_B    <- colSums(mat_B == 1, na.rm = TRUE) > 0
  shared_sites <- present_A & present_B
  
  # Test binomial contre H0 : p = 0.22
  test_binomial <- prop.test(
    x          = sum(shared_sites),
    n          = sum(root_state[[tree_age]]),
    p          = 0.22,
    conf.level = 0.95
  )
  
  shared_cognates$value[tree_age] <- sum(shared_sites) / sum(root_state[[tree_age]])
  shared_cognates$inf[tree_age]   <- test_binomial$conf.int[1]
  shared_cognates$sup[tree_age]   <- test_binomial$conf.int[2]
}

saveRDS(shared_cognates, here("output/results/shared_cognates.rds"))
#shared_cognates <- readRDS(here("output/results/shared_cognates.rds"))

# Proportion de cognats partagés entre une paire de feuilles fixe (avec homoplasie) ----
# Tirage d'une feuille dans chaque sous-arbre (fixe pour tous les âges)
tree_ref <- read.tree(glue(here(
  "data/simulated-2025-07-28/beast-data-sim-{N_sim}",
  "beast-data-sim-{N_sim}-1",
  "tree-sim-{N_sim}-1.tree"
)))
root_ref      <- find_root(tree_ref)
root_children <- tree_ref$edge[tree_ref$edge[, 1] == root_ref, 2]
tips_A <- tree_ref$tip.label[Descendants(tree_ref, root_children[1], type = "tips")[[1]]]
tips_B <- tree_ref$tip.label[Descendants(tree_ref, root_children[2], type = "tips")[[1]]]

tip_A <- sample(tips_A, 1)
tip_B <- sample(tips_B, 1)

# Pour chaque âge d'arbre, calcul de la proportion de cognats partagés par la paire
# (normalisée par le nombre de cognats présents à leur ancêtre commun)
prop_shared <- numeric(17)
for (tree_age in 1:17) {
  mat <- res[[tree_age]]
  # un cognat présent à la racine à au bout de 5000 ans, 60% de chance d'etre présent dans une des langues des deux coté
  
  prop_shared[tree_age] <- sum(mat[tip_A, ] == 1 & mat[tip_B, ] == 1, na.rm = TRUE) /
    sum(root_state[[tree_age]])
}

# save dataframe 
prop_shared_tip_pair <- tibble(
  tipA = tip_A,
  tipB = tip_B,
  prop = prop_shared
)

saveRDS(prop_shared_tip_pair, here("output/results/shared_cognate_tip_pair.rds"))
#prop_shared <- readRDS(here("output/results/shared_cognate_tip_pair.pdf"))


# Proportion de cognats partagés entre une paire de feuilles fixe (sans homoplasie) 







