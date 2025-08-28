library(here) 

path_seq1500 = here("data/simulated-2025-07-21/beast-data-sim-1/beast-data-sim-1-1/beast-simulated-seq-1-1.xml")

# Définir les paths
paths <- list(
  "1500" = here("data/simulated-2025-07-22-1500/beast-data-sim-1/beast-data-sim-1-8/beast-simulated-seq-1-8.xml"),
  #"3000" = here("data/simulated-2025-07-22/beast-data-sim-1/beast-data-sim-1-8/beast-simulated-seq-1-8.xml"),
  "6000" = here("data/simulated-2025-07-22-6000/beast-data-sim-1/beast-data-sim-1-8/beast-simulated-seq-1-8.xml"),
  "12000" = here("data/simulated-2025-07-22-12000/beast-data-sim-1/beast-data-sim-1-8/beast-simulated-seq-1-8.xml")
)

# Fonction pour analyser un fichier
analyze_sequences <- function(path) {
  # Lire le fichier
  content <- readLines(path, warn = FALSE)
  
  # Extraire les séquences
  sequences <- regmatches(content, gregexpr("value='[01]+'", content))
  sequences <- unlist(sequences)
  sequences <- gsub("value='|'", "", sequences)
  
  # Compter les "1" dans chaque séquence
  ones_counts <- sapply(sequences, function(x) sum(strsplit(x, "")[[1]] == "1"))
  
  # Retourner les statistiques
  list(
    n_sequences = length(sequences),
    seq_length = if(length(sequences) > 0) nchar(sequences[1]) else NA,
    mean_ones = round(mean(ones_counts), 2),
    percent_ones = round(mean(ones_counts) / nchar(sequences[1]) * 100, 2)
  )
}

# Analyser tous les fichiers
results <- lapply(paths, analyze_sequences)

# Créer un data frame pour l'affichage
df_results <- data.frame(
  Dataset = names(paths),
  N_Seq = sapply(results, function(x) x$n_sequences),
  Seq_Length = sapply(results, function(x) x$seq_length),
  Mean_Ones = sapply(results, function(x) x$mean_ones),
  Percent_Ones = sapply(results, function(x) x$percent_ones)
)

# Affichage concis
cat("=== Analyse des séquences BEAST ===\n\n")
print(df_results, row.names = FALSE)


tt=read.tree('/Users/kopp/Documents/phylogeny_trust/data/simulated-2025-07-22/beast-data-sim-8/beast-data-sim-8-17/tree-sim-8-17.tree')
plot(tt)
