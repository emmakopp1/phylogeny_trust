# ------------------------------------------------------------------------------
# Script Name: 00_compute_meaning_set.R
# Description: This script processes phylogenetic NEXUS files to extract and 
#              organize meaning groups from character state labels for ancestral 
#              reconstruction analysis. It performs the following steps:
#                - Reads NEXUS files containing character state labels
#                - Extracts meaning groups from Indo-European language data
#                - Extracts meaning groups from Sino-Tibetan language data
#                - Outputs organized meaning sets as CSV files for downstream analysis
# ------------------------------------------------------------------------------
library(ape)
library(here)

compute_meaning_set_ie <- function(input_file, output_file) {
  # Read the nexus file to get character labels
  nexus_content <- scan(input_file, what='character', sep='\n', quiet=TRUE)  
  
  # Find the charstatelabels section
  start_idx <- which(grepl("charstatelabels", nexus_content))
  if (length(start_idx) == 0) {
    stop("Could not find charstatelabels section in nexus file")
  }
  
  # Extract character labels
  char_labels <- c()
  for (i in (start_idx + 1):length(nexus_content)) {
    line <- trimws(nexus_content[i])
    if (line == ";" || line == "") break
    
    # Parse lines like "1 ant_group," or "2 ant_cognate_5007,"
    if (grepl("^[0-9]+ ", line)) {
      parts <- strsplit(line, " ")[[1]]
      if (length(parts) >= 2) {
        label <- gsub(",$", "", parts[2])  # Remove trailing comma
        char_labels <- c(char_labels, label)
      }
    }
  }
  
  # Extract meaning groups from character labels
  meaning_groups <- data.frame(
    meaning = character(),
    start = integer(),
    end = integer(),
    stringsAsFactors = FALSE
  )
  
  current_meaning <- ""
  start_pos <- 1
  
  for (i in 1:length(char_labels)) {
    label <- char_labels[i]
    
    # Extract meaning from character label (format: meaning_group)
    if (grepl("_group$", label)) {
      # This is a meaning group marker
      meaning_name <- gsub("_group$", "", label)
      
      # If we have a previous meaning, record its end position
      if (current_meaning != "") {
        meaning_groups <- rbind(meaning_groups, data.frame(
          meaning = current_meaning,
          start = start_pos,
          end = i - 1,
          stringsAsFactors = FALSE
        ))
      }
      
      # Start new meaning group
      current_meaning <- meaning_name
      start_pos <- i
    }
  }
  
  # Add the last meaning group
  if (current_meaning != "") {
    meaning_groups <- rbind(meaning_groups, data.frame(
      meaning = current_meaning,
      start = start_pos,
      end = length(char_labels),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by start position
  meaning_groups <- meaning_groups[order(meaning_groups$start), ]
  
  return(meaning_groups)
}

compute_meaning_set_st <- function(input_file, output_file) {
  # Read the nexus file to get character labels
  nexus_content <- scan(input_file, what='character', sep='\n', quiet=TRUE)  
  
  # Find the charstatelabels section
  start_idx <- which(grepl("CHARSTATELABELS", nexus_content))
  if (length(start_idx) == 0) {
    stop("Could not find charstatelabels section in nexus file")
  }
  
  # Extract character labels and their positions
  char_labels <- c()
  char_positions <- c()
  for (i in (start_idx + 1):length(nexus_content)) {
    line <- trimws(nexus_content[i])
    if (line == ";" || line == "") break
    
    # Parse lines like "1 I_first_person_singular_ascertainment,"
    if (grepl("^[0-9]+ ", line)) {
      parts <- strsplit(line, " ")[[1]]
      if (length(parts) >= 2) {
        position <- as.integer(parts[1])
        label <- gsub(",$", "", parts[2])  # Remove trailing comma
        char_labels <- c(char_labels, label)
        char_positions <- c(char_positions, position)
      }
    }
  }
  
  # Extract meaning groups from character labels
  meaning_groups <- data.frame(
    meaning = character(),
    start = integer(),
    end = integer(),
    stringsAsFactors = FALSE
  )
  
  current_meaning <- ""
  start_pos <- 1
  
  for (i in 1:length(char_labels)) {
    label <- char_labels[i]
    
    # Check if this is an ascertainment marker (new meaning group)
    if (grepl("_ascertainment$", label)) {
      # If we have a previous meaning, record its end position
      if (current_meaning != "") {
        meaning_groups <- rbind(meaning_groups, data.frame(
          meaning = current_meaning,
          start = start_pos,
          end = char_positions[i - 1],
          stringsAsFactors = FALSE
        ))
      }
      
      # Start new meaning group
      current_meaning <- gsub("_ascertainment$", "", label)
      start_pos <- char_positions[i]
    }
  }
  
  # Add the last meaning group
  if (current_meaning != "") {
    meaning_groups <- rbind(meaning_groups, data.frame(
      meaning = current_meaning,
      start = start_pos,
      end = char_positions[length(char_positions)],
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by start position
  meaning_groups <- meaning_groups[order(meaning_groups$start), ]
  
  return(meaning_groups)
}

# Main execution
# indo-european languages
input_file_ie <- here("data/real/iecor_ctmc-strict-M1/iecor.nex")
output_file_ie <- here("output/results/meanings_sets_ie.csv")

meaning_groups_ie <- compute_meaning_set_ie(input_file_ie, output_file_ie)
write.csv(meaning_groups_ie, output_file_ie, row.names = FALSE, quote = FALSE)

# sino-tibetan languages
input_file_st <- here("data/real/st_ctmc-strict-fbd-uni/st.nex")
output_file_st <- here("output/results/meanings_sets_st.csv")

meaning_groups_st <- compute_meaning_set_st(input_file_st, output_file_st)
write.csv(meaning_groups_st, output_file_st, row.names = FALSE, quote = FALSE)


