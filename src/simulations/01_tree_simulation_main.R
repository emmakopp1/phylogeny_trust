# ------------------------------------------------------------------------------
# Script Name: 01_tree_simulation.R
# Description: This script simulates phylogenetic trees and generates input XML 
#              files for BEAST 2.6.1 analyses. It performs the following steps:
#                - Simulates birth-death trees with calibrations
#                - Scales branch lengths to represent different ages
#                - Prepares BEAST XML files for sequence simulation
#                - Runs BEAST to simulate sequences
#                - Embeds simulated sequences into BEAST analysis XML files
#                - Applies node calibrations and modifies prior distributions
#                - Updates output file names for BEAST runs
# -----------------------------------------------------------------------------------------
library(here)
library(TreeSim)
library(ape)
library(purrr)
library(stringr)
library(xml2)
library(readr)
library(magrittr)
library(dplyr)
library(beastier)

# simulation of the initial tree --------------------------------------------
path <- here(sprintf("data/simulated_temp/beast-data-sim-%s", Sys.Date()))
dir_path <-sprintf("data/simulated-%s", Sys.Date())
unlink(dir_path, recursive = TRUE, force = T)
dir.create(here(dir_path))

# parameter initialization
N <- 50
div <- 0.245
r <- 0.293

pi1 <- 5.695E-2
pi0 <- 1 - pi1


# number of simulation per age
N_sim <- 50
# number of of different ages
N_rep <- 17
#N_trait <- round(150/pi1,0)
N_trait <- 3000


# create N_sim folders 
for (n_sim in 1:N_sim){
  dir.create(paste0(dir_path, sprintf("/beast-data-sim-%d", n_sim))) 
}


phylogeny <- read.nexus(here("data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees"))
M <- length(phylogeny)
burnin <- 0.1
phylogeny <- phylogeny[seq(burnin * M, M)]  
index <- sample(1:length(phylogeny), N_sim)
tree <- phylogeny[index]

# function which transform fossils as tips
standardize_to_max_depth <- function(tree) {
  # Compute depths and max in one pass
  tip_depths <- node.depth.edgelength(tree)[1:Ntip(tree)]
  max_depth <- max(tip_depths)

  # Vectorized computation of differences
  diff_to_add <- max_depth - tip_depths

  # Identify terminal edges (tips are always 1:Ntip in edge[,2])
  terminal_edges <- match(1:Ntip(tree), tree$edge[, 2])

  # Vectorized adjustment
  tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges] + diff_to_add
  
  return(tree)
}

tree = lapply(tree, function(t) standardize_to_max_depth(t))


# 3 calibrations, burmish, sinitic and tibetan
calib_chinese <- tree[[1]]$tip.label[grep('Sinitic',tree[[1]]$tip.label)]
calib_tibetan <- tree[[1]]$tip.label[grep('Tibetan',tree[[1]]$tip.label)]
calib_burmish <- c("BurmishOldBurmese","BurmishRangoon")
calib <- list(calib_chinese, calib_tibetan, calib_burmish)


# scaling ----------------------------------------------------------------------
l <- seq(1, N_rep, 1)
#l <- seq(N_rep, N_rep, 1)

for (i in 1:N_sim){
  tree_i <- tree[[i]]
  t_r <- max(node.depth.edgelength(tree_i)) 
  for (k in l){
    new_tree <- tree_i
    new_tree$edge.length <- k * as.numeric(new_tree$edge.length) / t_r
    dir.create(paste0(dir_path, sprintf("/beast-data-sim-%d/beast-data-sim-%d-%d",i, i, k)))
    write.tree(
      new_tree, 
      file = paste0(
        dir_path, 
        sprintf("/beast-data-sim-%d/beast-data-sim-%d-%d/tree-sim-%d-%d.tree", i, i, k, i, k))
    )
  }
}


# generate sequences -----------------------------------------------------------
# check the number of trait is accurate with N_trait
check_trait_in_template <- function(path_simulation_template, N_trait){
  # read the file: simulation_template
  doc <- read_xml(path_simulation_template)
  run_node <- xml_find_first(doc, ".//run")
  N_traits_in_file <- xml_attr(run_node, "sequencelength") |> as.integer()
  
  # modify the value
  if (N_traits_in_file != N_trait) {
    print("Number of traits changed in the template")
    xml_set_attr(run_node, "sequencelength", as.character(N_trait))
    write_xml(doc, path_simulation_template)
  }
}

check_trait_in_template(here("data/beast-data-sim.xml"), N_trait)

files <- list.files(dir_path, full.names = TRUE, recursive = TRUE)
target_text <- read_lines(here("data/beast-data-sim.xml")) %>% paste(collapse = "\n")
process_file <- function(file) {
  # get the tree
  file_text <- read_lines(file) |> paste(collapse = "\n")
  # age of the tree 
  tree_age <-  as.numeric(str_extract(file, "(\\d+)(?=\\.tree)"))
  tree_simulation_number <- as.numeric(str_match(file, "beast-data-sim-(\\d+)-\\d+")[,2])
  # replace xyz by tree & replace output_name by beast-simulated-seq-etc
  str_replace(target_text, "xyz", file_text) |>
    str_replace("output_name", sprintf("beast-simulated-seq-%d-%d.xml",tree_simulation_number, tree_age)) |> 
    write_lines(
      file |>
        str_replace("/[^/]*$", sprintf("/beast-data-sim-%d-%d.xml",tree_simulation_number, tree_age))
    )}


updated_texts <- files |>
  map_chr(~ process_file(.x))



# generate sequence with beast -------------------------------------------------
for(n_sim in 1:N_sim){
  for (i in  l) {
    # run beast to generate sequence
    system(paste0("../../../../Applications/BEAST2.6.7/bin/beast -overwrite ", 
                  dir_path, sprintf("/beast-data-sim-%d/beast-data-sim-%d-%d/beast-data-sim-%d-%d.xml", n_sim, n_sim, i, n_sim, i)))
    # change the emplacement of the output
    system(
      paste0(sprintf("mv beast-simulated-seq-%d-%d.xml ",n_sim, i),
             "/Users/kopp/Documents/phylogeny_trust/", 
             dir_path, 
             sprintf("/beast-data-sim-%d/beast-data-sim-%d-%d/beast-simulated-seq-%d-%d.xml", n_sim, n_sim, i, n_sim, i))
      )
  }
}

# generate xml files  ----------------------------------------------------------
# create ctmc-strict-bd-n_sim-1.xml for each simulation index
for (i in 1:N_sim){
  #cat('simu:',i, '\n')
  path_template_beauti <- here('data/ctmc-strict-bd-template.xml')
  beauti_template <- read_xml(path_template_beauti)
  
  tree_i <- read.tree(here(paste0(
    dir_path, 
    sprintf("/beast-data-sim-%d/beast-data-sim-%d-1/tree-sim-%d-1.tree", i, i, i))))
  
  # get the mrca node of calibrations for simulation i
  mrca_node <- list()
  for (j in 1:length(calib)){
    cal = unlist(calib[j])
    calib_tip <- match(cal, tree_i$tip.label)
    mrca_cal <- getMRCA(tree_i, calib_tip)
    mrca_node <- c(mrca_node,mrca_cal)
  }
  mrca_node <- unlist(mrca_node)
  
  # get ages
  t_R <- max(node.depth.edgelength(tree_i))
  t <- t_R - node.depth.edgelength(tree_i)[mrca_node]
  
  # set prior parameters
  calib_inf = sapply(t, function(x) max(0, x - 0.01))
  calib_sup = sapply(t, function(x) min(t_R, x + 0.01))
  
  # in the first folder : beast-data-sim-1
  # change the taxas of the calibration 
  calib_node_chinese <- beauti_template |> xml_find_all("//distribution[contains(@id, 'chinese.prior')]")
  calib_node_tibetan <- beauti_template |> xml_find_all("//distribution[contains(@id, 'tibetan.prior')]")
  calib_node_burmish <- beauti_template |> xml_find_all("//distribution[contains(@id, 'burmish.prior')]")
  
  # impose monophylecy
  xml_set_attr(calib_node_chinese[[1]], "monophyletic", "true")
  xml_set_attr(calib_node_tibetan[[1]], "monophyletic", "true")
  xml_set_attr(calib_node_burmish[[1]], "monophyletic", "true")
  
  # add the uniform prior calibrations
  # for the chinese
  prior_calibration_node_chinese = calib_node_chinese[[1]] |> xml_child(2) 
  xml_set_attr(prior_calibration_node_chinese, "lower", round(calib_inf,2)[1])
  xml_set_attr(prior_calibration_node_chinese, "upper", round(calib_sup,2)[1])
  
  # for tibetan
  prior_calibration_node_tibetan = calib_node_tibetan[[1]] |> xml_child(2) 
  xml_set_attr(prior_calibration_node_tibetan, "lower", round(calib_inf,2)[2])
  xml_set_attr(prior_calibration_node_tibetan, "upper", round(calib_sup,2)[2])
  
  # for burmish
  prior_calibration_node_burmish = calib_node_burmish[[1]] |> xml_child(2) 
  xml_set_attr(prior_calibration_node_burmish, "lower", round(calib_inf,2)[3])
  xml_set_attr(prior_calibration_node_burmish, "upper", round(calib_sup,2)[3])
  
  # write the xml
  output_path <- here(paste0(
    dir_path, 
    sprintf("/beast-data-sim-%d/beast-data-sim-%d-1/ctmc-strict-bd-%d-1.xml", i, i, i)
  ))
  
  write_xml(beauti_template, 
            output_path, 
            options = "format")
}


# Generate BEAST file for analysis ---------------------------------------------
sequences <- list.files(dir_path, full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.xml$")) |>
  keep(~ str_detect(.x, "simulated-seq"))


beast_files <- lapply(1:N_sim, function(i) {
  file_path <- here(paste0(
    dir_path, 
    sprintf("/beast-data-sim-%d/beast-data-sim-%d-1/ctmc-strict-bd-%d-1.xml", i, i, i)))
  read_xml(file_path)
})

# replace the sequence values of each taxon
replace_value <- function(xml_file, df, path_out) {
  xml_file |>
    xml_find_all("//sequence") |>
    walk(~ {
      taxa <- xml_attr(.x, "taxon")
      new_value <- df |>
        filter(taxon == taxa) |>
        pull(value)
      if (length(new_value) > 0) {
        xml_set_attr(.x, "value", new_value)
      }
    })
  
  # Save the modified XML file
  write_xml(xml_file, path_out)
}

# apply function replace_value
for (file in sequences) {
    
  # extract sequences
  taxons <- read_xml(file) |>
    xml_find_all("//sequence") |>
    map_df(~ {
      taxon <- xml_attr(.x, "taxon")
      value <- xml_attr(.x, "value")
      data.frame(
        taxon = taxon, value = value
      )
    })
    
  # output path of the xml file
  path_out <- str_replace(
    file,
    "beast-simulated-seq",
    paste0("ctmc-strict-bd")
  )
  
  simulation_number <- as.numeric(str_match(file, "beast-data-sim-(\\d+)-\\d+")[,2])
  
  replace_value(beast_files[[simulation_number]], taxons, path_out)
}

# Modify calibrations ----------------------------------------------------------
beast_inputs = list.files(dir_path, full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.xml$")) |>
  keep(~ str_detect(.x, "ctmc"))

# modify the prior calibration of a file. The coefficient is in the file_path name
modify_uniform_attributes <- function(file_path) {
  
  tree_age <- as.integer(str_extract(file_path, "(?<=-)(\\d+)(?=\\.xml)"))
  xml_file <- read_xml(file_path)
  uniform_node <- xml_file |>
    xml_find_all("//distribution[contains(@id, '.prior')]//Uniform")
  
  # Modify lower
  current_lower <- as.numeric(xml_attr(uniform_node, "lower"))
  new_lower <- current_lower * tree_age
  xml_set_attr(uniform_node, "lower", as.character(new_lower))

  # Modify upper
  current_upper <- as.numeric(xml_attr(uniform_node, "upper"))
  new_upper <- current_upper * tree_age
  xml_set_attr(uniform_node, "upper", as.character(new_upper))
  
  write_xml(xml_file, file_path)
}

# change the .trees and .log files (fileName attribute)
modify_filenames <- function(file_path) {
  
  tree_age <- as.integer(str_extract(file_path, "(?<=-)(\\d+)(?=\\.xml)"))
  xml_file <- read_xml(file_path)
  
  
  xml_find_all(xml_file, "//logger[contains(@id, 'treelog') or contains(@id, 'tracelog')]") |>
    walk(~ {
      current_fileName <- xml_attr(.x, "fileName")
      new_fileName <- str_replace(current_fileName, "\\d+", as.character(tree_age))
      #cat("Changing", current_fileName, "to", new_fileName, "\n")
      xml_set_attr(.x, "fileName", new_fileName)
    })
  
  write_xml(xml_file, file_path)
}

# Apply functions
walk(beast_inputs, modify_uniform_attributes)
walk(beast_inputs, modify_filenames)
 



