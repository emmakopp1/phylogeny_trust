rm(list=ls())
library(here)
library(TreeSim)
library(ape)
library(purrr)
library(stringr)
library(readr)
library(phangorn)
library(xml2)
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
s <- 0.226

lambda <- div / (1 - r)
mu <- (r * div) / (1 - r)
psi <- (s / (1 - s)) * (r * div / (1 - r))
t <- 1

# tree simulation
tree <- sim.bd.taxa.age(
  n = N,
  numbsim = 1,
  lambda = lambda,
  mu = mu,
  frac = 1,
  age = t,
  mrca = TRUE
)

tree <- tree[[1]]

# scaling ----------------------------------------------------------------------
l <- seq(1, 17, 1)

for (k in l) {
  new_tree <- tree
  new_tree$edge.length <- k * as.numeric(new_tree$edge.length)
  dir.create(paste0(dir_path, sprintf("/beast-data-sim-%d",k)))
  write.tree(
    new_tree,
    file = paste0(dir_path, sprintf("/beast-data-sim-%d/tree-sim-%d.tree", k, k))
  )
}

# generate sequences -----------------------------------------------------------
files <- list.files(dir_path, full.names = TRUE, recursive = TRUE)
target_text <- read_lines(here("data/beast-data-sim.xml")) %>% paste(collapse = "\n")

process_file <- function(file) {
  # get the tree
  file_text <- read_lines(file) |> paste(collapse = "\n")
  str_replace(target_text, "xyz", file_text) |>
    str_replace("output_name", sprintf("beast-simulated-seq-%d.xml", as.numeric(str_extract(file, "(\\d+)(?=\\.tree)")))) |>
    write_lines(
      file |>
        str_replace("/tree-sim-\\d+\\.tree$", "") %>% # Remove the tree file part
        str_replace("beast-data-sim-\\d+", "\\0/\\0.xml")
    )}

updated_texts <- files |>
  map_chr(~ process_file(.x))


# generate sequence with beast -------------------------------------------------
for (i in 1:17) {
  # run beast to generate sequence
  system(paste0("../../../../Applications/BEAST2.6.7/bin/beast -overwrite ", dir_path, sprintf("/beast-data-sim-%d/beast-data-sim-%d.xml", i, i)))
  # change the emplacement of the output
  system(
    paste0(sprintf("mv beast-simulated-seq-%d.xml ",i), "/Users/kopp/Documents/phylogeny_trust/", dir_path, sprintf("/beast-data-sim-%d/beast-simulated-seq-%d.xml",i,i))
    )
}

# GENERATE MANUALLY A BEAST FILE


# template ctmc
# import alignement
# calibration
#tt <- read.tree("/Users/kopp/Documents/phylogeny_trust/data/simulated-2025-03-17/beast-data-sim-1/tree-sim-1.tree")
#calib <- c('t22','t17','t19','t23')
#max(node.depth.edgelength(tt)) - node.depth.edgelength(tt)[findMRCA(tt,calib)]


# Generate BEAST file for analysis ---------------------------------------------
sequences <- list.files(dir_path, full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.xml$")) |>
  keep(~ str_detect(.x, "simulated-seq"))

beast_file <- read_xml(here(paste0(dir_path,"/beast-data-sim-1/ctmc-strict-bd-1.xml")))

# replace in .xml file the values of the generated sequences
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
  
  # Sauvegarder le fichier XML modifié
  write_xml(xml_file, path_out)
}

# apply function 
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
    "beast-simulated-seq-\\d+.xml",
    paste0("ctmc-strict-bd-", str_extract(file, "\\d+(?=\\.)"), ".xml")
  )

  replace_value(beast_file, taxons, path_out)
}

# Modify calibrations ----------------------------------------------------------
beast_inputs = list.files(dir_path, full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.xml$")) |>
  keep(~ str_detect(.x, "ctmc"))

# modify the prior calibration of a file. The coeficiant is in the file_path name
modify_uniform_attributes <- function(file_path) {
  
  scalar <- as.integer(str_extract(file_path, "(?<=-)(\\d+)(?=\\.xml)"))
  xml_file <- read_xml(file_path)
  uniform_node <- xml_file |>
    xml_find_all("//distribution[contains(@id, '.prior')]//Uniform")
  
  # Modifier lower
  current_lower <- as.numeric(xml_attr(uniform_node, "lower"))
  new_lower <- current_lower * scalar
  xml_set_attr(uniform_node, "lower", as.character(new_lower))
  
  # Modifier upper
  current_upper <- as.numeric(xml_attr(uniform_node, "upper"))
  new_upper <- current_upper * scalar
  xml_set_attr(uniform_node, "upper", as.character(new_upper))
  
  write_xml(xml_file, file_path)
}

# change the .trees and .log files (fileName attribute)
modify_filenames <- function(file_path) {
  scalar <- as.integer(str_extract(file_path, "(?<=-)(\\d+)(?=\\.xml)"))
  xml_file <- read_xml(file_path)
  
  
  xml_find_all(xml_file, "//logger[contains(@id, 'treelog') or contains(@id, 'tracelog')]") |>
    walk(~ {
      current_fileName <- xml_attr(.x, "fileName")
      new_fileName <- str_replace(current_fileName, "\\d+", as.character(scalar))
      #cat("Changing", current_fileName, "to", new_fileName, "\n")
      xml_set_attr(.x, "fileName", new_fileName)
    })
  
  write_xml(xml_file, file_path)
  #cat("File saved:", file_path, "\n")
}

# Apply functions
walk(beast_inputs, modify_uniform_attributes)
walk(beast_inputs, modify_filenames)

# Just need to lauch the files 












