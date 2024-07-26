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


# 1. Simulation of the initial tree
path <- here("data/simulated_temp/beast-data-sim-")
dir.create(here("data/simulated_temp"))

# Parameter initialization
N <- 50
div <- 0.245
r <- 0.293
s <- 0.226

lambda <- div / (1 - r)
mu <- (r * div) / (1 - r)
psi <- (s / (1 - s)) * (r * div / (1 - r))
t <- 1

# Tree simulation
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


# 2. Construction of the scaled trees
l <- seq(1, 17, 1)

for (k in l) {
  new_tree <- tree
  new_tree$edge.length <- k * as.numeric(new_tree$edge.length)
  dir.create(sprintf(here("data/simulated_temp/beast-data-sim-%d"), k))
  write.tree(new_tree,
    file = sprintf(here("data/simulated_temp/beast-data-sim-%d/tree-sim-%d.tree"), k, k)
  )
}

# Generate sequences
files <- list.files(here("data/simulated_temp"), full.names = TRUE, recursive = TRUE)
target_text <- read_lines(here("data/beast-data-sim.xml")) %>% paste(collapse = "\n")

process_file <- function(file) {
  file_text <- read_lines(file) |> paste(collapse = "\n")
  str_replace(target_text, "xyz", file_text) |>
    str_replace("output_name", sprintf("beast-simulated-seq-%d.xml", as.numeric(str_extract(file, "(\\d+)(?=\\.tree)")))) |>
    write_lines(
      file %>%
        str_replace("/tree-sim-\\d+\\.tree$", "") %>% # Remove the tree file part
        str_replace("beast-data-sim-\\d+", "\\0/\\0.xml")
    )
}

updated_texts <- files |>
  map_chr(~ process_file(.x))

for (i in 1:17) {
  # Run beast to generate sequence
  system(sprintf(("../../../../Applications/beast/bin/beast -overwrite /Users/kopp/Documents/phylogeny_trust/data/simulated_temp/beast-data-sim-%d/beast-data-sim-%d.xml"), i, i))
  # Change the emplacement of the output
  system(sprintf(("mv beast-simulated-seq-%d.xml  data/simulated_temp/beast-data-sim-%d/beast-simulated-seq-%d.xml"), i, i, i))
}


# Generate BEAST file for analysis
target_text <- read_lines(here("data/ctmc-strict-bd.xml")) |> paste(collapse = "\n")

sequences <- list.files(here("data/simulated_temp"), full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.xml$")) |>
  keep(~ str_detect(.x, "simulated-seq"))

beast_file <- read_xml(here("data/ctmc-strict-bd.xml"))

# Replace in .xml file the values of the generated sequences
replace_value <- function(xml_file, df, path_out) {
  xml_file |>
    xml_find_all("//sequence") |>
    walk(~ {
      taxon <- xml_attr(.x, "taxon")
      new_value <- df |>
        filter(taxon == taxon) |>
        pull(value)
      if (length(new_value) > 0) {
        xml_set_attr(.x, "value", new_value)
      }
    })
  
  # Sauvegarder le fichier XML modifié
  write_xml(xml_file, path_out)
}

# Appliquer la fonction
for (file in sequences) {
  taxons <- read_xml(file) |>
    xml_find_all("//sequence") |>
    map_df(~ {
      taxon <- xml_attr(.x, "taxon")
      value <- xml_attr(.x, "value")
      data.frame(
        taxon = taxon, value = value
      )
    })

  path_out <- str_replace(
    file,
    "beast-simulated-seq-\\d+.xml",
    paste0("ctmc-strict-bd-", str_extract(file, "\\d+"), ".xml")
  )


  replace_value(beast_file, taxons, path_out)
}

# Modify calibrations
beast_inputs = list.files(here("data/simulated_temp"), full.names = TRUE, recursive = TRUE) |>
  keep(~ str_detect(.x, "\\.xml$")) |>
  keep(~ str_detect(.x, "ctmc"))


modify_uniform_attributes <- function(file_path) {
  
  scalar = as.integer(str_extract(file_path, "(?<=-)(\\d+)(?=\\.xml)"))
  xml_file <- read_xml(file_path)
  uniform_node <- read_xml(file_path) |>
    xml_find_all("//distribution[contains(@id, '.prior')]//Uniform")

  # Modify lower
  current_lower <- as.numeric(xml_attr(uniform_node, "lower"))
  new_lower <- current_lower * scalar
  xml_set_attr(uniform_node, "lower", as.character(new_lower))
  
  # Modify upper
  current_upper <- as.numeric(xml_attr(uniform_node, "upper"))
  new_upper <- current_upper * scalar
  xml_set_attr(uniform_node, "upper", as.character(new_upper))

  write_xml(xml_file, file_path)
}


modify_filenames <- function(file_path) {
  
  scalar <- as.integer(str_extract(file_path, "(?<=-)(\\d+)(?=\\.xml)"))
  xml_file <- read_xml(file_path)
  
  # Modify file names
  xml_find_all(xml_file, "//logger[contains(@id, 'treelog') or contains(@id, 'tracelog')]") |>
    walk(~ {
      current_fileName <- xml_attr(.x, "fileName")
      new_fileName <- str_replace(current_fileName, "\\d+", as.character(scalar))
      xml_set_attr(.x, "fileName", new_fileName)
    })
  
  # Save 
  write_xml(xml_file, file_path)
}


# Apply functions
walk(beast_inputs, modify_filenames)



# BEAST analysis 
for (i in 1:17) {
  # Run beast to generate sequence
  system(
    sprintf(("../../../../Applications/beast/bin/beast -overwrite /Users/kopp/Documents/phylogeny_trust/data/simulated_temp/beast-data-sim-%d/ctmc-strict-bd-%d.xml"), i, i),
    ignore.stdout = TRUE, 
    ignore.stderr = TRUE)
  system(sprintf(("mv ctmc-strict-bd-%d.log  data/simulated_temp/beast-data-sim-%d/ctmc-strict-bd-%d.log"), i, i, i))
  system(sprintf(("mv ctmc-strict-bd-%d.trees  data/simulated_temp/beast-data-sim-%d/ctmc-strict-bd-%d.trees"), i, i, i))
  system(sprintf(("rm ctmc-strict-bd-%d.xml.0state"), i))
  system(sprintf(("rm ctmc-strict-bd-%d.xml.1state"), i))
  system(sprintf(("rm ctmc-strict-bd-%d.xml.2state"), i))
  system(sprintf(("rm ctmc-strict-bd-%d.xml.3state"), i))
  system(sprintf(("rm ctmc-strict-bd-%d.xml.state"), i))
}











