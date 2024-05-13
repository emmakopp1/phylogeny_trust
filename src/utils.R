library(here)
library(dplyr)
library(tibble)
source(here("src/init.R"))

remove_burnin = function(trees,burnin_rate){
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}

myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}

rwty::load.trees(file,trim=100)

path= here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees")
trees = read.nexus(path)
M = length(trees)


fossils_st = c("BurmishOldBurmese","Tangut","SiniticOldChinese","TibetanOldTibetan")
fossils_iecor = c(
  "Hittite",
  "Luvian",
  "TocharianA",
  "TocharianB",
  "MycenaeanGreek",
  "AncientGreek",
  "NTGreek",
  "TsakonianPropontis",
  "ClassicalArmenian",
  "VedicEarly",
  "Pali",
  "AvestanYounger",
  "Khwarazmian",
  "Sogdian",
  "Khotanese",
  "Bactrian",
  "Parthian",
  "OldPersian",
  "MiddlePersian",
  "OldPrussian",
  "OldChurchSlavonic",
  "SloveneEarlyModern",
  "OldPolish",
  "Polabian",
  "OldCzech6",
  "OldNovgorod5",
  "Gothic",
  "OldIcelandic",
  "OldSwedish",
  "OldSaxon5",
  "OldHighGerman",
  "MiddleHighGerman",
  "MiddleDutch",
  "OldFrisian5",
  "OldEnglish",
  "Latin",
  "DalmatianVegliote",
  "AngloNorman",
  "OldOccitan",
  "OldCatalan",
  "OldSpanish",
  "Oscan",
  "Umbrian",
  "Gaulish",
  "OldWelsh",
  "MiddleWelsh",
  "MiddleCornish",
  "LateCornish",
  "OldBreton",
  "MiddleBreton",
  "OldIrish",
  "GaelicManx"
  )


compute_min_age_fossil = function(trees,fossils){
  M = length(trees)
  min(rep(distRoot(trees[[M]],1),length(fossils)) - distRoot(trees[[M]],fossils))
}

df = tibble(
  age= compute_min_age_fossil(trees,fossils_st)
    )



