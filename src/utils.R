library(here)
library(dplyr)
library(tibble)
library(ape)
library(TreeTools)
library(phytools)
library(adephylo)
library(castor)

remove_burnin = function(trees,burnin_rate){
  n = as.numeric(length(trees))
  return(trees[as.integer(n*burnin_rate):n])
}

myconsensus = function(trees){
  consensus.edges(trees, consensus.tree = consensus(trees, p=.5), rooted=T)    
}


path= here("data/real/tea_ctmc-strict-fbd-constrained/tea_ctmc-strict-fbd-constrained.trees")
trees = read.nexus(path)
trees = remove_burnin(trees,0.2)
M = length(trees)
tree = trees[[M]]



# Age of MRCA japonic 
mrca_japonic = getMRCA(tree,tea_group_japonic$tip)
age_mrca_japonic = (distRoot(tree,1)[[1]] - distRoot(tree,mrca_japonic)[[1]])*0.1

# Age of MRCA koreanic
mrca_koreanic = getMRCA(tree,tea_group_koreanic$tip)
age_mrca_koreanic = (distRoot(tree,1)[[1]] - distRoot(tree,mrca_koreanic)[[1]])*0.1


# Age of MRCA turkic
mrca_turkic = getMRCA(tree,tea_group_turkic$tip)
age_mrca_turkic = (distRoot(tree,1)[[1]] - distRoot(tree,mrca_turkic)[[1]])*0.1

# Age of MRCA mongolian
mrca_mongolian = getMRCA(tree,tea_group_mongolian$tip)
age_mrca_mongolian = (distRoot(tree,1)[[1]] - distRoot(tree,mrca_mongolian)[[1]])*0.1

# Age of MRCA tungusic
mrca_tungusic = getMRCA(tree,tea_group_tungusic$tip)
age_mrca_tungusic = (distRoot(tree,1)[[1]] - distRoot(tree,mrca_tungusic)[[1]])*0.1



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


path= here("data/real/st_ctmc-strict-fbd/st_ctmc-strict-fbd.trees")
trees = read.nexus(path)
trees = remove_burnin(trees,0.99)

c = mean(sapply(trees, function(tree) distRoot(tree)[[1]]))


