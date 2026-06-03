library(phytools)
library(xml2)

library(ape)
library(phangorn)

trf = "/Users/kopp/Documents/phylogeny_trust/data/simulated-2025-07-28/beast-data-sim-1/beast-data-sim-1-1/tree-sim-1-1.tree"
tree = read.tree(trf)

datf = "/Users/kopp/Documents/phylogeny_trust/data/simulated-2025-07-28/beast-data-sim-1/beast-data-sim-1-1/beast-simulated-seq-1-1.xml"
dat = read_xml(datf)

seq_nodes <- xml_find_all(dat, ".//sequence")

taxa <- xml_attr(seq_nodes, "taxon")
seqs  <- xml_attr(seq_nodes, "value")

names(seqs) <- taxa
seqs_ordered <- seqs[match(tree$tip.label, taxa)]
names(seqs_ordered) <- tree$tip.label

# la racine a 2 enfants, qui donnent des sous-arbres A et B
root_children <- tree$edge[tree$edge[,1] == Ntip(tree) + 1, 2]


tips_A <- tree$tip.label[
  Descendants(tree, root_children[1], type = "tips")[[1]]
]

tips_B <- tree$tip.label[
  Descendants(tree, root_children[2], type = "tips")[[1]]
]

# sequences from each subtree
seqs_A <- seqs_ordered[tips_A]
seqs_B <- seqs_ordered[tips_B]


mat_A <- do.call(
  rbind,
  strsplit(seqs_A, "")
)

mat_B <- do.call(
  rbind,
  strsplit(seqs_B, "")
)

mat <- rbind(mat_A, mat_B)
N = sum(colSums(mat=="0") == ncol(mat))

present_A <- colSums(mat_A == "1") > 0
present_B <- colSums(mat_B == "1") > 0

shared_sites <- present_A & present_B

sum(shared_sites) / length(present_A) 


