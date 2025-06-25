# pour un seul path
tt = read.tree(here("data/simulated-2025-05-13/beast-data-sim-1/beast-data-sim-1-10/tree-sim-1-10.tree"))
posterior = read.nexus(here("data/simulated-2025-05-13/beast-data-sim-1/beast-data-sim-1-10/ctmc-strict-bd-10.trees"))
mcc = read.tree(here("data/simulated-2025-05-13/beast-data-sim-1/beast-data-sim-1-14/mcc-14.tree"))

cs = read.tree(here("data/simulated-2025-05-13/beast-data-sim-1/beast-data-sim-1-17/consensus-17.tree"))

plot(cs)
nodelabels(cex=0.8, frame='circle')

prob_first_split_cs = prob_first_split_resumed |> 
  select(- mcc_prob, - node_mcc, - node_cs) |> 
  group_by(age) |> 
  summarise(mean_cs_prob = mean(cs_prob, na.rm=T), .groups='drop') |> 
  ungroup()

tt = prob_first_split_resumed |> 
  filter(age == 17) |> 
  select(-mcc_prob, -node_mcc, -age)

prob_first_split_cs
