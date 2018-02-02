library(ape)
library(phytools)
library(geiger)
library(TreeSim)
source("scripts/functions.R")


n_tips <- c(10, 50, 100, 200, 400)
rates <- c(1,2,5,10,100,1000)
no_sim <- 50
birth = 1.5
death = 0.5
nsim=1

for (i in 1:no_sim){
  for (n in n_tips){
    sample <- sim_subclade(n, n/5, birth, death)
    tree <- sample[[1]]
    subtree <- sample[[2]]
    name = paste0("t", i, "_N",n, "_clade", n/5, "_b", birth, "_d", death)
    write.nexus(tree, file = paste0("data/", name, ".nex"))
    write.tree(tree, file = paste0("data/", name, ".nwk"))
    for (r in rates){
      trait <- simBM_2rates(tree, subtree, R1=1, R2=r, model = "BM", root = 1)
      write.table(trait, file = paste0("data/", name, "_rate", r, ".trait"),
                  quote = F, row.names = T, col.names = F)
    }
  }
}
