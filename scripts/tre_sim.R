library(ape)
library(phytools)
library(geiger)
library(TreeSim)
source("scripts/functions.R")

# script used in snakemake pipeline use the commented block if running without snakemake
# n_tips <- c(10, 50, 100, 200, 400)
# rates <- c(1,2,5,10,100,1000)
# n_sim <- 1
# birth = 1.5
# death = 0.5

n_sim  = snakemake@config[["n_sim"]]
n_tips = snakemake@config[["n_tips"]]
rates  = snakemake@config[["rates"]]
birth  = snakemake@config[["birth"]]
death  = snakemake@config[["death"]]
clade_size = snakemake@config[["clade_size"]]


for (i in 1:n_sim){
  for (n in n_tips){
    sample <- sim_subclade(n, n * clade_size, birth, death)
    tree <- sample[[1]]
    subtree <- sample[[2]]
    name = paste0("t", i, "_N",n, "_clade", n * clade_size, "_b", birth, "_d", death)
    write.nexus(tree, file = paste0("data/", name, ".nex"))
    write.tree(tree, file = paste0("data/", name, ".nwk"))
    for (r in rates){
      trait <- simBM_2rates(tree, subtree, R1=1, R2=r, model = "BM", root = 1)
      write.table(trait, file = paste0("data/", name, "_rate", r, ".trait"),
                  quote = F, row.names = T, col.names = F)
    }
  }
}

write.table(n_tips, "out/simulated_trees_info.tsv")

