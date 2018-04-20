library(ape)
library(phytools)
library(geiger)
library(TreeSim)
source("scripts/functions.R")

# script used in snakemake pipeline. use the commented block if running without snakemake

################################################################################
# Standalone settings

# n_tips <- c(10,50, 100, 200, 500)
# rates <- c(2,10,100)
# n_sim <- 1
# birth = 1.5
# death = 0.5
# clade_size = 0.2

################################################################################
# Settings for SNAKEMAKE runs

n_sim  = snakemake@config[["n_sim"]]
n_tips = snakemake@config[["n_tips"]]
rates  = snakemake@config[["rates"]]
birth  = snakemake@config[["birth"]]
death  = snakemake@config[["death"]]
clade_size = snakemake@config[["clade_size"]]

################################################################################


info_df <- data.frame(
  tree_name = rep(NA, n_sim * length(n_tips)),
  n_tips = NA, birth = NA, death = NA,
  clade_size = NA, shift_node = NA, subclade = NA)

k = 1
for (i in 1:n_sim){
  for (n in n_tips){
    for (c in clade_size){
      sample <- sim_subclade(n, n * c, birth, death)
      tree <- sample[[1]]
      subtree <- sample[[2]]
      name = paste0("t", i, "_N",n, "_clade", c, "_b", birth, "_d", death)
      # adding info to log file
      info_df[k,"tree_name"] = name
      info_df[k,"n_tips"] = n
      info_df[k,"birth"] = birth
      info_df[k,"death"] = death
      info_df[k,"clade_size"] = n * c
      info_df[k,"shift_node"] = getMRCA(tree, subtree$tip.label)
      info_df[k,"subclade"] = paste(subtree$tip.label, collapse = " ")
      # bamm needs newick and bayestrait needs nexus
      write.nexus(tree, file = paste0("data/", name, ".nex"))
      write.tree(tree, file = paste0("data/", name, ".nwk"))
      color_subclade(tree, subtree$tip.label, outpath = "out/plots/", outname = paste0(name, ".pdf"))
      k = k + 1
      for (r in rates){
        trait <- simBM_2rates(tree, subtree, R1=1, R2=r, model = "BM", root = 1)
        write.table(trait, file = paste0("data/", name, "_rate", r, ".trait"),
                    quote = F, row.names = T, col.names = F)
      }
    }
  }
}

write.csv(info_df, "out/simulated_trees_info.csv", quote = F, row.names = F)