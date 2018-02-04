library(ape)
library(geiger)
library(TreeSim)
library(phytools)


extract_clade <- function(tree, clade.size = 10) {
  # takes a phylogeny and checks whether there is a subtree with exactly clade.size tips
  # returns all subclades if there is more than one
  subtrees <- subtrees(tree)
  size_index <- lapply(subtrees, function(tre) {length(tre$tip.label)}) == clade.size
  res <- subtrees[size_index]
  return(res)
}


sim_subclade <- function(ntips, clade.size, birth, death) {
  # simulates a phylogeny with ntips that contains at least one subclade of size clade.size
  # returns list with tree and randomly chosen subclade
  subclade <- list()
  while(length(subclade) == 0) {
    tree = sim.bd.taxa(ntips, numbsim=1, lambda=birth, mu=death, complete=FALSE)[[1]]
    subclade <- extract_clade(tree = tree, clade.size = clade.size)
    i = 1
    if (length(subclade) > 1){
      i = sample(1:length(subclade), 1, replace = FALSE, prob = NULL)
    }
  }
  return(list(tree, subclade[[i]]))
}


simBM_2rates <- function(tree, subtree, R1=1, R2=10, model = "BM", root = 10) {
  # given a tree with a monophyletic subtree we simulate continuous trait data with two rates
  # returns named vector of trait values
  tree_p <- drop.clade(tree, tip = subtree$tip.label)
  # the pruned location is replaced by a tip called NA. I rename this for clarity.
  tree_p$tip.label[tree_p$tip.label == "NA"] = "shift"
  trait_1 <- sim.char(phy = tree_p, par = R1, nsim = 1, model = "BM", root = root)
  # the rest needs to be drawn with a different rate but with the appropriate starting point.
  shift_index = names(trait_1[,,1]) == "shift"
  trait_2 <- sim.char(phy = subtree, par = R2, nsim = 1, model = "BM",
                      root = trait_1[,,1][shift_index])
  # remove the shift value from the result
  trait <- c(trait_1[,,1][! shift_index], trait_2[,,1])
  return(trait)
}