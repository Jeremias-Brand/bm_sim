library(ape)
library(geiger)
library(phytools)
source("scripts/functions.R")
# simulation viz
# this script produces plots illustrating what output the simulation produces
n_tips <- c(10, 50, 100, 200, 400)
rates <- c(1,2,5,10,100,1000)
birth = 1.5
death = 0.5
nsim=1


for (n in n_tips){
  
  sample <- sim_subclade(n, n/5, birth, death)
  tree <- sample[[1]]
  subtree <- sample[[2]]
  pdf(paste0("fig/tree_sim_N",n, "_clade", n/5, "_b", birth, "_d", death ,".pdf"))
  
  for (r in rates){
    trait <- simBM_2rates(tree, subtree, R1=1, R2=r, model = "BM", root = 1)
    # paint the selcted clade red
    # to paint the selected clade we need to find the number of the root node
    # and the number of the MRCA of the selected clade
    tree_root <- getMRCA(tree, tree$tip.label)
    subclade_root <- getMRCA(tree, subtree$tip.label)
    # http://blog.phytools.org/2012/05/painting-different-clades-with.html
    tree<-paintSubTree(tree,node=tree_root,state="1")
    tree<-paintSubTree(tree,node=subclade_root,state="2")
    cols<-c("black","red"); names(cols)<-1:2
    # hack to not have labels in plotSimmap
    tmp <- tree
    tmp$tip.label <- rep("",length(tree$tip.label))
    
    par(mfrow=c(2,1), oma = c(0, 0, 2, 0))
    plotSimmap(tmp,cols,pts=F,lwd=1,node.numbers=F)
    # creat phenogram of resulting traits, not tip labels
    phenogram(tree,trait,colors=cols,fsize=0, spread.labels=F)
    mtext(paste("R2/R1 =", r), outer = TRUE, cex = 1.5)
  }
  dev.off()
}




