library(ape)
set.seed(1)
x = rtree(5,FALSE)
plot(x, show.tip.label=FALSE)
nodelabels(); # shows internal node numbers
tiplabels(); # shows numerical tip values (name collision w/ above?!)

x$edge
# this orders the rows in the edge matrix to follow an order used for post order traversal
x = reorder(x, "postorder")
x$edge

ntips = length(x$tip.label)
m = max(x$edge) # largest edge number
res = integer(m) # result list for each node
res[1:ntips] = 1L # the tips per definition always have one node below them at least
parent = x$edge[,1]
child = x$edge[,2]
for(i in 1:nrow(x$edge)){
  res[parent[i]] = res[parent[i]] + res[child[i]]  
}
res




for(i in 1:nrow(x$edge)){
  print(res[parent[i]] + res[child[i]])  
}
res
x = reorder(x, "cladewise")


# painting a subtree with different rate
t <- sim.bd.taxa(100, numbsim=1, lambda=1.5, mu=0.5, complete=FALSE)[[1]]
nodelabels()
t_oneR <- sim.char(phy = t, par = 1, root = 1, model = "BM")

sub_t <- subtrees(t)

l_subtree <- sapply(sub_t, function(x) {length(x$tip.label)})
# none of the subset have exactly half the tips
sum(l_subtree == 50)



sims <- sim.bd.taxa(100, numbsim=100, lambda=1.5, mu=0.5, complete=FALSE)


num_subtrees <- function(x, tips=50) {
  len = lapply(subtrees(x), function(y) {length(y$tip.label)})
  return(sum(len >= tips))
}

num_subtrees(sims)

l_subtree <- sapply(sims, num_subtrees)
# none of the subset have exactly half the tips
sum(l_subtree == 50)


len = lapply(subtrees(sims[[1]]), function(y) {length(y$tip.label)})
sum(len >= 50)

extract_clade <- function(tree, clade.size = 10) {
  # takes a phylogeny and checks whether there is a subtree with exactly clade.size tips
  # returns a random subclade if there is more than one
  subtrees <- subtrees(tree)
  size_index <- lapply(subtrees, function(tre) {length(tre$tip.label)}) == clade.size
  res <- subtrees[size_index]
  return(res)
}


extract_clade(sims[[1]], clade.size = 50)

sapply(sims, function(x) {length(extract_clade(tree = x, clade.size = 50)) > 0})

sim_subclade <- function(ntips, clade.size, birth, death) {
  subclade <- list()
  while(length(subclade) == 0) {
    tree = sim.bd.taxa(ntips, numbsim=1, lambda=birth, mu=death, complete=FALSE)[[1]]
    subclade <- extract_clade(tree = tree, clade.size = clade.size)
    i = 1
    if (length(subclade) > 1){
      i = sample(1:length(subclade), 1, replace = FALSE, prob = NULL)
    }
    return(list(tree, subclade[[i]]))
  }
}


sample <- sim_subclade(50, 10, 1.5, 0.5)
tree <- sample[[1]]
subtree <- sample[[2]]

# to paint the selected clade we need to find the number of the root node
# and the number of the MRCA of the selected clade
tree_root <- getMRCA(tree, tree$tip.label)
subclade_root <- getMRCA(tree, subtree$tip.label)
paintSubTree(tree = tree, node=21)

# http://blog.phytools.org/2012/05/painting-different-clades-with.html
tree<-paintSubTree(tree,node=tree_root,state="1")
tree<-paintSubTree(tree,node=subclade_root,state="2")
cols<-c("black","red"); names(cols)<-1:2
plotSimmap(tree,cols,pts=F,lwd=1,node.numbers=F)

# now we simulate traits since everything is multivariate normal I can simulate 
# for the whole tree and then replace the selected clade with another rate.
tree_p <- drop.clade(tree, tip = subtree$tip.label)
# the pruned location is replaced by a tip called NA. I rename this for clarity.
tree_p$tip.label[tree_p$tip.label == "NA"] = "shift"
trait_1 <- sim.char(phy = tree_p, par = 1, nsim = 1, model = "BM", root = 1000)
# the rest needs to be drawn with a different rate but with the appropriate starting point.
shift_index = names(trait_1[,,1]) == "shift"
trait_2 <- sim.char(phy = subtree, par = 2, nsim = 1, model = "BM",
                    root = trait_1[,,1][shift_index])
# remove the shift value from the result
trait <- c(trait_1[,,1][! shift_index], trait_2[,,1])
# creat phenogram of resulting traits
par(mfrow=c(1,1))
phenogram(tree,trait,colors=cols,fsize=0.8)

# lables in the subtree correspond to the node numbers in the tree
par(mfrow=c(3,1))
plot(tree)
nodelabels()
plot(subclade[[1]], show.node.label = T)
# coloring the tips
colors <- ifelse(tree$tip.label %in% subclade[[1]]$tip.label, "red", "black")
# names(colors) <- tree$tip.label
# colors <- colors[tree$tip.label]
plot(tree, tip.color = colors)






par(mfrow=c(1,1))
plot(tree, tip.color = trait[tree$tip.label])






par(mfrow=c(3,1))
plot(tree)
nodelabels()
plot(subclade[[1]], show.node.label = T)
plot(tree_pruned )
nodelabels()

trait_1 <- sim.char(phy = tree_pruned, par = 1, nsim = 1, model = "BM", root = 10)


library(ape)
data(bird.orders)
phy<-bird.orders

node<-43 #Node with taxa of interest (Columbiformes:Passeriformes)
site<-21 #Where above group is expected to be (In a clade with Struthioniformes:Anseriformes)

j<-allDesc(node) #See function definition below
grp<-extract.clade(phy, node)
rest<-drop.tip(phy, j)

tr<-bind.tree(rest, grp, where=site)

plot(tr) #Should get an error: "Error in plot.phylo(tr) : NA/NaN/Inf in foreign function call (arg 6)"

write.tree(tr) #Can see Upupiformes as sister to NODE11, among others.


