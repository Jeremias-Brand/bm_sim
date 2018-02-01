library(ape)
library(phytools)
library(geiger)
library(TreeSim)

n_tips <- c(10, 50, 100, 200, 400)
birth = 1.5
death = 0.5
nsim=2

sims <- lapply(n_tips, FUN=sim.bd.taxa, numbsim=nsim, lambda=birth, mu=death, complete=FALSE)

for (s_class in sims){
  for (tree in s_class){
    print(plot(tree))
    print(sim.char(tree, 1, 1, root = 10))
  }  
}

tree <- sims[[1]][[1]]
trait <- sim.char(tree, 1, 1, root = 10)
par(mar=c(2,0,0,1))
plot(tree)

write.nexus(tree, file = "t1.nex")
write.tree(tree, file = "t1.nwk")
write.table(x = trait, file = "t1_trait.tsv",
            quote = F, row.names = T, col.names = F)

tiplabels(pch=21, bg = trait, adj = 0.5)

make.simmap(tree, trait)
lapply(dd, plot)

trees <- sim.bd.taxa(13, numbsim=1, lambda=1.5, mu=0.5, complete=FALSE)

# getting decendants of node
t <- sim.bd.taxa(13, numbsim=1, lambda=1.5, mu=0.5, complete=FALSE)[[1]]
nodelabels()
getDescendants(t, 14, curr=NULL)

par(mfrow = c(2,3))
tre <- sim.bd.taxa(n = 50, numbsim = 9, lambda = 1.5, mu = 0.5, frac = 1, complete = F)
lapply(X = tre, plot)
s_c <- sim.char(tre[[1]], 1, 1)
plot(s_c)

s_c <- sim.char(tre[[1]], 100000000, 1)
hist(s_c)

tre <- sim.bd.taxa(n = 50, numbsim = 6, lambda = 1, mu = 1, frac = 1, complete = F)
lapply(X = tre, plot)

geo <- get(data(geospiza))
s <- ratematrix(geo$phy, geo$dat)
csims <- sim.char(geo$phy, s, 100)


# painting a subtree with different rate
t <- sim.bd.taxa(100, numbsim=1, lambda=1.5, mu=0.5, complete=FALSE)[[1]]
nodelabels()
t_oneR <- sim.char(phy = t, par = 1, root = 1, model = "BM")

sub_t <- subtrees(t)

l_subtree <- sapply(sub_t, function(x) {length(x$tip.label)})
sum(l_subtree > 50)


plot.phylo(tre[[1]])
male.gen.sim.slow <- lapply(1:100, function(i) sapply(male.gen, function(x) sim.char(tree[[i]], par=0.5, model="BM", root=x) ) )
male.gen.sim.fast <- lapply(1:100, function(i) sapply(male.gen, function(x) sim.char(tree[[i]], par=1, model="BM", root=x) ) )
female.gen.sim <- lapply(1:100, function(i) sapply(female.gen, function(x) sim.char(tree[[i]], par=0.5, model="BM", root=x) ) )