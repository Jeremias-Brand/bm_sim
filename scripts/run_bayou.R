library(bayou)
################################################################################

n_gen  = snakemake@config[["n_gen"]]
burnin = snakemake@config[["burnin"]]
bayou_no_stones = snakemake@config[["bayou_no_stones"]]
bayou_iter_stones = snakemake@config[["bayou_iter_stones"]]
# tree_f = snakemake@input[["tree_f"]]
trait_f = snakemake@input[["trait_f"]]
prior_f = snakemake@input[["prior_f"]]
tree_f = paste0(strsplit(trait_f, "_rate")[[1]][1], ".nwk")
checkpoint = snakemake@output[["checkpoint"]]

# READ DATA
# tree_f = "data/t10_N100_clade20_b1.5_d0.5.nwk"
# trait_f = "data/t10_N100_clade20_b1.5_d0.5_rate10.trait"
source(prior_f)
prefix = strsplit( strsplit(trait_f, ".trait")[[1]][1], "/")[[1]][2]
pr = strsplit( strsplit(prior_f, ".bayou")[[1]][1], "/")[[1]][2]

outdir <- "./out/bayou/"
tempdir=outdir
# READ PARAMETER
tree <- read.newick(tree_f)
trait_tmp <- read.table(trait_f)
trait <- trait_tmp$V2
names(trait) <- trait_tmp$V1; rm(trait_tmp)

if (!name.check(tree, trait) == "OK") {
  stop(paste("Names in", tree_f, "do not match names in", trait_f))
}

prior.fixed <- eval(prior_q)


chain_cmd <- quote(
  bayou.mcmc(
    tree = tree, dat = trait, SE = 0, model = "OU", prior = prior.fixed, startpar=startpar,
    ngen = n_gen, samp = 10, chunk = 100, control = NULL, tuning = NULL,
    new.dir = outdir, plot.freq = NULL, ticker.freq = 10000,
    outname = prefix
  )
)

ss1_cmd <- quote(
    steppingstone(
        Bk=seq(0,1, length.out=bayou_no_stones), chain = ch1, tree = tree, dat = traits,
        SE=SE, prior=e_prior, startpar = startpar, burnin = burnin,
        new.dir = tempdir, ngen=, parallel = FALSE
        )
    )

ss2_cmd <- quote(
    steppingstone(
        Bk=seq(0,1, length.out=bayou_no_stones), chain = ch2, tree = tree, dat = traits,
        SE=SE, prior=e_prior, startpar = startpar, burnin = burnin,
        new.dir = tempdir, ngen=ss_ngen, parallel = FALSE
        )
    )
  


t_ch1 <- eval(chain_cmd)
t_ch2 <- eval(chain_cmd)
# loading the chains
ch1 <- load.bayou(
  t_ch1, cleanup = F, file = paste0(outdir, prefix, "_", pr, "_1.rds"))
ch1 <- set.burnin(ch1, burnin)

ch2 <- load.bayou(
  t_ch2,cleanup = F, file = paste0(outdir, prefix, "_", pr, "_2.rds"))
ch2 <- set.burnin(ch2, burnin)


print(paste("Stepping Stone chain 1", msg))
ss1 <- eval(ss1_cmd)
print(paste("Stepping Stone chain 2", msg))
ss2 <- eval(ss2_cmd)
ss1 <- set.burnin(chain = ss1, burnin = burnin)
lh1 <- ss1$lnr
ss2 <- set.burnin(chain = ss2, burnin = burnin)
lh2 <- ss2$lnr

sink(paste0(outdir, prefix, "_", pr, ".summary.txt"))
print("# Chain 1")
summary(ch1)
print("# Chain 2")
summary(ch2)  
print("# Stepping stone chain 1")
print(ss1)
print("# Stepping stone chain 2")
print(ss2)
sink()


G_lnl <- gelman.R("lnL", chain1=ch1, chain2=ch2,
                  plot=F, type="n", ylim=c(0.9, 2))

G_sig2 <- gelman.R("sig2", chain1=ch1, chain2=ch2,
         plot=F, type="n", ylim=c(0.9, 2))


# Output statistics about convergence
write.table(G_lnl, paste0(outdir, prefix, "_", pr, ".gelman.lnl.txt"), quote = F, row.names = F)
write.table(G_sig2, paste0(outdir, prefix, "_", pr, ".gelman.sig2.txt"), quote = F, row.names = F)

# plotting the parameters
pdf(paste0(outdir, prefix, "_", pr, ".ss.pdf"))
  eval(prior_q)
  par(mfrow=c(1,1))
    gelman.R("lnL", chain1=ch1, chain2=ch2,
                        plot=TRUE, type="n", ylim=c(0.9, 2))
    gelman.R("sig2", chain1=ch1, chain2=ch2,
                        plot=TRUE, type="n", ylim=c(0.9, 2))
    plot(ss1)
    plot(ss2)
dev.off()
write.csv(x = "", file = checkpoint)
