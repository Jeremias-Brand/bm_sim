#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
      stop("USAGE: run_bayou_test.R --vanilla trait_f prior_f tree_f", call.=FALSE)
} else if (length(args)==1) {
      # default output file
      args[2] = "out.txt"
}
library(bayou)
################################################################################

n_gen  = 10000
burnin = 0.1
# tree_f = snakemake@input[["tree_f"]]
trait_f = args[1]
prior_f = args[2]
tree_f = paste0(strsplit(trait_f, "_rate")[[1]][1], ".nwk")
checkpoint = "bayou_test_dummy"

# READ DATA
# tree_f = "data/t10_N100_clade20_b1.5_d0.5.nwk"
# trait_f = "data/t10_N100_clade20_b1.5_d0.5_rate10.trait"
source(prior_f)
prefix = strsplit( strsplit(trait_f, ".trait")[[1]][1], "/")[[1]][2]
pr = strsplit( strsplit(prior_f, ".bayou")[[1]][1], "/")[[1]][2]

outdir <- "./out/testdir/"

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


t_ch1 <- eval(chain_cmd)
t_ch2 <- eval(chain_cmd)
# loading the chains
ch1 <- load.bayou(
  t_ch1, cleanup = F, file = paste0(outdir, prefix, "_1.rds"))
ch1 <- set.burnin(ch1, burnin)

ch2 <- load.bayou(
  t_ch2,cleanup = F, file = paste0(outdir, prefix, "_2.rds"))
ch2 <- set.burnin(ch2, burnin)

sink(paste0(outdir, prefix, "_", pr, ".summary.txt"))
print("# Chain 1")
summary(ch1)
print("# Chain 2")
summary(ch2)  
sink()


G_lnl <- gelman.R("lnL", chain1=ch1, chain2=ch2,
                  plot=F, type="n", ylim=c(0.9, 2))

G_sig2 <- gelman.R("sig2", chain1=ch1, chain2=ch2,
         plot=F, type="n", ylim=c(0.9, 2))

write.table(G_lnl, paste0(outdir, prefix, ".gelman.lnl.txt"), quote = F, row.names = F)
write.table(G_sig2, paste0(outdir, prefix, ".gelman.sig2.txt"), quote = F, row.names = F)
write.csv(x = "", file = checkpoint)
