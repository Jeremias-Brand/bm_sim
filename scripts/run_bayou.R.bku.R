library(bayou)
################################################################################

n_gen  = snakemake@config[["n_gen"]]
burnin = snakemake@config[["burnin"]]
tree_f = snakemake@input[["tree_f"]]
trait_f = snakemake@input[["trait_f"]]
prior_f = snakemake@input[["prior_f"]]
tree_f = paste0(strsplit(trait_f, "_rate")[[1]][1], ".nwk")

# READ DATA
# tree_f = "data/t10_N100_clade20_b1.5_d0.5.nwk"
# trait_f = "data/t10_N100_clade20_b1.5_d0.5_rate10.trait"
source(prior_f)
prefix = strsplit( strsplit(dd, ".trait")[[1]][1], "/")[[1]][2]
outdir <- "./out/bayou/"

# READ PARAMETER
runid="1"
burnin=0.3

tree <- read.newick(tree_f)
trait_tmp <- read.table(trait_f)
trait <- trait_tmp$V2
names(trait) <- trait_tmp$V1; rm(trait_tmp)

if (!name.check(tree, trait) == "OK") {
  stop(paste("Names in", tree_f, "do not match names in", trait_f))
}

ngen = 10000
fixed.pars <- list(alpha= 1e-10)
prior.fixed <- quote(
              make.prior(tree,
                          dists=list(
                            dalpha="fixed", dsig2="dhalfcauchy",
                            dsb="dsb", dk="cdpois", dtheta="dnorm"
                          ), 
                          param=list(
                            dsig2=list(scale=1), dk=list(lambda=15, kmax=200),
                            dsb=list(bmax=1,prob=1), dtheta=list(mean=mean(trait), sd=2)
                          ),
                          fixed=fixed.pars,
                         plot.prior = F, model = "OU")
)

startpar <- list(alpha=1e-10, sig2=1, k=1, ntheta=2, theta=c(1,2), 
                 sb=1, t2=2, loc=0)


prior.fixed <- eval(prior.fixed)


chain_cmd <- quote(
  bayou.mcmc(
    tree = tree, dat = trait, SE = 0, model = "OU", prior = prior.fixed, startpar=startpar,
    ngen = ngen, samp = 10, chunk = 100, control = NULL, tuning = NULL,
    new.dir = outdir, plot.freq = NULL, ticker.freq = 10000,
    outname = prefix
  )
)



t_ch1 <- eval(chain_cmd)
t_ch2 <- eval(chain_cmd)
# loading the chains and cleanup
ch1 <- load.bayou(
  t_ch1, save.Rdata = T, cleanup = F, file = paste0(outdir,runid,"_1.rds"))
ch1 <- set.burnin(ch1, burnin)
ch2 <- load.bayou(
  t_ch2, save.Rdata = T, cleanup = F, file = paste0(outdir,runid,"_2.rds"))
ch2 <- set.burnin(ch2, burnin)



gelman.R("lnL", chain1=ch1, chain2=ch2,
         plot=TRUE, type="n", ylim=c(0.9, 2))
gelman.R("alpha", chain1=ch1, chain2=ch2,
         plot=TRUE, type="n", ylim=c(0.9, 2))
gelman.R("sig2", chain1=ch1, chain2=ch2,
         plot=TRUE, type="n", ylim=c(0.9, 2))
