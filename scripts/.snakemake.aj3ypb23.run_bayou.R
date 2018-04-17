
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('data/t10_N50_clade0.2_b1.5_d0.5_rate10.trait', 'config_files/pr1.bayou', "trait_f" = 'data/t10_N50_clade0.2_b1.5_d0.5_rate10.trait', "prior_f" = 'config_files/pr1.bayou'),
    output = list('t10_N50_clade0.2_b1.5_d0.5_rate10_pr1.bayou', 'out/bayou/t10_N50_clade0.2_b1.5_d0.5_rate10_pr1.summary.txt', "checkpoint" = 't10_N50_clade0.2_b1.5_d0.5_rate10_pr1.bayou', "ch_sum" = 'out/bayou/t10_N50_clade0.2_b1.5_d0.5_rate10_pr1.summary.txt'),
    params = list(),
    wildcards = list('t10_N50_clade0.2_b1.5_d0.5_rate10', '1', "sample" = 't10_N50_clade0.2_b1.5_d0.5_rate10', "n" = '1'),
    threads = 1,
    log = list('log/t10_N50_clade0.2_b1.5_d0.5_rate10_pr1.bayou.log'),
    resources = list(),
    config = list("n_sim" = 100, "n_tips" = c(10, 50, 100, 200, 400), "rates" = c(1, 2, 5, 10, 100, 1000), "clade_size" = c(0.2), "birth" = 1.5, "death" = 0.5, "sel_n_sim" = 10, "sel_n_tips" = c(10, 50, 100, 200, 400), "sel_clade" = c(0.2), "sel_rate" = c(10), "sel_prior_set" = c(1), "sel_tool" = c('bayou', 'bt', 'bamm'), "n_gen" = 10000000, "burnin" = 0.3, "n_chains" = 2, "mcmc_write_freq" = 1000, "mcmc_print_freq" = 1000000, "mcmc_event_freq" = 10),
    rule = 'run_bayou'
)
######## Original script #########
library(bayou)
################################################################################

n_gen  = snakemake@config[["n_gen"]]
burnin = snakemake@config[["burnin"]]
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
  t_ch1, save.Rdata = T, cleanup = F, file = paste0(outdir, prefix, "_1.rds"))
ch1 <- set.burnin(ch1, burnin)

ch2 <- load.bayou(
  t_ch2, save.Rdata = T, cleanup = F, file = paste0(outdir, prefix, "_2.rds"))
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
