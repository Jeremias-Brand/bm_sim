# plot_shift_bamm.R
library(ape)
library(geiger)
library(phytools)
library(BAMMtools)
library(data.table)

# we need both the tree and the event data

runs <- unlist(
  strsplit(
    list.files(path = "data/", pattern = "t1_.+nwk"),
                        split = ".nwk"))



tree_name = "t1_N400_clade0.2_b1.5_d0.5"
rate_par = "rate10"
prior_par = "pr1"

get_mcmc_prob <- function(tree_name, rate_par, prior_par) {
  library(ape)
  library(geiger)
  library(phytools)
  library(BAMMtools)
  library(data.table)
  N_tips <- strsplit(tree_name, "_")[[1]][2]
  tree <- read.tree(paste0("data/", tree_name, ".nwk"))
  edata <- getEventData(tree, eventdata = paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_event_data.txt")
                        ,type= "trait"   , burnin=0.1)
  css <- credibleShiftSet(edata, expectedNumberOfShifts=1)
  
  # pdf(paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_credibleshiftset.pdf"), width=14, height = 14)
  # plot.credibleshiftset(css, lwd=1.7, plotmax=9)
  # dev.off()
  postfile <- paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_mcmc_out.txt")
  pdf(paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_prior_post.pdf"), width=14, height = 14)
  pripost_tbl <- plotPrior(postfile , expectedNumberOfShifts=1)
  dev.off()
  
  pripost_df <- as.data.frame(pripost_tbl)
  pripost_df$N_tips <- N_tips
  pripost_df$rate_par <- rate_par
  pripost_df$prior_par <- prior_par
  pripost_df$tree_name <- tree_name
  
  return(pripost_df)
  }

library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)


t10_r1_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate1", prior_par)
t10_r2_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate2", prior_par)
t10_r10_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate10", prior_par)
t10_r100_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate100", prior_par)

all_list <- c(t10_r1_runs, t10_r2_runs, t10_r10_runs, t10_r100_runs)
all_df <- rbindlist(all_list)

library(ggplot2)
pdf("out/plots/bamm_Pre_post.pdf", width = 14, height = 14)
ggplot(all_df) +
  geom_line(aes(x=N_shifts, y=priorProbs), color= "blue",fill= "blue") +
  geom_line(aes(x=N_shifts, y=postProbs), color= "red") +
  facet_grid(rate_par~N_tips) + 
  geom_vline(xintercept = 1) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))
dev.off()


get_mcmc_prob(tree_name, rate_par, prior_par)

t10_runs <- lapply(runs, get_mcmc_prob, rate_par, prior_par)

all_df <- rbindlist(t10_runs )



shift_probs <- summary(edata)


postfile <- paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_mcmc_out.txt")
computeBayesFactors(postfile , expectedNumberOfShifts=1, burnin=0.1)


mcmcout <- read.csv(postfile , header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

plot(postburn$logLik ~ postburn$generation)

post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
post_probs['1'] / post_probs['2']

for (run in runs[2]){
  print(run)
  tree <- read.tree(paste0("data/", run, ".nwk"))
  edata <- getEventData(tree, eventdata = paste0("out/bamm/", run, "_", rate_par, ".trait_", prior_par, "_event_data.txt")
                        ,type= "trait"   , burnin=0.99)
  css <- credibleShiftSet(edata, expectedNumberOfShifts=1)
  plot.credibleshiftset(css, lwd=1.7, plotmax=4)
}



