# plot_shift_bamm.R
source("scripts/functions.R")

threads <- snakemake@threads
rates <- snakemake@config[["sel_rate"]]
priors_num <- snakemake@config[["sel_prior_set"]]
priors <- paste0("pr", priors_num)
burnin <- snakemake@config[["burnin"]]
no_cores <- threads

cl <- makeCluster(no_cores)
setDefaultCluster(cl)
# we need to evaluate the source script on all the clusters
clusterEvalQ(cl = NULL, source("scripts/functions.R"))
print(no_cores)
# we need both the tree and the event data

runs <- unlist(
  strsplit(
    list.files(path = "data/", pattern = "t.+nwk"),
                        split = ".nwk"))

# TODO get this from config file 
sel_n_sim = snakemake@config[["sel_n_sim"]]

# rates <- c("rate1", "rate2", "rate5", "rate10", "rate100")
# priors <- c("pr1", "pr2")
# postfile <- paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_mcmc_out.txt")
# burnin = 0.1
expectedNumberOfShifts <- 1


sims <- paste0("^t", 1:sel_n_sim, "_.+")
selected_sims = c()
for (s in sims){
  selected_sims <- c(selected_sims,(grep(s, runs, value=TRUE)))
}
# run <- grep("^t1_.+", runs, value=TRUE)
fixed_rates <- paste0("_rate", rates)
all_mcmc_out_path <- apply(expand.grid("out/bamm/", selected_sims, fixed_rates, ".trait_", priors, "_mcmc_out.txt"), 1, paste, collapse = "", sep = "")
print(all_mcmc_out_path)
basic_tmp <- parLapply(cl = cl, all_mcmc_out_path, basic_stats_wrapper, expectedNumberOfShifts, burnin)
mcmc_stats_df <- rbindlist(basic_tmp) 
write.table(mcmc_stats_df, "out/bamm/basic_stats_mcmc_out.txt", row.names = FALSE)

c = 1
for (sim in sims){
  # subset by tree
  run <- grep(sim, selected_sims, value=TRUE)

  all_list <- parLapply(cl = cl, run, get_mcmc_prob_wrapper, rates, priors)
#  all_list <- list()
#  for (p in priors){
#    for (r in rates){
#      print(run)
#      res <- parLapply(cl=cl, run, get_mcmc_prob, r, p)
#      all_list <- c(all_list, res)
#  }}

  all_df <- rbindlist(all_list)
  
  pdf(paste0("out/plots/bamm_pre_post_sim", c, "_", priors, ".pdf"), width = 14, height = 14)
  ggplot(all_df) +
    geom_line(aes(x=N_shifts, y=priorProbs), color= "blue",fill= "blue") +
    geom_line(aes(x=N_shifts, y=postProbs), color= "red") +
    facet_grid(rate_par~N_tips) + 
    geom_vline(xintercept = 1) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8))
  dev.off()
  c = c + 1
  }


write.table("placeholder", file="bamm.report")
# 
# tree_name <- 
# 
# t10_r1_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate1", prior_par)
# t10_r2_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate2", prior_par)
# t10_r10_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate10", prior_par)
# t10_r100_runs <- parLapply(cl=cl,runs, get_mcmc_prob, "rate100", prior_par)
# 
# all_list <- c(t10_r1_runs, t10_r2_runs, t10_r10_runs, t10_r100_runs)
# all_df <- rbindlist(all_list)
# 
# 
# pdf(paste0("out/plots/bamm_pre_post_", prior_par, ".pdf"), width = 14, height = 14)
# ggplot(all_df) +
#   geom_line(aes(x=N_shifts, y=priorProbs), color= "blue",fill= "blue") +
#   geom_line(aes(x=N_shifts, y=postProbs), color= "red") +
#   facet_grid(rate_par~N_tips) + 
#   geom_vline(xintercept = 1) +
#   theme_bw() +
#   theme(strip.text.x = element_text(size = 8))
# dev.off()
# 
# 
# 
# 
# tree_name = "t1_N400_clade0.2_b1.5_d0.5"
# rate_par = "rate10"
# prior_par = "pr1"
# 
# 
# 
# 
# 
# 
# 
# 
# 
# get_mcmc_prob(tree_name, rate_par, prior_par)
# 
# t10_runs <- lapply(runs, get_mcmc_prob, rate_par, prior_par)
# 
# all_df <- rbindlist(t10_runs )
# 
# 
# 
# shift_probs <- summary(edata)
# 
# 
# postfile <- paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_mcmc_out.txt")
# computeBayesFactors(postfile , expectedNumberOfShifts=1, burnin=0.1)
# 
# 
# mcmcout <- read.csv(postfile , header=T)
# plot(mcmcout$logLik ~ mcmcout$generation)
# burnstart <- floor(0.1 * nrow(mcmcout))
# postburn <- mcmcout[burnstart:nrow(mcmcout), ]
# library(coda)
# effectiveSize(postburn$N_shifts)
# effectiveSize(postburn$logLik)
# 
# plot(postburn$logLik ~ postburn$generation)
# 
# post_probs <- table(postburn$N_shifts) / nrow(postburn)
# names(post_probs)
# post_probs['1'] / post_probs['2']
# 
# for (run in runs[2]){
#   print(run)
#   edata_file <- paste0("out/bamm/", run, "_", rate_par, ".trait_", prior_par, "_event_data.txt")
#   tree <- read.tree(paste0("data/", run, ".nwk"))
#   edata <- getEventData(tree, eventdata = edata_file
#                         ,type= "trait"   , burnin=0.99)
#   css <- credibleShiftSet(edata, expectedNumberOfShifts=1)
#   plot.credibleshiftset(css, lwd=1.7, plotmax=4)
# }



