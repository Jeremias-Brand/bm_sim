library(ape)
library(geiger)
library(TreeSim)
library(phytools)
library(coda)
library(phytools)
library(BAMMtools)
library(parallel)
library(ggplot2)
library(data.table)


extract_clade <- function(tree, clade.size = 10) {
  # takes a phylogeny and checks whether there is a subtree with exactly clade.size tips
  # returns all subclades if there is more than one
  subtrees <- subtrees(tree)
  size_index <- lapply(subtrees, function(tre) {length(tre$tip.label)}) == clade.size
  res <- subtrees[size_index]
  return(res)
}


sim_subclade <- function(ntips, clade.size, birth, death) {
  require(TreeSim)
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


color_subclade <- function(tree, subclade, outpath, outname) {
  # given a tree and a monophyletic subclade this produces a plot with the sublaced colored red
  require(ape)
  if(!is.monophyletic(tree, subclade)) stop("the sublcade you specify is not monophyletic")
  subclade_edges <-  which.edge(tree, subclade)
  edge_col <- rep("black", dim(tree$edge)[1]) #defining a color vector; look up rep() and dim() functions!
  edge_col[subclade_edges] <- "red" #specifying color of branches
  nTips = length(tree$tip.label)
  if(nTips < 20){
    pdf(file = paste0(outpath, outname), width = 7, height = 7)
  } else {
    # for trees with many tips we need to scale the pdf
    pdf(file = paste0(outpath, outname),
        width = 2 * log10(length(tree$tip.label)),
        height = 0.2 * length(tree$tip.label))
  }
  plot(tree, edge.color=edge_col)
  dev.off()
}


estimate_mode <- function(s) {
  d <- density(s)
  d$x[which.max(d$y)]
}


basic_stats_bamm_mcmc <- function(mcmc_file, expectedNumberOfShifts, burnin) {
  mcmc_out <- read.table(mcmc_file, sep = ",", header = TRUE)
  burnstart <- floor(burnin * nrow(mcmc_out))
  postburn <- mcmc_out[burnstart:nrow(mcmc_out), ]
  ch <- mcmc(mcmc_out)
  # Collect some basic stats about the run
  HPD_int <- HPDinterval(ch, 0.95)
  HPD_0.95 <- HPD_int["N_shifts",]
  k_mode <- estimate_mode(mcmc_out$N_shifts)
  k_median <- median(mcmc_out$N_shifts)
  k_mean <- mean(mcmc_out$N_shifts)
  k_sd <- sd(mcmc_out$N_shifts)
  q_int <- quantile(mcmc_out$N_shifts, c(0.025, 0.975))
  
  BF_modeOver0 <- BF_pairwise(postburn, expectedNumberOfShifts = expectedNumberOfShifts, m0 = 0, m1 = k_mode)
  BF_1Over0 <- BF_pairwise(postburn, expectedNumberOfShifts = expectedNumberOfShifts, m0 = 0, m1 = 1)
  BF_modeOver1 <- BF_pairwise(postburn, expectedNumberOfShifts = expectedNumberOfShifts, m0 = 1, m1 = k_mode)
  
  stats_df <- data.frame(expectedNumberOfShifts = expectedNumberOfShifts, samples = nrow(postburn), k_mode = k_mode, k_mean = k_mean, k_median = k_median, k_sd = k_sd,
                         k_HPD95_low = HPD_0.95[1], k_HPD95_high = HPD_0.95[2],
                         k_quantile95_low = q_int[1], k_quantile95_high = q_int[2],
                         BF_modeOver0 = BF_modeOver0, BF_1Over0 = BF_1Over0, BF_modeOver1 = BF_modeOver1)
  return(stats_df)
}



basic_stats_wrapper <- function(mcmc_file, expectedNumberOfShifts, burnin) {
  # here we need a label based on the file name for some applications
  # we get this by parsing the regular file name this only works if the name follows
  # this pattern:
  # out/bamm/t1_N10_clade0.2_b1.5_d0.5_rate1.trait_pr1_mcmc_out.txt
  f_name <- strsplit(mcmc_file, "/")[[1]][3]
  tree <- strsplit(f_name, "_rate1")[[1]][1]
  rate <- strsplit(f_name, "_")[[1]][6]
  prior <- strsplit(f_name, "_")[[1]][7]
  N_tips <- strsplit(f_name, "_")[[1]][2]
  res <- basic_stats_bamm_mcmc(mcmc_file, expectedNumberOfShifts, burnin)
  res$tree <- tree
  res$rate <- rate
  res$prior <- prior
  return(res)
  }


get_mcmc_prob_wrapper <- function(run, rate_pars, prior_pars){
# this version allows more efficient parallelization
  temp_list = list()
  for (p in prior_pars) {
    for (r in rate_pars){
      res <- lapply(run, get_mcmc_prob, r, p)
      temp_list = c(temp_list, res)
    }
  }
  return(temp_list)
}


get_mcmc_prob <- function(tree_name, rate_par, prior_par) {
  library(ape)
  library(geiger)
  library(phytools)
  library(BAMMtools)
  library(data.table)
  N_tips <- strsplit(tree_name, "_")[[1]][2]
  tree <- read.tree(paste0("data/", tree_name, ".nwk"))
  edata <- getEventData(tree, eventdata = paste0("out/bamm/", tree_name, "_rate", rate_par, ".trait_", prior_par, "_event_data.txt")
                        ,type= "trait"   , burnin=0.1)
  css <- credibleShiftSet(edata, expectedNumberOfShifts=1)
  
  # pdf(paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_credibleshiftset.pdf"), width=14, height = 14)
  # plot.credibleshiftset(css, lwd=1.7, plotmax=9)
  # dev.off()
  postfile <- paste0("out/bamm/", tree_name, "_rate", rate_par, ".trait_", prior_par, "_mcmc_out.txt")
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

BF_pairwise <- function(mcmc_out, expectedNumberOfShifts = 1, m0 = 0, m1 = 1) {
  # The input must be an integer
  m1 = round(m1)
  # I copy what the BAMMtools formula does to calculate bayes factors
  post_probs <- table(mcmc_out$N_shifts) / nrow(mcmc_out)
  post <- data.frame(N_shifts = as.numeric(names(post_probs)), prob = as.numeric(post_probs))
  configurations <- as.numeric(names(post_probs)) # this lists all the shifts that occured
  prior_prob <- (1/(1 + expectedNumberOfShifts)) # that is the prior prob for 0 shifts
  prior <- dgeom(configurations, prob = prior_prob) # under the prior higher shift configurations just follow a poison 
  # distribution and therefore we can calculate the probability to fail x times to get the prior.
  names(prior) <- configurations # we need to name it gaing so we remember what shifts hat what prior probability.
  
  prior_odds <- prior[m1 + 1]/prior[m0 + 1] # odd for 1 vs 2 shift
    
  if (sum(post$N_shifts == m1) == 0){
    return(NA)
  }
  # if we don't have 0 sampled they recommend we take just 1/(number of samples)
  # as the probability. 
  if (m0 == 0 & sum(post$N_shifts == m0) == 0){
    post_odds <- post$prob[post$N_shifts == m1]/(1/nrow(mcmc_out))
  } else if (sum(post$N_shifts == m0) == 0){
    return(NA)
  } else {
    post_odds <- post$prob[post$N_shifts == m1]/post$prob[post$N_shifts == m0]
  }
  BF <- post_odds * (1/prior_odds)
  names(BF) <- paste0(m1, " over ", m0)
  return(BF)
}
