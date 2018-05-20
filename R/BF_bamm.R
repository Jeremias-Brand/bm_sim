library(coda)
library(BAMMtools)
library(ggplot2)
# test to calcultate 95% CI for bamm and bt
# https://www.r-bloggers.com/probable-points-and-credible-intervals-part-1/
# there are multiple ways we can calculate a credible interval,
# most useful is the higest Density Interval because it is the narrowest for a given
# level of plausibility. It incluedes the most plausible values.

# we need a credible interval, posterior credibility

mcmc_out <- read.table("out/bamm/t1_N100_clade0.2_b1.5_d0.5_rate100.trait_pr2_mcmc_out.txt", sep = ",", header = TRUE)
event_dat <- read.table("out/bamm/t1_N100_clade0.2_b1.5_d0.5_rate100.trait_pr2_event_data.txt", sep = ",", header = TRUE)



expectedNumberOfShifts = 1

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
  
  stats_df <- data.frame(k_mode = k_mode, k_mean = k_mean, k_median = k_median, k_sd = k_sd,
                         k_HPD95_low = HPD_0.95[1], k_HPD95_high = HPD_0.95[2],
                         k_quantile95_low = q_int[1], k_quantile95_high = q_int[2],
                         BF_modeOver0 = BF_modeOver0, BF_1Over0 = BF_1Over0, BF_modeOver1 = BF_modeOver1)
  return(stats_df)
}

basic_stats_bamm_mcmc("out/bamm/t1_N100_clade0.2_b1.5_d0.5_rate100.trait_pr2_mcmc_out.txt", 1, 0.1)



ggplot(mcmc_out) +
  geom_density(aes(x=N_shifts), fill = "magenta") +
  geom_point(aes(x=mode,y=-0.02), color="blue") +
  geom_segment(aes(x=HPD_0.95[1], xend=HPD_0.95[2],y=-0.02, yend=-0.02), color="blue") +
  geom_point(aes(x=median,y=-0.03), color="green") +
  geom_segment(aes(x=q_int[1], xend=q_int[2],y=-0.03, yend=-0.03), color="green") +
  geom_point(aes(x=mean,y=-0.04), color="red") +
  geom_segment(aes(x=mean - 1.96 * SD, xend=mean + 1.96 * SD,y=-0.04, yend=-0.04), color="red") +
  theme_bw()



bfmat <- computeBayesFactors(mcmc_out, expectedNumberOfShifts=1, burnin=0.1)




pripost_tbl <- plotPrior(mcmc_out , expectedNumberOfShifts=1)

length(mcmc_out$N_shifts)
sum(mcmc_out$N_shifts == 1)

# the posterior probability of rate shifts is simply the frequencig this configuration
# appears in the posterior


expectedNumberOfShifts = 1
# I copy what the BAMMtools formula does to calculate bayes factors
post_probs <- table(mcmc_out$N_shifts) / nrow(mcmc_out)
configurations <- as.numeric(names(post_probs)) # this lists all the shifts that occured
prior_prob <- (1/(1 + expectedNumberOfShifts)) # that is the prior prob for 0 shifts
pp
prior <- dgeom(configurations, prob = prior_prob) # under the prior higher shift configurations just follow a poison 
# distribution and therefore we can calculate the probability to fail x times to get the prior.
names(prior) <- configurations # we need to name it gaing so we remember what shifts hat what prior probability.

# create a results matrix to be filled in.
mm <- matrix(NA, nrow = length(prior), ncol = length(prior))
rownames(mm) <- names(prior)
colnames(mm) <- names(prior)
post <- data.frame(N_shifts = as.numeric(names(post_probs)), prob = as.numeric(post_probs))

prior_odds <- prior[3]/prior[2] # odd for 1 vs 2 shift
post_odds <- post$prob[post$N_shifts == 2]/post$prob[post$N_shifts == 1] # this is the odds in the posterior. THis is much smaller

post_odds * (1/prior_odds)





burnstart <- floor(0.1 * nrow(mcmc_out))
postburn <- mcmc_out[burnstart:nrow(mcmc_out), ]

BF_pairwise(postburn, expectedNumberOfShifts = 1, m0 = 1, m1 = 8)

bfmat

for (i in 1:length(prior)) {
  mi <- ux[i]
  for (j in 1:length(prior)) {
    mj <- ux[j]
    prior_odds <- prior[i]/prior[j]
    post_odds <- post$prob[post$N_shifts == mi]/post$prob[post$N_shifts == 
                                                            mj]
    mm[i, j] <- post_odds * (1/prior_odds)
  }
}



sum(mcmc_out$N_shifts == 0)/length(mcmc_out$N_shifts)
sum(mcmc_out$N_shifts == 1)/length(mcmc_out$N_shifts)
sum(mcmc_out$N_shifts == 2)/length(mcmc_out$N_shifts)
sum(mcmc_out$N_shifts == 3)/length(mcmc_out$N_shifts)
sum(mcmc_out$N_shifts == 4)/length(mcmc_out$N_shifts)


