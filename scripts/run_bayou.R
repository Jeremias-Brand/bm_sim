library(bayou)
################################################################################

# READ DATA
tree_f = "data/t10_N100_clade20_b1.5_d0.5.nwk"
trait_f = "data/t10_N100_clade20_b1.5_d0.5_rate10.trait"


# READ PARAMETER


tree <- read.newick(tree_f)
trait_tmp <- read.table(trait_f)
trait <- trait_tmp$V2
names(trait) <- trait_tmp$V1; rm(trait_tmp)

if (!name.check(tree, trait) == "OK") {
  stop(paste("Names in", tree_f, "do not match names in", trait_f))
}

fixed.pars <- list(alpha= 1e-10)
prior.fixed <- make.prior(tree,
                          dists=list(
                            dalpha="fixed", dsig2="dhalfcauchy",
                            dsb="dsb", dk="cdpois", dtheta="dnorm"
                          ), 
                          param=list(
                            dsig2=list(scale=1), dk=list(lambda=15, kmax=200),
                            dsb=list(bmax=1,prob=1), dtheta=list(mean=mean(trait), sd=2)
                          ),
                          fixed=fixed.pars)

startpar <- list(alpha=1e-10, sig2=1, k=1, ntheta=2, theta=c(1,2), 
                 sb=1, t2=2, loc=0)

bayou.mcmc(
  tree = tree, dat = trait, SE = 0.05, model = "OU", prior = prior.fixed, startpar=startpar,
  ngen = 10000, new.dir = ".", plot.freq = NULL, ticker.freq = 1000
  )












# startpar <- list(sig2=1, k=13, ntheta=3, theta=c(0.7,0.4,0.3), 
#                 lambda=15, kmax=200, bmax=1,prob=1)

startpar <- list(alpha=0, sig2=1, k=3)


bayou.mcmc(
  tree = tree, dat = trait, SE = 0.05, model = "OU", prior = prior.fixed, startpar=fixed.pars,
  ngen = ngen, new.dir = ".", plot.freq = NULL, ticker.freq = 1000)




startpar <- list(alpha=1, sig2=1, k=3, ntheta=3, theta=c(0.7,0.4,0.3), 
                 sb=c(12,23,31), t2=c(2,2,3))

bayou.mcmc(
  tree = tree, dat = trait, SE = 0.05, model = "OU", prior = prior.fixed, startpar=startpar,
  ngen = ngen, new.dir = ".", plot.freq = NULL, ticker.freq = 1000
  )














fixed.hypothesisK3 <- list(k=3)
prior.fixedK3<- make.prior(tree, dists=list(
    dalpha="dhalfcauchy", dsig2="dhalfcauchy",
    dsb="dsb", dk="fixed", dtheta="dnorm", dloc="dloc"
  ),
  param=list(
    dalpha=list(scale=1), dsig2=list(scale=1),
    dtheta=list(mean=mean(trait), sd=2)),
  fixed=fixed.hypothesisK3
  )


startparK3 <- list(
  alpha=0, sig2=1, ntheta=4, k=3, t2 = c(2,3,4),
  theta=c(0.7, 0.3, 0.5, 1), sb = c(38, 102, 98), loc = c(0,0,0)
)


chain_cmd <- 
  
  debug( 
    bayou.mcmc(
  tree = tree, dat = trait, SE = 0.01, model = "OU", prior = fixed.hypothesisK3,
  startpar = startparK3,
  ngen = 100000, new.dir = tempdir, plot.freq = NULL, ticker.freq = 1000)
  )





# fixed parameters need both a prior as well as start parameters
shift_loc <- c(132, 6, 8, 16, 37, 45, 59, 75, 79, 87, 100, 110, 89)
fixed.pars <- list(sb = shift_loc, loc = rep(0,13))
## Setting up fixed priors
prior.fixed13 <- make.prior(tree, dists=list(
    dalpha="dhalfcauchy",
    dsig2="dhalfcauchy",
    dsb="fixed",
    dk="fixed",
    dtheta="dnorm",
    dloc="dloc"),
    param=list(dalpha=list(scale=1),
               dsig2=list(scale=1),
               dtheta=list(mean=mean(trait), sd=2)),
    fixed=fixed.hypothesis13)


# For fixed hypothesis we need to define startparameters
startpar13 <- list(alpha=1, sig2=1, k=13, ntheta=3, theta=c(0.7,0.4,0.3))
startpar13$sb <- fixed.hypothesis13$sb
startpar13$t2 <- c(2, rep(startpar13$ntheta, 12))
startpar13$loc <- fixed.hypothesis13$loc





prior <- make.prior(tree, dists=list(
    dalpha="dhalfcauchy", dsig2="dhalfcauchy",
    dsb="dsb", dk="cdpois", dtheta="dnorm"
  ),
  param=list(
    dalpha=list(scale=1), dsig2=list(scale=1),
    dk=list(lambda=15, kmax=200), dsb=list(bmax=1,prob=1),
    dtheta=list(mean=mean(trait), sd=2))
  )



run_bayou_2chains_no_chain_summary <- function(
  tree, traits, SE, model="OU", prior, startpar=NULL, tempdir=getwd(), outdir, runid,
  ngen=100000, ss_ngen=50000, ticker.freq=10000, burnin=0.3) {
  # this function runs two bayou mcmc chains and saves a session file
  # loads the chains as runid_ch1 and runid_ch2
  # prior must be quoted: quote(make.prior(...))
  
  # WHEN RUNNING A FIXED MODEL WE NEED TO SET START PARAMETERS!
  require(bayou)
  e_prior <- eval(prior)
  if (is.null(startpar)) {
    msg <- "NO STARPAR FILE"
    chain_cmd <- quote(bayou.mcmc(
      tree = tree, dat = traits, SE = SE, model = model, prior = e_prior,
      ngen = ngen, new.dir = tempdir, plot.freq = NULL, ticker.freq = ticker.freq))
    ss1_cmd <- quote(steppingstone(
      Bk=seq(0,1,length.out=5), chain = ch1, tree = tree, dat = traits,
      SE=SE, prior=e_prior, new.dir = tempdir, ngen=ss_ngen, parallel = FALSE))
    ss2_cmd <- quote(steppingstone(
      Bk=seq(0,1,length.out=5), chain = ch2, tree = tree, dat = traits,
      SE=SE, prior=e_prior, new.dir = tempdir, ngen=ss_ngen, parallel = FALSE))
  } else {
    msg <- "WITH STARPAR FILE"
    chain_cmd <- quote(bayou.mcmc(
      tree = tree, dat = traits, SE = SE, model = model, prior = e_prior,
      startpar = startpar,
      ngen = ngen, new.dir = tempdir, plot.freq = NULL, ticker.freq = ticker.freq))
    ss1_cmd <- quote(steppingstone(
      Bk=seq(0,1, length.out=5), chain = ch1, tree = tree, dat = traits,
      SE=SE, prior=e_prior, startpar = startpar,
      new.dir = tempdir, ngen=ss_ngen, parallel = FALSE))
    ss2_cmd <- quote(steppingstone(
      Bk=seq(0,1, length.out=5), chain = ch2, tree = tree, dat = traits,
      SE=SE, prior=e_prior, startpar = startpar,
      new.dir = tempdir, ngen=ss_ngen, parallel = FALSE))
  }
  print(paste("MCMC chain 1", msg))
  t_ch1 <- eval(chain_cmd)
  print(paste("MCMC chain 2", msg))
  t_ch2 <- eval(chain_cmd)
  print("CHAINS FINISHED. SAVING.")
  # loading the chains and cleanup
  ch1 <- load.bayou(
    t_ch1, save.Rdata = T, cleanup = T, file = paste0(outdir,runid,"_1.rds"))
  ch1 <- set.burnin(ch1, burnin)
  ch2 <- load.bayou(
    t_ch2, save.Rdata = T, cleanup = T, file = paste0(outdir,runid,"_2.rds"))
  ch2 <- set.burnin(ch2, burnin)
  
  print(paste("Stepping Stone chain 1", msg))
  ss1 <- eval(ss1_cmd)
  print(paste("Stepping Stone chain 2", msg))
  ss2 <- eval(ss2_cmd)
  print("1")
  ss1 <- set.burnin(chain = ss1, burnin = burnin)
  print("2")
  lh1 <- ss1$lnr
  print("3")
  ss2 <- set.burnin(chain = ss2, burnin = burnin)
  print("4")
  lh2 <- ss2$lnr
  
  # there seems to be a bug in the code when only having t2 and theta to be estimated
  # bayou creates a list instead of avector which makes plotting difficult
  # as a workaround I just unlist here.
  if (typeof(ch1$theta) == "list") {
    print("YOLO")
    ch1$theta = unlist(ch1$theta)
    ch2$theta = unlist(ch2$theta)
    ch1$t2 = unlist(ch1$t2)
    ch2$t2 = unlist(ch2$t2)
    ss1$chains$theta = unlist(ss1$chains$theta)
    ss2$chains$theta = unlist(ss2$chains$theta)
    ss1$chains$t2 = unlist(ss1$chains$t2)
    ss2$chains$t2 = unlist(ss2$chains$t2)
    
    ch1$sb = rep(length(ch1$lnL), 0)
    ch2$sb = rep(length(ch2$lnL), 0)
    ch1$loc = rep(length(ch1$lnL), 0)
    ch2$loc = rep(length(ch2$lnL), 0)
  } 
  sink(paste0(outdir, runid, "_chain_summary.txt"))
  print(ss1)
  print(ss2)
  sink()
  
  # plot of prior followed by the two posteriors
  pdf(paste0(outdir, runid, "_prior_ss.pdf"))
  eval(prior)
  par(mfrow=c(1,1))
  gelman.R("lnL", chain1=ch1, chain2=ch2,
           plot=TRUE, type="n", ylim=c(0.9, 2))
  gelman.R("alpha", chain1=ch1, chain2=ch2,
           plot=TRUE, type="n", ylim=c(0.9, 2))
  gelman.R("sig2", chain1=ch1, chain2=ch2,
           plot=TRUE, type="n", ylim=c(0.9, 2))
  plot(ss1)
  plot(ss2)
  dev.off()
  
  # return the results as a list
  res = list(ch1=ch1, ch2=ch2, ss1=ss1, ss2=ss2, lh1=lh1, lh2=lh2)
  return(res)
}
