# plot_shift_bamm.R
library(ape)
library(geiger)
library(phytools)
library(BAMMtools)

# we need both the tree and the event data

runs <- unlist(
  strsplit(
    list.files(path = "data/", pattern = "*nwk"),
                        split = ".nwk"))

tree_name = "t1_N100_clade0.2_b1.5_d0.5"
rate_par = "rate10"
prior_par = "pr1"

tree <- read.tree(paste0("data/", tree_name, ".nwk"))
edata <- getEventData(tree, eventdata = paste0("out/bamm/", tree_name, "_", rate_par, ".trait_", prior_par, "_event_data.txt")
                     ,type= "trait"   , burnin=0.1)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1)
plot.credibleshiftset(css, lwd=1.7, plotmax=4)

shift_probs <- summary(edata)
 

for (run in runs[2]){
  print(run)
  tree <- read.tree(paste0("data/", run, ".nwk"))
  edata <- getEventData(tree, eventdata = paste0("out/bamm/", run, "_", rate_par, ".trait_", prior_par, "_event_data.txt")
                        ,type= "trait"   , burnin=0.99)
  css <- credibleShiftSet(edata, expectedNumberOfShifts=1)
  plot.credibleshiftset(css, lwd=1.7, plotmax=4)
}



