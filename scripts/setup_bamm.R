#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Incorrect number of arguments. \n 
       USAGE: setup_bamm.R tree_file trait_file", call.=FALSE)
}
library(BAMMtools)
tree = read.tree(args[1])
trait_tmp = read.table(args[2])

trait = trait_tmp$V2
names(trait) = trait_tmp$V1
setBAMMpriors(phy = tree,
              traits = trait,
              outfile = "config_files/temp_priors.bamm")