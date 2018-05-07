library(ggplot2)  
library(ape)
library(tidyr)
# we need this package toprevent masking from screwing everything up
library(dplyr)
library(dtplyr)
library(data.table)
################################################################################

# Summarize the number of shifts
params <- 100
col_names <- c("file", "tree", "No_tips", "clade_proportion", "birth",
               "death", "rate_delta", "tool", "prior_set" ,"chain", "It", "Lh", "Lh_and_Prior", "No_Pram", "Alpha", "Sigma_square", "Scale_Prior",
               paste(rep(c("Node_ID", "Scaler", "Create_It", "Node_Branch"),params),
                     rep(1:params, each = 4)), sep = "_")

rates_f <- list.files("out/bt", pattern = "*VarRatesCombined.txt")

# combine all the files into one dataframe
# using the data.table version of this is actually at least one third faster!
res_ls <- lapply(paste0("out/bt/", rates_f),
       read.table, fill = TRUE, sep = "\t", col.names = col_names)
All_VarRes <- rbindlist(res_ls)


# All_VarRes <- do.call(rbind,
#                            lapply(paste0("out/bt/", rates_f),
#                                   read.table, fill = TRUE, sep = "\t", col.names = col_names))
All_VarRes$No_tips <- factor(All_VarRes$No_tips,
                   levels = c("N10", "N50", "N100", "N200", "N400"))

All_VarRes$rate_delta <- factor(All_VarRes$rate_delta,
                                levels = c("rate1.trait", "rate2.trait",
                                           "rate5.trait", "rate10.trait",
                                           "rate100.trait"))

var_sum <-  unite(All_VarRes, run, tree:tool, sep="_")

################################################################################

# First analysis plots

width_multiplier <- length(unique(All_VarRes$rate_delta))
pdf(file = "out/plots/bt_NoShifts_vs_NoTips.pdf", width = 7 * width_multiplier,
    height = 7)
p  <- ggplot(data = All_VarRes, aes(No_Pram, x = No_tips, color = factor(chain)))
p  <- p + geom_boxplot(alpha=0.55) + facet_grid(~ rate_delta) + theme_bw()
p
dev.off()

wrap_multiplier <- length(unique(var_sum$run))/2
pdf(file = "out/plots/bt_NoShifts_vs_Lh.pdf",
    width = 7 * wrap_multiplier, height = 7 * wrap_multiplier)
p  <- ggplot(data = var_sum, aes(x = factor(No_Pram), y = Lh, color=No_Pram))
p  <- p + geom_point(alpha=0.55) + facet_wrap(~ run, scales = "free") + theme_bw()
p
dev.off()

################################################################################


lh_cols <- c("file", "tree", "No_tips", "clade_proportion", "birth",
             "death", "rate_delta", "tool", "prior_set" ,"chain", "Lh")
marginal_lh <- read.table("out/bt/bt.marginalLH.tsv", col.names =  lh_cols)
marginal_lh$No_tips <- factor(marginal_lh$No_tips,
                              levels = c("N10", "N50", "N100", "N200", "N400"))


width_multiplier <- length(unique(All_VarRes$rate_delta))
pdf(file = "out/plots/bt_LhDelta_vs_NoTips.pdf", width = 7 * width_multiplier,
    height = 7)
p <- ggplot(marginal_lh)
p + geom_boxplot(aes(x=No_tips, y=Lh, color = prior_set), width = 0.25) +
  ylab("Marginal likelihood") + xlab("Tree size")+
  facet_grid(~ rate_delta) + theme_bw()
dev.off()

################################################################################

# Bayes Factor

################################################################################
# Summary table that combines marginal likelihood
marginal_lh_combined <- marginal_lh %>%
  unite(run, tree:tool, sep="_") %>%
  group_by(run, prior_set) %>%
    summarise(Lh = mean(Lh))

# calculate the BF for each contrast
df = marginal_lh_combined
BF_table <- data.frame(run = rep(NA,length(unique(df$run))),
                       BF = rep(NA,length(unique(df$run))))
i = 1
for (r in unique(df$run)){
  simpleMod = filter(df, run == r, prior_set == "prNull1")$Lh
  complexMod = filter(df, run == r, prior_set == "pr1")$Lh
  BF = 2*(complexMod - simpleMod)
  BF_table[i, 1] = r
  BF_table[i, 2] = BF
  i = i + 1
}

# output the BF table
write.table(x = BF_table, file = "out/bt/BF_bt_chainsCombined.tsv", row.names = F)


pdf(file = "out/plots/bt_BF_vs_NoTips.pdf", width = 7 * width_multiplier,
    height = 7)
BF_table %>%
  separate(run, sep = "_", c("tree", "No_tips", "clade_proportion", "birth",
                             "death", "rate_delta", "tool")) %>%
  mutate(No_tips = factor(No_tips, levels = c("N10", "N50", "N100", "N200", "N400")),
         rate_delta = factor(rate_delta, levels = c("rate1.trait", "rate2.trait",
                                                    "rate5.trait", "rate10.trait",
                                                    "rate100.trait"))) %>%
  ggplot() +
  geom_jitter(aes(x=No_tips, y=BF), width=0.1) +
  ylab("Bayes Factor") + 
  xlab("Tree size") +
  facet_grid(~ rate_delta, scales = "free") + 
  geom_hline(yintercept=10, linetype="dashed", color = "red", size=2) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=2) + 
  theme_bw()
dev.off()

pdf(file = "out/plots/bt_BF_vs_NoTips_scalefree.pdf", width = 7 * width_multiplier,
    height = 7)
BF_table %>%
  separate(run, sep = "_", c("tree", "No_tips", "clade_proportion", "birth",
                             "death", "rate_delta", "tool")) %>%
  mutate(No_tips = factor(No_tips, levels = c("N10", "N50", "N100", "N200", "N400")),
         rate_delta = factor(rate_delta, levels = c("rate1.trait", "rate2.trait",
                                                    "rate5.trait", "rate10.trait",
                                                    "rate100.trait"))) %>%
  ggplot() +
  geom_jitter(aes(x=No_tips, y=BF), width=0.1) +
  ylab("Bayes Factor") + 
  xlab("Tree size") +
  facet_wrap(~ rate_delta, scales = "free", nrow = 1) + 
  geom_hline(yintercept=10, linetype="dashed", color = "red", size=2) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=2) + 
  theme_bw()
dev.off()


################################################################################

# Tree Parsing
merge_branch_ls <- function(l){
  res = c(l[1], l[2], paste(l[3:length(l)], collapse = "_"))
  return(res)
}

translate <- function(x, dictionary){
  res = dictionary[as.numeric(x)+1,2]
  if (is.na(res)){
    return("")
    }
  return(res)
}

translate_branch <- function(s, dictionary, sep = "_"){
  res = s
  cand = unlist(strsplit(s[3], "_"))
  res[3] = paste(sapply(cand, translate, dictionary), collapse="_")
  return(res)
}
################################################################################
# The VarRates file contains information about the tree that we need to determin where bt puts shifts

lines <- readLines("out/bt/t1_N10_clade0.2_b1.5_d0.5_rate2.trait_bt_pr1_ch1.VarRates.txt")

n_tips = as.numeric(lines[[1]])
n_nodes = as.numeric(lines[[n_tips + 2]])

strsplit(lines[[2]],split = "\t")

tip_ls <- sapply(X = lines[1:n_tips+1],FUN = strsplit, "\t")
tip_df <- data.frame(matrix(unlist(tip_ls), nrow=n_tips, byrow=T))
#names(tip_df) <- c("bt_name", "sim_name")
tip_df <- transmute(tip_df, bt_name = as.character(X1), sim_name = as.character(X2))

node_ls <- sapply(X = lines[(n_tips + 3):(n_tips + n_nodes + 2)], FUN = strsplit, "\t")

node_ls2 <- lapply(node_ls, merge_branch_ls)
node_ls3 <- lapply(node_ls2, translate_branch, tip_df)
node_df <- data.frame(matrix(unlist(node_ls3), nrow=n_nodes, byrow=T))
names(node_df) <- c("bt_node", "length", "clade")


tree <- read.tree("data/t1_N10_clade0.2_b1.5_d0.5.nwk")
write.table("", file="bt.report")
