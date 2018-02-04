n_tips <- snakemake@config[["n_tips"]]

print(n_tips)

write.table(n_tips, "out/simulated_trees_info.tsv")