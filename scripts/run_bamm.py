from shutil import copyfile
import subprocess

n_gen  = snakemake.config["n_gen"]
n_chains  = snakemake.config["n_chains"]
burnin = snakemake.config["burnin"]
mcmc_write_freq = snakemake.config["mcmc_write_freq"]
mcmc_print_freq = snakemake.config["mcmc_print_freq"]
mcmc_event_freq = int(n_gen/5000)
trait_f = snakemake.input["trait_f"]
prior_f = snakemake.input["prior_f"]
tree_f = trait_f.split("_rate")[0] + ".nwk" 
checkpoint = snakemake.output["checkpoint"]

prefix = trait_f.split("/")[1]
pr = prior_f.split("/")[1].split(".")[0]
prefix = prefix + "_" + pr


subprocess.call("scripts/setup_bamm.R " + tree_f + " "  + trait_f, shell=True)

outfile_n = "out/bamm/" + prefix + ".bamm.cmd"
copyfile("./config_files/pr1.bamm", outfile_n)
with open(outfile_n, "a") as o_f:
    o_f.write("treefile = " + tree_f + "\n")
    o_f.write("traitfile = " + trait_f + "\n")
    o_f.write("mcmcWriteFreq =  " + str(mcmc_write_freq) + "\n")
    o_f.write("printFreq =  " + str(mcmc_print_freq) + "\n")
    o_f.write("eventDataWriteFreq =  " + str(mcmc_event_freq) + "\n")
    o_f.write("numberOfGenerations = " + str(n_gen) + "\n")
    o_f.write("numberOfChains = " + str(n_chains) + "\n")
    o_f.write("mcmcOutfile = out/bamm/" + prefix + "_mcmc_out.txt" + "\n")
    o_f.write("eventDataOutfile = out/bamm/" + prefix + "_event_data.txt" + "\n")
    o_f.write("runInfoFilename = out/bamm/" + prefix + "_run_info.txt" + "\n")
    o_f.write("chainSwapFileName = out/bamm/" + prefix + "_chain_swap.txt" + "\n")
    with open("config_files/temp_priors.bamm", "r") as pr_f:
    	for line in pr_f:
    		o_f.write(line)
	    

subprocess.call("bamm -c " + outfile_n, shell=True)

with open(checkpoint, "w") as f:
	f.write("Placeholder for better things")
