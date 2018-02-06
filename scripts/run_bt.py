from shutil import copyfile
import subprocess

n_gen  = snakemake.config["n_gen"]
n_chains  = snakemake.config["n_chains"]
burnin = snakemake.config["burnin"]
trait_f = snakemake.input["trait_f"]
prior_f = snakemake.input["prior_f"]
tree_f = trait_f.split("_rate")[0] + ".nex" 
checkpoint = snakemake.output["checkpoint"]

prefix = trait_f.split("/")[1]
pr = prior_f.split("/")[1].split(".")[0]


print(tree_f)

bt_burnin = int(n_gen * burnin)

for i in range(1,n_chains + 1):
	outfile_n = "out/bt/" + prefix + ".bt.cmd"
	copyfile("/Users/jeremias/Dropbox/Bioinformatics/Reports/bm_sim/config_files/pr1.bt", outfile_n)

	with open(outfile_n, "a") as o_f:
	    o_f.write("Iterations " + str(n_gen) + "\n")
	    o_f.write("Burnin " + str(bt_burnin) + "\n")
	    o_f.write("logFile " + "./out/bt/" + prefix + "_bt_" + pr + "_ch" + str(i) + "\n")
	    o_f.write("SaveModels" + "\n")
	    o_f.write("Run")
	    
	subprocess.call("BayesTraitsV3 " + tree_f + " "  + trait_f + " < " + outfile_n, shell=True)

with open(checkpoint, "w") as f:
	f.write("Placeholder for better things")
