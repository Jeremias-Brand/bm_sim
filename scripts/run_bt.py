from shutil import copyfile
import subprocess

n_gen  = snakemake.config["n_gen"]
n_chains  = snakemake.config["n_chains"]
bt_no_stones = snakemake.config["bt_no_stones"]
bt_iter_stones = snakemake.config["bt_iter_stones"]
burnin = snakemake.config["burnin"]
trait_f = snakemake.input["trait_f"]
prior_f = snakemake.input["prior_f"]
prior_null = snakemake.input["prior_null"]
tree_f = trait_f.split("_rate")[0] + ".nex" 
checkpoint = snakemake.output["checkpoint"]

prefix = trait_f.split("/")[1]
pr = prior_f.split("/")[1].split(".")[0]
pr_null = prior_null.split("/")[1].split(".")[0]


print(tree_f)

bt_burnin = int(n_gen * burnin)

for i in range(1,n_chains + 1):
    # we create the config files to run bt
    # for each run we have a file with the varrates turned on and corresponding null model run
	pr_file = "out/bt/" + prefix + "_bt_" + pr + "_ch" + str(i) +".cmd"
	null_file = "out/bt/" + prefix + "_bt_" + pr_null + "_ch" + str(i) +".cmd"

	copyfile("./" + prior_f, pr_file)
	copyfile("./" + prior_null, null_file)

	with open(pr_file, "a") as o_f:
	    o_f.write("Iterations " + str(n_gen) + "\n")
	    o_f.write("Burnin " + str(bt_burnin) + "\n")
	    o_f.write("logFile " + "./out/bt/" + prefix + "_bt_" + pr + "_ch" + str(i) + "\n")
	    o_f.write("Stones " + str(bt_no_stones) + " " + str(bt_iter_stones) + "\n")
	    o_f.write("SaveModels" + "\n")
	    o_f.write("Run")
	    
	with open(null_file, "a") as o_f:
	    o_f.write("Iterations " + str(n_gen) + "\n")
	    o_f.write("Burnin " + str(bt_burnin) + "\n")
	    o_f.write("logFile " + "./out/bt/" + prefix + "_bt_" + pr_null + "_ch" + str(i) + "\n")
	    o_f.write("Stones " + str(bt_no_stones) + " " + str(bt_iter_stones) + "\n")
	    o_f.write("SaveModels" + "\n")
	    o_f.write("Run")

	subprocess.call("BayesTraitsV3 " + tree_f + " "  + trait_f + " < " + pr_file, shell=True)
	subprocess.call("BayesTraitsV3 " + tree_f + " "  + trait_f + " < " + null_file, shell=True)

with open(checkpoint, "w") as f:
	f.write("Placeholder for better things")
