from shutil import copyfile
import subprocess

n_gen  = snakemake.config["n_gen"]
burnin = snakemake.config["burnin"]
trait_f = snakemake.input["trait_f"]
prior_f = snakemake.input["prior_f"]
tree_f = trait_f.split("_rate")[0] + ".nex" 
checkpoint = snakemake.output["checkpoint"]

prefix = trait_f.split("/")[1]


outfile_n = "out/bt/" + prefix + ".bt.cmd"

print(tree_f)

bt_burnin = int(n_gen * burnin)

copyfile("/Users/jeremias/Dropbox/Bioinformatics/Reports/bm_sim/config_files/pr1.bt", outfile_n)

with open(outfile_n, "a") as o_f:
    o_f.write("Iterations " + str(n_gen) + "\n")
    o_f.write("Burnin " + str(bt_burnin) + "\n")
    o_f.write("logFile " + "./out/bt/" + prefix + "_bt" + "\n")
    o_f.write("SaveModels" + "\n")
    o_f.write("Run")
    
  
print(tree_f, trait_f, outfile_n)    
subprocess.call("BayesTraitsV3 " + tree_f + " "  + trait_f + " < " + outfile_n, shell=True)
