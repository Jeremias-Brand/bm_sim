
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/jeremias/.pyenv/versions/anaconda3-4.3.0/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X,\x00\x00\x00data/t1_N100_clade0.2_b1.5_d0.5_rate10.traitq\x06X\x13\x00\x00\x00config_files/pr1.btq\x07e}q\x08(X\x06\x00\x00\x00_namesq\t}q\n(X\x07\x00\x00\x00trait_fq\x0bK\x00N\x86q\x0cX\x07\x00\x00\x00prior_fq\rK\x01N\x86q\x0euh\x0bh\x06h\rh\x07ubX\x06\x00\x00\x00outputq\x0fcsnakemake.io\nOutputFiles\nq\x10)\x81q\x11X(\x00\x00\x00t1_N100_clade0.2_b1.5_d0.5_rate10_pr1.btq\x12a}q\x13(h\t}q\x14X\n\x00\x00\x00checkpointq\x15K\x00N\x86q\x16sh\x15h\x12ubX\x06\x00\x00\x00paramsq\x17csnakemake.io\nParams\nq\x18)\x81q\x19}q\x1ah\t}q\x1bsbX\t\x00\x00\x00wildcardsq\x1ccsnakemake.io\nWildcards\nq\x1d)\x81q\x1e(X!\x00\x00\x00t1_N100_clade0.2_b1.5_d0.5_rate10q\x1fX\x01\x00\x00\x001q e}q!(h\t}q"(X\x06\x00\x00\x00sampleq#K\x00N\x86q$X\x01\x00\x00\x00nq%K\x01N\x86q&uX\x06\x00\x00\x00sampleq\'h\x1fh%h ubX\x07\x00\x00\x00threadsq(K\x01X\t\x00\x00\x00resourcesq)csnakemake.io\nResources\nq*)\x81q+(K\x01K\x01e}q,(h\t}q-(X\x06\x00\x00\x00_coresq.K\x00N\x86q/X\x06\x00\x00\x00_nodesq0K\x01N\x86q1uh.K\x01h0K\x01ubX\x03\x00\x00\x00logq2csnakemake.io\nLog\nq3)\x81q4X0\x00\x00\x00log/t1_N100_clade0.2_b1.5_d0.5_rate10_pr1.bt.logq5a}q6h\t}q7sbX\x06\x00\x00\x00configq8}q9(X\x05\x00\x00\x00n_simq:KdX\x06\x00\x00\x00n_tipsq;]q<(K\nK2KdK\xc8M\x90\x01eX\x05\x00\x00\x00ratesq=]q>(K\x01K\x02K\x05K\nKdM\xe8\x03eX\n\x00\x00\x00clade_sizeq?]q@G?\xc9\x99\x99\x99\x99\x99\x9aaX\x05\x00\x00\x00birthqAG?\xf8\x00\x00\x00\x00\x00\x00X\x05\x00\x00\x00deathqBG?\xe0\x00\x00\x00\x00\x00\x00X\t\x00\x00\x00sel_n_simqCK\x01X\n\x00\x00\x00sel_n_tipsqD]qE(K\nK2KdK\xc8M\x90\x01eX\t\x00\x00\x00sel_cladeqF]qGG?\xc9\x99\x99\x99\x99\x99\x9aaX\x08\x00\x00\x00sel_rateqH]qIK\naX\r\x00\x00\x00sel_prior_setqJ]qKK\x01aX\x08\x00\x00\x00sel_toolqL]qM(X\x05\x00\x00\x00bayouqNX\x02\x00\x00\x00btqOX\x04\x00\x00\x00bammqPeX\x05\x00\x00\x00n_genqQM\x10\'X\x06\x00\x00\x00burninqRG?\xd3333333X\x08\x00\x00\x00n_chainsqSK\x02X\x0f\x00\x00\x00mcmc_write_freqqTM\xe8\x03X\x0f\x00\x00\x00mcmc_print_freqqUJ@B\x0f\x00X\x0f\x00\x00\x00mcmc_event_freqqVK\nuX\x04\x00\x00\x00ruleqWX\x06\x00\x00\x00run_btqXub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
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
