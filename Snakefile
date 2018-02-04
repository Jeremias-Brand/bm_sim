"""
Author: Jeremias N. Brand
Affiliation: Unibas
Aim: Workflow for transcriptome assembly, annotation and assessment
Date: 20. Nov. 2016
Run: snakemake   -s Snakefile
Latest modification:
  - add more flexible input
  - add reporting
  - add checks for gzip
  - integrate gzip of intermediate results
  - add a block to run if fail or succeed
  - add evironment file for bioconda
  - plotting of busco automated
  - move transcriptomes to better places
  - add trinity version to software.versions file
"""
# import statements
from snakemake.utils import report
from os.path import join
import os
import time


# Uses a yaml file as input
configfile: "config.yaml"

# Preparation------------------------------------------------------------------
#TIMESTAMP = time.strftime("%Y%m%d")
TIMESTAMP = "rcorrtest"
# check if the necessary dirs exist and if not creates them
def save_mkdir( dirs ):
	for d in dirs:
		if not os.path.isdir(d):
			os.mkdir(d)
			print("Creating directory: " + d)

def check_trailing_slash( path ):
    if path[-1] != "/":
        path = path + "/"
    return path



# Globals ---------------------------------------------------------------------

# define global things here


# Full path to a folder that holds all of your FASTQ files.
# FASTQ_DIR = check_trailing_slash(config["fastqdir"])
# # loads all of the path into variable
# TRANSRATE_DIR = check_trailing_slash(config["transrateDir"])
# HOME_DIR = check_trailing_slash(config["homedir"])
# VOUCHER_DIR = check_trailing_slash(config["voucher_dir"])


# A Snakemake regular expression matching the forward mate FASTQ files.
# we get a list with the names of all the files in the working directory



# Rules -----------------------------------------------------------------------



rule all:
    input:
        "report.html"


# TODO: the tree simulation must output a file that specifies where the shift is and also stores other params
rule tre_sim:
    output:
        "out/simulated_trees_info.tsv",
        str(config["n_sim"]) + "lelfile"
    log:
        "log/tre_sim.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    benchmark:
        "bench/tre_sim.txt"
    script:
        "scripts/tre_sim.R"





rule report:
    input:
       "out/simulated_trees_info.tsv"
    output:
        report_name         = "report.html"
    run:
        from snakemake.utils import report
        run_name = "testRUN"
        report("""
        Results of the macpipe_trans run for transcriptome assembly and assessment.
        Run timestamp/ID: {run_name}
        """, output[0], **input)



# Finishing up --------------------------------------------------------------
onsuccess:
    print("Workflow finished, no error")
onerror:
    print("An error occurred with the snakemake run")
    # here it would be good to include timstamping
    # shell("mail -s 'Error in Snakemake run' jeremias.brand@unibas.ch < {log}")

# TODO add separate snakefile for annotation
