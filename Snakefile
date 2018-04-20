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

sel_n_sim = range(1, config["sel_n_sim"] + 1)
sel_n_tips = config["sel_n_tips"]
sel_clade = config["sel_clade"]
sel_rate = config["sel_rate"]
sel_prior_set = config["sel_prior_set"]
sel_tool = config["sel_tool"]

# I don't expect that I will need to change these often
birth = config["birth"]
death = config["death"]


analysis_list = []
for s in sel_n_sim:
    for n in sel_n_tips:
        for c in sel_clade:
            for r in sel_rate:
                for p in sel_prior_set:
                    for t in sel_tool:
                        analysis_list.append(
                            "t"      +str(s)+ 
                            "_N"     +str(n)+ 
                            "_clade" +str(c)+
                            "_b"     +str(birth)+ 
                            "_d"     +str(death)+ 
                            "_rate"  +str(r)+ 
                            "_pr"    +str(p)+ 
                            "."      +str(t)
                            )


tree_list = []
for s in sel_n_sim:
    for n in sel_n_tips:
        for c in sel_clade:
            for r in sel_rate:
                tree_list.append("data/" +
                    "t"      +str(s)+ 
                    "_N"     +str(n)+ 
                    "_clade" +str(c)+
                    "_b"     +str(birth)+ 
                    "_d"     +str(death)+ 
                    "_rate"  +str(r)+
                    ".trait" 
                    )


bayou_runs = []
for i in analysis_list:
    if "bayou" in i:
        bayou_runs.append(i)

bt_runs = []
for i in analysis_list:
    if "bt" in i:
        bt_runs.append(i)

bamm_runs = []
for i in analysis_list:
    if "bamm" in i:
        bamm_runs.append(i)
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
        "out/simulated_trees_info.csv",
        tree_list
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

rule tre_sim_illustration:
    output:
        "fig/tree_sim_N100_clade0.2_b" +str(birth)+ "_d" +str(death)+ ".pdf" 
    script:
        "scripts/sim_illustration.R"


rule run_bayou:
    input:
        trait_f = "data/{sample}.trait",
        prior_f = "config_files/pr{n}.bayou"
    output:
        checkpoint = "{sample}_pr{n}.bayou",
        ch_sum = "out/bayou/{sample}_pr{n}.summary.txt"
    log:
        "log/{sample}_pr{n}.bayou.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    benchmark:
        "bench/{sample}_pr{n}.bayou.txt"
    script:
        "scripts/run_bayou.R"

rule bayou:
    input:
        bayou_runs
    output:
        "bayou.report"
    log:
        "log/bayou.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    shell:
        """
        touch bayou.report
        """

rule run_bt:
    input:
        trait_f = "data/{sample}.trait",
        prior_f = "config_files/pr{n}.bt",
        prior_null = "config_files/pr_null{n}.bt"
    output:
        checkpoint = "{sample}_pr{n}.bt"
        #ch_sum = "out/bt/{sample}_pr{n}.summary.txt"
    log:
        "log/{sample}_pr{n}.bt.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    benchmark:
        "bench/{sample}_pr{n}.bt.txt"
    script:
        "scripts/run_bt.py"


rule bt:
    input:
        bt_runs
    output:
        "bt.report"
    log:
        "log/bt.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    shell:
        """
        touch bt.report
        """


rule run_bamm:
    input:
        trait_f = "data/{sample}.trait",
        prior_f = "config_files/pr{n}.bamm"
    output:
        checkpoint = "{sample}_pr{n}.bamm"
        #ch_sum = "out/bt/{sample}_pr{n}.summary.txt"
    log:
        "log/{sample}_pr{n}.bamm.log"
    threads: 1
    benchmark:
        "bench/{sample}_pr{n}.bamm.txt"
    script:
        "scripts/run_bamm.py"


rule bamm:
    input:
        bamm_runs
    output:
        "bamm.report"
    log:
        "log/bamm.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    shell:
        """
        touch bamm.report
        """


rule report:
    input:
        "bt.report",
        "bayou.report",
        "bamm.report",
        "out/simulated_trees_info.csv"
        # "fig/tree_sim_N100_clade0.2_b" +str(birth)+ "_d" +str(death)+ ".pdf" 
    output:
        report_name         = "report.html"
    run:
        from snakemake.utils import report
        run_name = "testRUN"
        report("""
        Template for the report file {run_name}
        """, output[0], **input)



# Finishing up --------------------------------------------------------------
onsuccess:
    print("Workflow finished, no error")
onerror:
    print("An error occurred with the snakemake run")
    # here it would be good to include timstamping
    # shell("mail -s 'Error in Snakemake run' jeremias.brand@unibas.ch < {log}")

# TODO add separate snakefile for annotation
