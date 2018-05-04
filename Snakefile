"""
Author: Jeremias N. Brand
Affiliation: Unibas
Aim: simulation study comparing the use and behaviour of bayestraits, bayou and bamm when detecting shifts of rates in BM
Date: 20. Dez 2017
Run: snakemake -s Snakefile
Latest modification:
"""
# import statements
from snakemake.utils import report
from os.path import join
import os
import time


# Uses a yaml file as input
configfile: "config.yaml"

# Preparation------------------------------------------------------------------
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

# build the directories
dir_ls = ["data", "out", "out/bt", "out/bayou", "out/bamm", "out/plots"]
save_mkdir(dir_ls)

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

# Models for each tool
bayou_model_set = config["bayou_model_set"]

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


tmp_bayou_runs = []
for i in analysis_list:
    if "bayou" in i:
        tmp_bayou_runs.append(i)

# we have several models we run for each prior set so we need to modify the runs list

bayou_runs = []
for a in tmp_bayou_runs:
    for m in bayou_model_set:
        bayou_runs.append(a.split("_pr")[0] + "_pr" + str(m) + a.split("_pr")[1])


bt_runs = []
for i in analysis_list:
    if "bt" in i:
        bt_runs.append(i)

print(bayou_runs)

bamm_runs = []
for i in analysis_list:
    if "bamm" in i:
        bamm_runs.append(i)
# define global things here


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
        prior_f = "config_files/pr{prior}{n}.bayou",
        trait_f = "data/{sample}.trait"
    output:
        checkpoint =  "{sample}_pr{prior}{n}.bayou",
        gelman_lnl =  "out/bayou/{sample}_pr{prior}{n}.gelman.lnl.txt", 
        gelman_sig =  "out/bayou/{sample}_pr{prior}{n}.gelman.sig2.txt", 
        plot       =  "out/bayou/{sample}_pr{prior}{n}.ss.pdf",
        sum        =  "out/bayou/{sample}_pr{prior}{n}.summary.txt"
    log:
        "log/{sample}_pr{prior}{n}.bayou.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    benchmark: "bench/{sample}_pr{prior}{n}.bayou.txt"
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
        prior_null = "config_files/prNull{n}.bt"
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
        runs = bt_runs
    output:
        marginal = "out/bt/bt.marginalLH.tsv",
        log      = "out/bt/combined_logs.tsv",
        lhDelta  = "out/plots/bt_LhDelta_vs_NoTips.pdf",
        shiftLh  = "out/plots/bt_NoShifts_vs_Lh.pdf",
        shiftTip = "out/plots/bt_NoShifts_vs_NoTips.pdf",
        report   = "bt.report"
    log:
        "log/bt.log"
    params:
        #nsim = config["nsim"],
        #n_tips = config["n_tips"]
    threads: 1
    # using curly braces in snakemake requires escaping them by repeating
    shell:
        """
        for run in $( find ./out/bt -name *Stones.txt );do
            name=${{run##*/}}
            name=${{name%.Stones.txt}}
            echo -n ${{run##*/}}"\t"
            echo -n $name | tr "_" "\t"
            echo -n "\t"
            tail -1 $run | awk '{{print $NF}}'
        done > {output.marginal}

        # now we concatenate all the chains into one fuile for easy plotting in R
        for run in $( find ./out/bt -name *ch1.VarRates.txt );do
            name=${{run##*/}}
            echo $name
            base_name=${{name%_ch*.VarRates.txt}}
            echo $base_name
            name=${{name%.VarRates.txt}}
            for chain in $( find ./out/bt -name $base_name*VarRates.txt  );do
                chain_name=${{chain%.VarRates.txt}}
                chain_name=${{chain_name##*/}}
                chain_detail=$( echo -n $chain_name | tr "_" "\t" )
                sed -n -e '/It/,$p' $chain  | awk -F "\t" -v chain="$chain_name" -v chain_detail="$chain_detail" '{{printf chain "\t" chain_detail "\t" $i"\t"; print ""}}' | tail -n+2  >> "out/bt/"$base_name"_VarRatesCombined.txt"
                wait
            done
        done

        # combine the scheduels of all the files so we can plot the alpha and sigma values
        for run in $( find ./out/bt -name *ch1.Log.txt );do
            name=${{run##*/}}
            chain_name=${{name%.Log.txt}}
            chain_detail=$( echo -n $chain_name | tr "_" "\t" )
            sed -n -e '/Iteration\t/,$p' $run  | awk -F "\t" -v chain="$chain_name" -v chain_detail="$chain_detail" '{{printf chain "\t" chain_detail "\t" $i"\t"; print ""}}' | tail -n+2  >> {output.log}
        done

        touch {output.report}
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
        "out/simulated_trees_info.csv",
        marginal = "out/bt/bt.marginalLH.tsv",
        log      = "out/bt/combined_logs.tsv",
        lhDelta  = "out/plots/bt_LhDelta_vs_NoTips.pdf",
        shiftLh  = "out/plots/bt_NoShifts_vs_Lh.pdf",
        shiftTip = "out/plots/bt_NoShifts_vs_NoTips.pdf",
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
