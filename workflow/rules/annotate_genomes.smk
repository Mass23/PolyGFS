# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: bin_collect_mantis, bin_link_mantis, mantis_config, mantis_reformat_consensus, bin_folder_sample_mantis, bin_folder_sample

###########################
# default
rule annotate_genomes:
    input:
        os.path.join(RESULTS_DIR, "logs/mantis.done")
        #os.path.join(RESULTS_DIR, "logs/mantis_move.done"),
        #os.path.join(RESULTS_DIR, "mantis_bins")
    output:
        touch("status/annotate_genomes.done")


########################################
# checkpoint rules for collecting bins #
########################################
rule bin_collect_mantis:
    input:
        os.path.join(RESULTS_DIR, "MAGs"),
        os.path.join(RESULTS_DIR, "Genomes")
    output:
        directory(os.path.join(RESULTS_DIR, "mantis_bins"))
    shell:
        "mkdir {output} &&"
        "cp -v {input[0]}/*.fasta {output} && "
        "cp -v {input[1]}/*.fna {output}"

checkpoint move_bins:
    input:
        os.path.join(RESULTS_DIR, "mantis_bins")
    output:
        os.path.join(RESULTS_DIR, "logs/mantis_move.done")
    shell:
        "cd {input} && "
        "for f in *.fna; do mv \"$f\" \"${{f%.fna}}.fa\" ;done && "
        "for f in *.fasta; do mv \"$f\" \"${{f%.fasta}}.fa\" ; done && "
        "touch {output}"

#rule bin_link_mantis:
#    input:
#        os.path.join(RESULTS_DIR, "mantis_links/{i}.fa"),
#    output:
#        os.path.join(RESULTS_DIR, "mantis_bins/{i}.fa")
#    wildcard_constraints:
#        #sample="|".join(SAMPLES)
#        i="(\w\.)+"
#    shell:
#        "ln -vs {input} {output}"


####################
# rules for MANTIS #
####################
# Prokka 
rule prokka:
    input:
        rules.checkpoint move_bins.output,
#        os.path.join(RESULTS_DIR, "logs/mantis_move.done"),
        os.path.join(RESULTS_DIR, "mantis_bins/{i}.fa")
    output:
        FAA=os.path.join(RESULTS_DIR, "prokka/{i}.faa"),
        FFN=os.path.join(RESULTS_DIR, "prokka/{i}.ffn"),
        GFF=os.path.join(RESULTS_DIR, "prokka/{i}.gff")
    log:
        os.path.join(RESULTS_DIR, "logs/prokka_{i}.log")
    threads:
        config['prokka']['threads']
    conda:
        os.path.join(ENV_DIR, "prokka.yaml")
    message:
        "Running Prokka on mags"
    shell:
        "(date && prokka --outdir $(dirname {output.FAA}) --prefix $(echo ${basename {input[1]}} | cut -d. -f1) {input[1]} --cpus {threads} --force && date) &> >(tee {log})"

##########
# MANTIS #
rule list_mantis_bins:
    input:
        os.path.join(RESULTS_DIR, "mantis_bins/")
    output:
        os.path.join(DATA_DIR, "mantis_bins.txt")
    shell:
        "ls {input}/*.fa > {output}"

rule mantis_metadata:
    input:
        txt=rules.list_mantis_bins.output,
        FAA=rules.prokka.output
    output:
        os.path.join(RESULTS_DIR, "data/mantis_metadata.tsv")
    shell:
        "for fname in {input.txt} ; do echo echo \"${{fname}}\"\"    \"$(echo {input.FAA}) ; done > {output}"

# Mantis: create config file from IMP config
rule mantis_config:
    output:
        "mantis.config"
    message:
        "Creating config for MANTIS"
    run:
        with open(output[0], "w") as ofile:
            # default HMMs
            for hmm_name, hmm_path in config["mantis"]["default"].items():
                ofile.write("%s=%s\n" % (hmm_name, hmm_path))
            # weights
            for weights_name, weights_value in config["mantis"]["weights"].items():
                ofile.write("%s=%f\n" % (weights_name, weights_value))

                
 # for hmm_path in config["mantis"]["custom"]:
 #     ofile.write("custom_hmm=%s\n" % hmm_path)
 
# Mantis: protein annotation
# NOTE: check installation before use: python submodules/mantis/ check_installation
rule mantis_run:
    input:
        FAA=rules.prokka.output,
        config="mantis.config"
    output:
        os.path.join(RESULTS_DIR, "mantis/{i}/consensus_annotation.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{i}.analysis_mantis.log")
    threads:
        config['mantis']['cores']
    params:
         # The below acquires the number of cores to use from the snakemake launch command
#        CORES=int(os.environ.get("CORES"))
        cores=config['mantis']['hmmer_threads']
    conda:
        os.path.join(ENV_DIR, "mantis.yaml")
    message:
        "Annotation with MANTIS for mags"
    shell:
        "(date && python {config[mantis][path]}/ run_mantis -t {input.FAA} --output_folder $(dirname {output}) --mantis_config {input.config} --hmmer_threads {params.cores} --cores {threads} --memory {config[mantis][single_mem]} --kegg_matrix && date) &> >(tee {log})"

# Mantis: reformat consensus output (to be imported in Python/R)
rule mantis_reformat_consensus:
    input:
        rules.mantis_run.output
    output:
        os.path.join(RESULTS_DIR, "mantis/{i}/consensus_annotation.reformatted.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{i}.reformat.log")
    message:
        "Reformatting MANTIS output for mags"
    shell:
        # merge annotations after "|", remove "|" separator
        "(date && paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2 {input} | sed 's/^\\t//;s/\\t/;/g') > {output} && date) &> >(tee {log})"

# def bins_mantis(wildcards):
#    checkpoint_output = checkpoints.move_bins.get(**wildcards).output[0]
#    return expand(os.path.join(RESULTS_DIR, "mantis/{i}/consensus_annotation.reformatted.tsv"),i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa")).i)

def bins_mantis(wildcards):
    checkpoint_output = checkpoints.move_bins.get(**wildcards).output[0]
    return expand(os.path.join(RESULTS_DIR, "prokka/{i}.faa"),
    i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa")).i
    )

rule bin_folder_sample_mantis:
    input:
        bins_mantis
    output:
        touch(os.path.join(RESULTS_DIR, "logs/mantis.done"))

