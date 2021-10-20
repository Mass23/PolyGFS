# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
rule phylo_analysis:
    input:
        os.path.join(RESULTS_DIR,"gtotree_output")
    output:
        touch("status/phylo_analysis.done")

rule run_gtotree:
    input:
        #os.path.join(DATA_DIR, "status/drep_mags.done")
    output:
        directory(os.path.join(RESULTS_DIR, "gtotree_output"))
    conda:
        os.path.join(ENV_DIR, "gtotree.yaml")
    params:
        hmm_gtotree=config['gtotree']['hmm'],
        threads_gtotree=config['gtotree']['threads']
    shell:
        "gunzip {RESULTS_DIR}/Genomes/*.fna.gz &&"
        "ls {RESULTS_DIR}/drep_mags/dereplicated_genomes/*.fasta > drep_mags_paths.txt &&"
        "ls {RESULTS_DIR}/Genomes/*.fna > genomes_paths.txt &&"
        "cat drep_mags_paths.txt genomes_paths.txt > genomes_drep_mags_paths.txt &&"
        "GToTree -o {RESULTS_DIR}/gtotree_output -f genomes_drep_mags_paths.txt -H {params.hmm_gtotree} -j {params.threads_gtotree}"


