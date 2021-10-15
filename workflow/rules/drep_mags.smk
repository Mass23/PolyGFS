# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
# default        
rule drep_mags:
    input:
        os.path.join(RESULTS_DIR, "drep_mags")
    output:
        touch("status/drep_mags.done")

rule run_drep:
    input:
        os.path.join(RESULTS_DIR, "MAGs")
    output:
        directory(os.path.join(RESULTS_DIR, "drep_mags"))
    conda:
        os.path.join(ENV_DIR, "drep.yaml")
    threads:
        config["drep"]["threads"]
    params:
        config["checkm"]["db_path"]
    shell:
        "checkm data setRoot {params} &&"
        "dRep dereplicate {output} -p {threads} -g {input}/*.fasta -sa 0.99 -comp 90 -con 5"


