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
        os.path.join(RESULTS_DIR, "MAGs")
    output:
        directory(os.path.join(RESULTS_DIR, "PanPolyGFS")),
        os.path.join(DATA_DIR, "status/drep_mags.done")
    conda:
        os.path.join(ENV_DIR, "drep.yaml")
    threads:
        config["drep"]["threads"]
    shell:
        "checkm data setRoot /mnt/esb-storage-01/NOMIS/databases &&"
        "dRep dereplicate {RESULTS_DIR}/drep_mags -p {threads} -g {input}/*.fasta -sa 0.99 -comp 90 -con 5 &&"
        "touch status/drep_mags.done"


