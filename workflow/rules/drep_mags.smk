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
    output:
        os.path.join(RESULTS_DIR, "")
    shell:
        "cat genome_list.txt mags_list.txt > merged_genomes_mags_list.txt"
        "dRep dereplicate PanPolyGFS -p CORES -g merged_genomes_mags_list.txt -sa 0.95 -comp 75 -con 25"


