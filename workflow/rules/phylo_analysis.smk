# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
# default

hmm_gtotree = 'Gammaproteobacteria'

# checkpoint to collect representative sequences
checkpoint rep_genomes_collect_drep:
    input:
    output:
    shell:
        "ls DREP_OUT_DIR/dereplicated_genomes/*.fasta > rep_genomes.txt"

rule phylo_analysis:
    input:
        os.path.join(RESULTS_DIR, "logs/phylo_analysis.done")
    output:
    shell:
        "GToTree -a alteromonas_refseq_accessions.txt -f rep_genomes_list.txt -H {hmm_gtotree} -j 4 -o PanPolyGFS_tree"


