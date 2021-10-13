# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
# default 
rule dn_genomes:
    input:
        os.path.join(RESULTS_DIR,"Genomes")
    output:
        os.path.join(RESULTS_DIR, "genomes_list.txt")

# end it by writing the file "genomes_list.txt" that contains the ones that were collected

########################################
# rules for ncbi-download-genomes #
########################################

rule download_genomes:
    input:
        GENOMES
    output:
        directory(os.path.join(RESULTS_DIR,"Genomes")),
        os.path.join(RESULTS_DIR,"Genomes/download_genomes.done")
    threads:
        1
    conda:
        os.path.join(ENV_DIR, "ncbi-g-d.yaml")
    script:
        "dn_genomes_script.py"

rule create_mags_dir_file:
    input:
        os.path.join(RESULTS_DIR,"Genomes/download_genomes.done")
    output:
        os.path.join(RESULTS_DIR, "genomes_list.txt")
    shell:
        "$(basename -s '.fasta' $(ls $(dirname {input}))/*.fasta) > genomes_list.txt"

