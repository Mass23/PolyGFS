# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
# default 
repos = ["genbank","refseq"]

rule dn_genomes:
    input:
        expand('Genomes/{repo}', genus=genera, repo=repos),
        gunzip_genomes(),
        copy_genomes()
    output:
        touch("status/annotate_genomes.done")

########################################
# checkpoint rules for collecting bins #
########################################
checkpoint bin_collect_mantis:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Binning/selected_DASTool_bins")
    output:
        directory(os.path.join(RESULTS_DIR, "mantis_links/{sample}_bins"))
    shell:
        "ln -vs {input} {output}"

rule bin_link_mantis:
    input:
        os.path.join(RESULTS_DIR, "mantis_links/{sample}_bins/{i}.contigs.fa"),
    output:
        os.path.join(RESULTS_DIR, "mantis_bins/{sample}__{i}.contigs.fa")
    wildcard_constraints:
        sample="|".join(SAMPLES)
        #i="(\w\.)+"
    shell:
        "ln -vs {input} {output}"

####################
# rules for ncbi-download-genomes #
####################

rule download_genomes:
    input:
        GENOMES
    output:
        directory(os.path.join(RESULTS_DIR,"Genomes/"))
    threads:
        1
    run:
        for line in open(input, 'r').readline():
            if line.startswith('GCA'):
                args = ['ncbi-genome-download','-o',output,'--formats','fasta','-s','genbank','-A',line,'all']
                subprocess.call(' '.join(args), shell = True)
            elif input.startswith('GCF'):
                args = ['ncbi-genome-download','-o',output,'--formats','fasta','-s','refseq','-A',line,'all']
                subprocess.call(' '.join(args), shell = True)
            else:
                print(str(input) + ' is not a refseq or genbank accession!')
