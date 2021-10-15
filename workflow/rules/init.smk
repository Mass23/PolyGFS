# Pipeline for mapping reads to bins and then identify SNPs
#
# Example call: snakemake -s Snakefile --configfile config.yaml --use-conda --cores 1 -rpn

##############################
# MODULES
import os, re
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO

##############################
# CONFIG
# can be overwritten by using --configfile <path to config> when calling snakemake
configfile:"config/config.yaml"

##############################
# Relevant directories
DATA_DIR = config["work_dir"]
RESULTS_DIR = config["results_dir"]
ENV_DIR = config["env_dir"]
SRC_DIR = config["scripts_dir"]
MAGS_DIR = config["mags_dir"]
FASTQ_DIR = config["fastq_dir"]

##############################
# Input
#MAGS = [line.strip() for line in open("mags_list.txt").readlines()]
SAMPLES = [line.strip() for line in open("samples_list.txt").readlines()]
GENOMES = "accessions_list.txt"

##############################
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]


##############################
#rule all:
#    input:
#        expand(os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam"), mag=SED_GI, sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "merged_bam/{mag}_merged.bam"), mag=MAG),
#        expand(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz"), mag=MAG),
#        expand(os.path.join(RESULTS_DIR, "coverage/{mag}_depth.cov"), mag=SED_GI),
#        expand(os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.genome.txt"), sed_gi=SED_GI, sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "sed_gi_coverage/{sed_gi}_{sample}_depth.cov"), sed_gi=SED_GI, sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "sed_gi_vcf/{sed_gi}_{sample}_filtered.bcf.gz"), sed_gi=SED_GI, sample=SAMPLES)

