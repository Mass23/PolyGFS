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
configfile:"config.yaml"

##############################
# Relevant directories
DATA_DIR = config["work_dir"]
RESULTS_DIR = config["results_dir"]
ENV_DIR = config["env_dir"]
SRC_DIR = config["scripts_dir"]
BIN_DIR = config["bin_dir"]
FASTQ_DIR = config["fastq_dir"]

##############################
# Input
MAG = [line.strip() for line in open("polaromonas_mag_list").readlines()]
SAMPLES = [line.strip() for line in open("sample_list").readlines()]

##############################
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]


##############################
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam"), mag=MAG, sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "merged_bam/{mag}_merged.bam"), mag=MAG),
        expand(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz"), mag=MAG),
        expand(os.path.join(RESULTS_DIR, "coverage/{mag}_coverage.txt"), mag=MAG)


###########
# Mapping #
rule fa_index:
    input:
        ref_genome=os.path.join(BIN_DIR, "{mag}.fa"),
    output:
        expand(os.path.join(BIN_DIR, "{{mag}}.{ext}"), mag=MAG, ext=BWA_IDX_EXT)
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_index.log")
    threads:
        config["bwa"]["threads"]
    params:
        idx_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Indexing {wildcards.mag}"
    shell:
        "(date && bwa index {input} -p {params.idx_prefix} && date) &> {log}"

rule bwa_map:
    input:
        ref_genome=os.path.join(BIN_DIR, "{mag}.fa"),
        r1=os.path.join(FASTQ_DIR, "{sample}_R1.fastp.fastq.gz"),
        r2=os.path.join(FASTQ_DIR, "{sample}_R2.fastp.fastq.gz"),
        idx=expand(os.path.join(BIN_DIR, "{mag}.{ext}"), mag=MAG, ext=BWA_IDX_EXT)
    output:
        os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}.log")
    threads:
        config["bwa"]["threads"]
    params:
        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    wildcard_constraints:
        mag="|".join(MAG)
    message:
        "Mapping {wildcards.sample} onto {wildcards.mag}"
    shell:
        """(date && bwa mem -t {threads} -R @RG\\tID:{wildcards.sample}\\SM:{wildcards.sample} {input.ref_genome} {input.r1} {input.r2} | samtools view -Sb - > {output} && date) &> {log}"""

rule merge_bam:
    input:
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{{mag}}_{sample}.bam"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "merged_bam/{mag}_merged.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_merged_bam.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Merging the bam files for all samples for {wildcards.mag}"
    shell:
        "(date && samtools merge -rh {output} {input} && date) &> {log}"

rule bam2vcf:
    input:
        ref=os.path.join(BIN_DIR, "{mag}.fa"),
        bam=rules.merge_bam.output
    output:
        calls=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_calls.bcf")),
        filtered=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf")),
        gz=os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz")
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    message:
        "Calling SNPs on {wildcards.mag}"
    shell:
        "(date && bcftools mpileup --max-depth 10000 -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {output.calls} && "
        "bcftools view -i '%QUAL>=20' {output.calls} > {output.filtered} && "
        "bgzip {output.filtered} > {output.gz} date) &> {log}"

###########
# Depth 
rule depth:
    input:
        rules.merge_bam.output
    output:
        os.path.join(RESULTS_DIR, "coverage/{mag}_coverage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_coverage.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Estimating the coverage for all samples on {wildcards.mag}"
    shell:
        "(date && samtools depth {input} > {output} && date) &> {log}"
