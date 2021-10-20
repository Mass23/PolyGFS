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
# SAMPLES = [line.strip() for line in open("sample_list").readlines()]
SED_GI = [line.strip() for line in open("sed_gi_list").readlines()]
SAMPLES = [line.strip() for line in open("sample_list").readlines()]
MICRO_FQ = [line.strip() for line in open("microthrix_fastq_list").readlines()]
ISOLATES = [line.strip() for line in open("microthrix_isolates_list").readlines()]

##############################
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]


##############################
rule all:
    input:
#        expand(os.path.join(BIN_DIR, "{mag}.bwt"), mag=ISOLATES),
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam"), mag=ISOLATES, sample=MICRO_FQ),
#        expand(os.path.join(RESULTS_DIR, "merged_bam/{mag}_merged.bam"), mag=MAG),
#        expand(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz"), mag=MAG),
#        expand(os.path.join(RESULTS_DIR, "coverage/{mag}_depth.cov"), mag=SED_GI),
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.genome.txt"), mag=ISOLATES, sample=MICRO_FQ),
        expand(os.path.join(RESULTS_DIR, "microthrix_coverage/{mag}_{sample}_depth.cov"), mag=ISOLATES, sample=MICRO_FQ),
        expand(os.path.join(RESULTS_DIR, "microthrix_vcf/{mag}_{sample}_filtered.bcf.gz"), mag=ISOLATES, sample=MICRO_FQ)


###########
# Mapping #
rule fa_index:
    input:
        ref_genome=os.path.join(BIN_DIR, "{mag}.fa"),
    output:
        expand(os.path.join(BIN_DIR, "{{mag}}.{ext}"), ext=BWA_IDX_EXT)
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
        r1=os.path.join(FASTQ_DIR, "{sample}/mg.r1.preprocessed.fq"),
        r2=os.path.join(FASTQ_DIR, "{sample}/mg.r2.preprocessed.fq"),
        idx=expand(os.path.join(BIN_DIR, "{mag}.{ext}"), mag=ISOLATES, ext=BWA_IDX_EXT)
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
        mag="|".join(ISOLATES),
        sample="|".join(MICRO_FQ)
    message:
        "Mapping {wildcards.sample} onto {wildcards.mag}"
    shell:
        """(date && bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' $(echo {input.ref_genome} | sed 's/.fa//g') {input.r1} {input.r2} | samtools view -Sb -F 4 - > {output} && date) &> {log}"""

rule sort_index_bam:
    input:
        os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam")
    output:
        sorted=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam"),
        index=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam.bai")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_bam_index.log")
    threads:
        config["bwa"]["sort"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    params:
        chunk_size=config["bwa"]["sort"]["chunk_size"]
    message:
        "Sorting and indexing bam files for {wildcards.mag} and {wildcards.sample}"
    shell:
        "(date && samtools sort --threads {threads} -m {params.chunk_size} {input} > {output.sorted} && "
        "samtools index {output.sorted} && date) &> {log}"

rule merge_bam:
    input:
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{{mag}}_{sample}.sorted.bam"), sample=MICRO_FQ)
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
        gz=os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_bam2vcf.log")
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    threads:
        config["bwa"]["vcf"]["threads"]
    params:
        calls=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_calls.bcf")),
        filtered=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf"))
    message:
        "Calling SNPs on {wildcards.mag}"
    shell:
        "(date && mkdir -p $(dirname {output.gz}) && bcftools mpileup -a FORMAT/AD,FORMAT/DP --threads {threads} --max-depth 10000 -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {params.calls} && "
        "bcftools view -i '%QUAL>=20' {params.calls} > {params.filtered} && "
        "bgzip {params.filtered} > {output.gz} && date) &> {log}"

rule sed_GI_bam2vcf:
    input:
        ref=os.path.join(BIN_DIR, "{mag}.fa"),
        bam=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam")
    output:
        gz=os.path.join(RESULTS_DIR, "microthrix_vcf/{mag}_{sample}_filtered.bcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_sed_GI_vcf.log")
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    threads:
        config["bwa"]["vcf"]["threads"]
    wildcard_constraints:
        mag="|".join(ISOLATES)
    params:
        calls=temp(os.path.join(RESULTS_DIR, "microthrix_vcf/{mag}_{sample}_calls.bcf")),
        filtered=temp(os.path.join(RESULTS_DIR, "microthrix_vcf/{mag}_{sample}_filtered.bcf"))
    message:
        "Calling SNPs on {wildcards.mag} and {wildcards.sample}"
    shell:
        "(date && mkdir -p $(dirname {output.gz}) && bcftools mpileup -a FORMAT/AD,FORMAT/DP --threads {threads} --max-depth 10000 -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {params.calls} && "
        "bcftools view -i '%QUAL>=20' {params.calls} > {params.filtered} && "
        "bgzip -f {params.filtered} > {output.gz} && date) &> {log}"


###########
# Depth
rule genome_file:
    input:
        os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam")
    output:
        os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.genome.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_genome.log")
    threads:
        config["bwa"]["sort"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Creating GENOME file for bedtools for {wildcards.mag} and {wildcards.sample}"
    shell:
        "(date && samtools view -H {input} | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > {output} && date) &> {log}"

rule depth:
    input:
        gen_file=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.genome.txt"),
        bam=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam")
    output:
        wind=temp(os.path.join(RESULTS_DIR, "microthrix_coverage/{mag}_{sample}_SlidingWindows1kb.bed")),
        depth=os.path.join(RESULTS_DIR, "microthrix_coverage/{mag}_{sample}_depth.cov")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_coverage.log")
    threads:
        config["bwa"]["vcf"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bedtools.yaml")
    wildcard_constraints:
        mag="|".join(ISOLATES)
    message:
        "Estimating the coverage for {wildcards.mag} and {wildcards.sample}"
    shell:
        "(date && mkdir -p $(dirname {output.depth}) && bedtools makewindows -g {input.gen_file} -w 1000 > {output.wind} && "
        "bedtools multicov -bams {input.bam} -bed {output.wind} > {output.depth} && date) &> {log}"



############ updated rules ##########
# rule picard_mark_dup:
#     input:
#         os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam")
#     output:
#         BAM=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.marked_dup.bam"),
#         METRICS=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.marked_dupMetrics.txt")
#     log:
#         os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_picard.log"
#     conda:
#         os.path.join(ENV_DIR, "gatk.yaml")
#     message:
#         "Marking duplicate reads in the BAM files for {wildcards.mag} with {wildcards.sample}"
#     shell:
#         "(date && java -jar picard.jar MarkDuplicates I={input} O={output.BAM} M={output.METRICS} && date) &> {log}"

# rule gatk_realign_indels:
#     input:
#         ref=os.path.join(BIN_DIR, "{mag}.fa"),
#         BAM=rules.picard_mark_dup.output.BAM, 
#     output:
#         intervals=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}_realigner.intervals")
#         gatkBAM=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.marked_realigned.bam")
#     log:
#         os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_GATK.log"
#     conda:
#         os.path.join(ENV_DIR, "gatk.yaml")
#     message:
#         "Realigning indels based on GATK for {wildcards.mag} with {wildcards.sample}"
#     shell:
#         "(date && "
#         "java –jar GenomeAnalysisTK.jar –T RealignerTargetCreator –R {input.ref} –I {input.BAM} –o {output.intervals} && "
#         "java –jar GenomeAnalysisTK.jar –T IndelRealigner –R {input.ref} –I {input.BAM} –targetIntervals {output.intervals} –o {output.gatkBAM} && "
#         "&& date) &> {log}"


# ###########
# # Depth
# rule genome_file:
#     input:
#         rules.gatk_realign_indels.output.gatkBAM
#     output:
#         os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.genome.txt")
#     log:
#         os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_genome.log")
#     threads:
#         config["bwa"]["sort"]["threads"]
#     conda:
#         os.path.join(ENV_DIR, "bwa.yaml")
#     message:
#         "Creating GENOME file for bedtools for {wildcards.mag} and {wildcards.sample}"
#     shell:
#         "(date && samtools view -H {input} | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > {output} && date) &> {log}"

# rule depth:
#     input:
#         gen_file=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.genome.txt"),
#         bam=rules.gatk_realign_indels.output.gatkBAM
#     output:
#         wind=temp(os.path.join(RESULTS_DIR, "microthrix_coverage/{mag}_{sample}_SlidingWindows1kb.bed")),
#         depth=os.path.join(RESULTS_DIR, "microthrix_coverage/{mag}_{sample}_depth.cov")
#     log:
#         os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_coverage.log")
#     threads:
#         config["bwa"]["vcf"]["threads"]
#     conda:
#         os.path.join(ENV_DIR, "bedtools.yaml")
#     wildcard_constraints:
#         mag="|".join(ISOLATES)
#     message:
#         "Estimating the coverage for {wildcards.mag} and {wildcards.sample}"
#     shell:
#         "(date && mkdir -p $(dirname {output.depth}) && bedtools makewindows -g {input.gen_file} -w 1000 > {output.wind} && "
#         "bedtools multicov -bams {input.bam} -bed {output.wind} > {output.depth} && date) &> {log}"
