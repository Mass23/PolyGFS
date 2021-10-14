# mags workflow
# taxonomy of the mags, and then mags selection 

localrules:

###########################
# default

rule collect_mags:
    input:
        os.path.join(RESULTS_DIR, "status/download_genomes.done"),
        os.path.join(RESULTS_DIR, "gtdbtk_output")
    output:
        touch("status/collect_mags.done")


################################################
# rule for collecting mags with taxo. matching #
################################################

rule gtdbtk:
    input:
        os.path.join(DATA_DIR, "genomes_list.txt"),
        os.path.join(MAGS_DIR)
    output:
        directory(os.path.join(RESULTS_DIR, "gtdbtk_output"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk_updated.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB toolkit on MAGs"
    shell:
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fasta --genome_dir {input[1]} --out_dir {output} && date) &> >(tee {log})"

