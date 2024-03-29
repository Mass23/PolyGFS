# mags workflow
# taxonomy of the mags, and then mags selection 

localrules:

###########################
# default

rule collect_mags:
    input:
        #os.path.join(DATA_DIR, "status/dn_genomes.done"),
        os.path.join(RESULTS_DIR, "gtdbtk_output"),
        os.path.join(DATA_DIR, "mags_list.txt"),
        os.path.join(RESULTS_DIR, "MAGs")
    output:
        touch("status/collect_mags.done")


################################################
# rule for collecting mags with taxo. matching #
################################################

rule gtdbtk:
    input:
        #os.path.join(DATA_DIR, "status/dn_genomes.done"),
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
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fasta --genome_dir {input[1]} --out_dir {output} && date) &> {log}"

rule list_target_mags:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk_output/gtdbtk.bac120.summary.tsv")
    output:
        os.path.join(DATA_DIR, "mags_list.txt")
    run:
        tax_string = 'g__Polaromonas'
        gtdbtk_file = pd.read_csv(input[0], sep='\t')
        gtdbtk_sub = gtdbtk_file[gtdbtk_file['classification'].str.contains(tax_string)].user_genome
        gtdbtk_sub.to_csv(output[0], sep='\t', index=False)

rule copy_target_mags:
    input:
        os.path.join(DATA_DIR, "mags_list.txt")
    output:
        directory(os.path.join(RESULTS_DIR, "MAGs"))
    shell:
        "mkdir {output} &&"
#        "for line in $(cat {input}); do cp -v {MAGS_DIR}/$line.fasta {output} ;done"
        "for line in $(tail -n+2 {input}); do cp -v {MAGS_DIR}/$line.fasta {output} ;done"

