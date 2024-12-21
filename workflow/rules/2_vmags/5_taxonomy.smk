rule taxonomy_viruses:
    output: os.path.join(RESULTS_DIR, "viruses", "genomad", "HQ_viruses_dereplicated_annotate", "HQ_viruses_dereplicated_taxonomy.tsv")
    input: 
        viruses = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "HQ_viruses_dereplicated.fna"),
        db = os.path.join(RESULTS_DIR, "dbs", "genomad_db", "genomad_marker_metadata.tsv"),
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: config['genomad']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/coverm_filter.log")
    message: "Finding the taxonomies of the viruses with geNomad"
    shell:
        "(date && genomad end-to-end --threads {threads} --disable-nn-classification {input.viruses} $(dirname $(dirname {output})) $(dirname {input.db}) && date) &> {log}"
