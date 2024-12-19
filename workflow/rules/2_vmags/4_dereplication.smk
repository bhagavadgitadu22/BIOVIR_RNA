# Dereplication of the viruses
rule blast_before_dereplication:
    output: 
        blast_results = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "blastn_not_derep_{completeness}_viruses_vs_all.tsv")
    input:
        fna_viruses_not_derep = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "{completeness}_viruses_no_duplicates.fna")
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    threads: config['blast']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "blast_all_vs_all_{completeness}_viruses.log")
    message: "BLAST all-vs-all of the viruses"
    shell:
        """
        (date && blastn -num_threads {threads} -query {input.fna_viruses_not_derep} -subject {input.fna_viruses_not_derep} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out {output.blast_results} &&
        date) &> {log}
        """

rule ani_for_dereplication:
    output: 
        ani_results = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "ani_not_derep_{completeness}_viruses_vs_all.tsv"),
        clustering_results = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "clusters_not_derep_{completeness}_viruses.tsv"),
    input:
        blast_results = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "blastn_not_derep_{completeness}_viruses_vs_all.tsv"), 
        fna_viruses_not_derep = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "{completeness}_viruses_no_duplicates.fna"),
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "ani_for_dereplication_{completeness}_viruses.log")
    message: "ANI for dereplication of the viral dataset"
    shell:
        """
        (date && python scripts/ani_calc.py -i {input.blast_results} -o {output.ani_results} &&
        python scripts/ani_clust.py --fna {input.fna_viruses_not_derep} --ani {output.ani_results} --out {output.clustering_results} --min_ani 95 --min_tcov 85 --min_qcov 0 && date) &> {log}
        """

rule viruses_dereplicated:
    output:
        list_viruses_derep = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "{completeness}_viruses_dereplicated.txt"),
        fna_viruses_derep = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "{completeness}_viruses_dereplicated.fna"),
    input:
        fna_viruses_not_derep = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "{completeness}_viruses_no_duplicates.fna"),
        clustering_results = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "clusters_not_derep_{completeness}_viruses.tsv"),
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "dereplication_{completeness}_viruses.log")
    message: "Dereplication of the viral dataset"
    shell:
        """
        (date && cut -f 1 {input.clustering_results} > {output.list_viruses_derep} &&
        seqtk subseq {input.fna_viruses_not_derep} {output.list_viruses_derep} > {output.fna_viruses_derep} && date) &> {log}
        """

# Checking similarities between viruses kept in case for instance some viruses are included in others
rule blast_dereplicated_viruses:
    output:
        blast_results = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "blastn_dereplicated_HQ_viruses_vs_all.tsv"),
    input: 
        fna_viruses_derep = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "HQ_viruses_dereplicated.fna")
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    threads: config['blast']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "blast_dereplicated_viruses.log")
    message: "Checking similarities of dereplicated viruses"
    shell:
        """(date && blastn -num_threads {threads} -query {input.fna_viruses_derep} -subject {input.fna_viruses_derep} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out {output.blast_results} && date) > {log}"""
