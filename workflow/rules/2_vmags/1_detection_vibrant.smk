# VIBRANT: output list of viruses including some proviruses whose names were changed compared to original contigs
rule db_vibrant:
    output: os.path.join(RESULTS_DIR, "logs/vibrant_db_downloaded.txt")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: config['vibrant']['threads']
    log: os.path.join(RESULTS_DIR, "logs/vibrant_db.log")
    message: "Downloading the VIBRANT database"
    shell: """(date && download-db.sh && echo "Database downloaded" > {output} && date) &> {log}"""

rule vibrant:
    output: os.path.join(RESULTS_DIR, "viruses", "vibrant", "VIBRANT_{assembly}.renamed.contigs", "VIBRANT_phages_{assembly}.renamed.contigs", "{assembly}.renamed.contigs.phages_combined.fna")
    input: 
        assembly = os.path.join(RESULTS_DIR, "megahit", "{assembly}_megahit", "{assembly}.renamed.contigs.fa"),
        #db = rules.db_vibrant.output,
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: config['vibrant']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_vibrant.log")
    message: "Running VIBRANT"
    shell:
        """(date && VIBRANT_run.py -t {threads} -i {input.assembly} -folder $(dirname $(dirname {output})) && date) &> {log}"""

# First CheckV
rule db_checkv:
    output: os.path.join(RESULTS_DIR, "dbs", "checkv-db-v1.5", "genome_db", "checkv_reps.dmnd")
    log: os.path.join(RESULTS_DIR, "logs/checkv_db.log")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    message: "Downloading the CheckV database"
    shell:
        """
        (date && mkdir -p $(dirname $(dirname {output})) && cd $(dirname $(dirname $(dirname {output}))) &&
        wget -nc https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz &&
        tar --skip-old-files -zxvf checkv-db-v1.5.tar.gz && 
        cd checkv-db-v1.5/genome_db && diamond makedb --in checkv_reps.faa --db checkv_reps && date) &> {log}
        """

rule first_checkv:
    output: 
        checkv_quality = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "{assembly}_first_checkv", "quality_summary.tsv"),
        checkv_proviruses = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "{assembly}_first_checkv", "proviruses.fna"),
    input: 
        vibrant = os.path.join(RESULTS_DIR, "viruses", "vibrant", "VIBRANT_{assembly}.renamed.contigs", "VIBRANT_phages_{assembly}.renamed.contigs", "{assembly}.renamed.contigs.phages_combined.fna"),
        db = rules.db_checkv.output,
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: config['checkv']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_first_checkv.log")
    message: "Running the first CheckV per assembly"
    shell:
        """(date && checkv end_to_end -t {threads} -d $(dirname $(dirname {input.db})) {input.vibrant} $(dirname {output.checkv_quality}) && date) &> {log}"""

# Extracting HQ (pro)viruses from checkV
rule list_proviruses_checkv:
    output: os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "proviruses.txt")
    input: os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "proviruses.fna")
    shell:
        """grep ">" {input} | sed 's/>//g' > {output} || touch {output}"""

rule list_HQ_checkv:
    output:
        list_HQ = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "HQ_viruses.txt"),
        temp_proviruses = temp(os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "temp_HQ_proviruses.txt")),
        list_HQ_proviruses = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "HQ_proviruses.txt"),
    input:
        quality = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "quality_summary.tsv"),
        proviruses = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "proviruses.txt"),
    message: "Extracting list of HQ (pro)viruses from checkV"
    shell:
        """
        grep -q -E "No.*High" {input.quality} && grep -E "No.*High" {input.quality} | cut -f 1 > {output.list_HQ} || echo "" > {output.list_HQ};
        
        grep -q -E "Yes.*High" {input.quality} && grep -E "Yes.*High" {input.quality} | cut -f 1,2 > {output.temp_proviruses} || echo "No contigs" > {output.temp_proviruses}
        sed -i 's/\t/|.*\//g' {output.temp_proviruses}
        grep -q -E -f {output.temp_proviruses} {input.proviruses} && grep -E -f {output.temp_proviruses} "$(dirname {input.quality})/proviruses.txt" > {output.list_HQ_proviruses} || echo "" > {output.list_HQ_proviruses};
        """

rule seq_viruses_checkv:
    output:
        fna_HQ = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "{completeness}_viruses.fna"),
        fna_HQ_proviruses = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "{completeness}_proviruses.fna"),
    input:
        quality = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "quality_summary.tsv"),
        list_HQ = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "{completeness}_viruses.txt"),
        list_HQ_proviruses = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "{assembly}_{number}_checkv", "{completeness}_proviruses.txt")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    message: "Extracting sequences of (pro)viruses from checkV"
    shell:
        """
        seqtk subseq "$(dirname {input.quality})/viruses.fna" {input.list_HQ} > {output.fna_HQ};
        sed -i 's/ /_/g' "$(dirname {input.quality})/proviruses.fna" && sed -i 's/ /_/g' {input.list_HQ_proviruses} &&
        seqtk subseq "$(dirname {input.quality})/proviruses.fna" {input.list_HQ_proviruses} > {output.fna_HQ_proviruses};
        """

rule regroup_viruses_checkv:
    output:
        all_viruses = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "all_{completeness}_viruses.fna"),
        all_proviruses = os.path.join(RESULTS_DIR, "viruses", "{number}_checkv", "all_{completeness}_proviruses.fna"),
    input:
        list_viruses = expand(os.path.join(RESULTS_DIR, "viruses", "{{number}}_checkv", "{assembly}_{{number}}_checkv", "{{completeness}}_viruses.fna"), assembly=SAMPLES_VIR),
        list_proviruses = expand(os.path.join(RESULTS_DIR, "viruses", "{{number}}_checkv", "{assembly}_{{number}}_checkv", "{{completeness}}_proviruses.fna"), assembly=SAMPLES_VIR),
    message: "Concatenating viruses and proviruses from CheckV"
    shell:
        """
        cat {input.list_viruses} > {output.all_viruses}
        cat {input.list_proviruses} > {output.all_proviruses}
        """

# Removal of complete viruses and proviruses from VIBRANT results for COBRA
rule prep_cobra:
    output:
        list_cobra = os.path.join(RESULTS_DIR, "viruses", "prep_cobra", "{assembly}_for_COBRA.txt"),
        fna_cobra = os.path.join(RESULTS_DIR, "viruses", "prep_cobra", "{assembly}_for_COBRA.fna"),
    input: 
        quality = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "{assembly}_first_checkv", "quality_summary.tsv"),
        fna = os.path.join(RESULTS_DIR, "viruses", "vibrant", "VIBRANT_{assembly}.renamed.contigs", "VIBRANT_phages_{assembly}.renamed.contigs", "{assembly}.renamed.contigs.phages_combined.fna"),
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    message: "Preparing COBRA input contigs"
    shell:
        """
        grep -v "_fragment" {input.quality} | grep -v "Complete" | tail +2 | cut -f 1 > {output.list_cobra} && 
        seqtk subseq {input.fna} {output.list_cobra} > {output.fna_cobra}
        """

# Statistics of viruses obtained with first CheckV
rule statistics_first_checkv:
    output: os.path.join(RESULTS_DIR, "statistics", "statistics_post_VIBRANT.csv")
    input:
        vibrant_list = expand(os.path.join(RESULTS_DIR, "viruses", "vibrant", "VIBRANT_{assembly}.renamed.contigs", "VIBRANT_phages_{assembly}.renamed.contigs", "{assembly}.renamed.contigs.phages_combined.fna"), assembly=SAMPLES_VIR),
        checkv_list = expand(os.path.join(RESULTS_DIR, "viruses", "first_checkv", "{assembly}_first_checkv", "quality_summary.tsv"), assembly=SAMPLES_VIR),
        cobra_list = expand(os.path.join(RESULTS_DIR, "viruses", "prep_cobra", "{assembly}_for_COBRA.fna"), assembly=SAMPLES_VIR),
    params: 
        assemblies = SAMPLES_VIR,
        results_dir = RESULTS_DIR,
    message: "Computing statistics of the first CheckV"
    shell:
        """
        echo "sample;#all_vibrant;#proviruses_vibrant;#complete_viruses_vibrant;#complete_proviruses_vibrant;#HQ_viruses_vibrant;#HQ_proviruses_vibrant;#contigs_for_COBRA" > {output}
        for ass in {params.assemblies}; do
            vibrant="{params.results_dir}/viruses/vibrant/VIBRANT_$(echo $ass).renamed.contigs/VIBRANT_phages_$(echo $ass).renamed.contigs/$(echo $ass).renamed.contigs.phages_combined.fna"
            checkv="{params.results_dir}/viruses/first_checkv/$(echo $ass)_first_checkv/quality_summary.tsv"
            cobra="{params.results_dir}/viruses/prep_cobra/$(echo $ass)_for_COBRA.fna"

            n_vibrant=$(grep -q ">" $vibrant && grep ">" $vibrant | wc -l || echo 0);
            n_vibrant_proviruses=$(grep -q "_fragment" $vibrant && grep "_fragment" $vibrant | wc -l || echo 0); 
            n_complete=$(grep -q "Complete" $checkv && grep "Complete" $checkv | wc -l || echo 0);
            n_complete_proviruses=$(grep -q -E "(_fragment|\|provirus)" $checkv && grep -E "(_fragment|\|provirus)" $checkv | grep -q "Complete" && grep "Complete" | wc -l || echo 0);
            n_high_quality=$(grep -q -P "High-quality\tHigh-quality" $checkv && grep -P "High-quality\tHigh-quality" $checkv | wc -l || echo 0);
            n_high_quality_proviruses=$(grep -q -E "(_fragment|\|provirus)" $checkv | grep -q -P "High-quality\tHigh-quality" && grep -P "High-quality\tHigh-quality" | wc -l || echo 0); 
            ncobra=$(grep -q ">" $cobra && grep ">" $cobra | wc -l || echo 0); 
            echo "$ass;$n_vibrant;$n_vibrant_proviruses;$n_complete;$n_complete_proviruses;$n_high_quality;$n_high_quality_proviruses;$ncobra"
        done >> {output}
        """
