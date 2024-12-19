# Merge COBRA results for a CheckV quality check
rule concat_viruses_cobra:
    output: os.path.join(RESULTS_DIR, "viruses", "cobra", "{assembly}_viruses_post_cobra.fna")
    input: os.path.join(RESULTS_DIR, "viruses", "cobra", "{assembly}_cobra", "COBRA_joining_summary.txt")
    shell:
        """
        dir=$(dirname {input})
        cat $dir/COBRA_category_i_self_circular.fasta > {output}
        cat $dir/COBRA_category_ii-a_extended_circular_unique.fasta >> {output}
        cat $dir/COBRA_category_ii-b_extended_partial_unique.fasta >> {output}
        cat $dir/COBRA_category_ii-c_extended_failed.fasta >> {output}
        cat $dir/COBRA_category_iii_orphan_end.fasta >> {output}
        """

rule second_checkv:
    output:
        checkv_quality = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "quality_summary.tsv"),
        checkv_proviruses = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "proviruses.fna"),
    input:
        cobra = os.path.join(RESULTS_DIR, "viruses", "cobra", "{assembly}_viruses_post_cobra.fna"),
        db = rules.db_checkv.output,
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: config['checkv']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_second_checkv.log")
    message: "Running the second CheckV per assembly"
    shell:
        """(date && checkv end_to_end -t {threads} -d $(dirname $(dirname {input.db})) {input.cobra} $(dirname {output.checkv_quality}) && date) &> {log}"""

# Constitution of final viral datasets
rule final_vmag_dataset:
    output:
        viruses_with_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "HQ_viruses_with_duplicates.fna"),
        viruses_no_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "HQ_viruses_no_duplicates.fna"),
    input:
        viruses_checkv1 = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "all_HQ_viruses.fna"),
        proviruses_checkv1 = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "all_HQ_proviruses.fna"),
        viruses_checkv2 = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "all_HQ_viruses.fna"),
        proviruses_checkv2 = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "all_HQ_proviruses.fna")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    message: "Establishing final list of viruses (not dereplicated yet)"
    shell:
        """
        cat {input.viruses_checkv1} {input.proviruses_checkv1} {input.viruses_checkv2} {input.proviruses_checkv2} > {output.viruses_with_duplicates} &&
        seqkit rmdup -n {output.viruses_with_duplicates} > {output.viruses_no_duplicates}
        """

# Extract contigs >= 5kb with at least one viral gene according to checkV for second viral dataset
rule list_viral_contigs_checkv2:
    output:
        list_all = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "all_viruses.txt"),
        list_all_proviruses = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "all_proviruses.txt"),
    input:
        quality = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "quality_summary.tsv"),
        proviruses = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "proviruses.txt"),
    message: "Extracting list of HQ (pro)viruses from second checkV"
    shell:
        """
        awk '$1>=2' {input.quality} | awk '$2>=5000' | awk '$3=="No"' > {output.list_all};
        awk '$1>=2' {input.quality} | awk '$2>=5000' | awk '$3=="Yes"' > {output.list_all_proviruses};
        """

rule final_viral_contigs_dataset:
    output:
        viruses_with_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "all_viruses_with_duplicates.fna"),
        viruses_no_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "all_viruses_no_duplicates.fna"),
    input:
        viruses_checkv1 = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "all_HQ_viruses.fna"),
        proviruses_checkv1 = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "all_HQ_proviruses.fna"),
        viruses_checkv2 = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "all_HQ_viruses.fna"),
        proviruses_checkv2 = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "all_HQ_proviruses.fna")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    message: "Establishing final list of viruses (not dereplicated yet)"
    shell:
        """
        cat {input.viruses_checkv1} {input.proviruses_checkv1} {input.viruses_checkv2} {input.proviruses_checkv2} > {output.viruses_with_duplicates} &&
        seqkit rmdup -n {output.viruses_with_duplicates} > {output.viruses_no_duplicates}
        """

# Statistics of viruses obtained with second CheckV
rule statistics_second_checkv:
    output: os.path.join(RESULTS_DIR, "statistics", "statistics_post_COBRA.csv")
    input:
        vibrant_list = expand(os.path.join(RESULTS_DIR, "viruses", "vibrant", "VIBRANT_{assembly}.renamed.contigs", "VIBRANT_phages_{assembly}.renamed.contigs", "{assembly}.renamed.contigs.phages_combined.fna"), assembly=SAMPLES_VIR),
        cobra_list = expand(os.path.join(RESULTS_DIR, "viruses", "cobra", "{assembly}_viruses_post_cobra.fna"), assembly=SAMPLES_VIR),
        checkv_list = expand(os.path.join(RESULTS_DIR, "viruses", "second_checkv", "{assembly}_second_checkv", "quality_summary.tsv"), assembly=SAMPLES_VIR),
    params:
        assemblies = SAMPLES_VIR,
        results_dir = RESULTS_DIR,
    message: "Computing statistics of the second CheckV"
    shell:
        """
        echo "sample;#all_vibrant;#all_cobra;#proviruses_cobra;#complete_viruses_cobra;#complete_proviruses_cobra;#HQ_viruses_cobra;#HQ_proviruses_cobra" > {output}
        for ass in {params.assemblies}; do
            vibrant="{params.results_dir}/viruses/vibrant/VIBRANT_$(echo $ass).renamed.contigs/VIBRANT_phages_$(echo $ass).renamed.contigs/$(echo $ass).renamed.contigs.phages_combined.fna"
            cobra="{params.results_dir}/viruses/cobra/$(echo $ass)_viruses_post_cobra.fna"
            checkv="{params.results_dir}/viruses/second_checkv/$(echo $ass)_second_checkv/quality_summary.tsv"

            n_vibrant=$(grep -q ">" $vibrant && grep ">" $vibrant | wc -l || echo 0);
            n_cobra=$(grep -q ">" $cobra && grep ">" $cobra | wc -l || echo 0);
            n_cobra_proviruses=$(grep -q "_fragment" $cobra && grep "_fragment" $cobra | wc -l || echo 0);
            n_complete=$(grep -q "Complete" $checkv && grep "Complete" $checkv | wc -l || echo 0);
            n_complete_proviruses=$(grep -q -E "(_fragment|\|provirus)" $checkv && grep -E "(_fragment|\|provirus)" $checkv | grep -q "Complete" && grep "Complete" | wc -l || echo 0);
            n_high_quality=$(grep -q -P "High-quality\tHigh-quality" $checkv && grep -P "High-quality\tHigh-quality" $checkv | wc -l || echo 0);
            n_high_quality_proviruses=$(grep -q -E "(_fragment|\|provirus)" $checkv | grep -q -P "High-quality\tHigh-quality" && grep -P "High-quality\tHigh-quality" | wc -l || echo 0);
            echo "$ass;$n_vibrant;$n_cobra;$n_cobra_proviruses;$n_complete;$n_complete_proviruses;$n_high_quality;$n_high_quality_proviruses"
        done >> {output}
        """

rule statistics_final_viral_dataset:
    output: os.path.join(RESULTS_DIR, "statistics", "statistics_final_viral_dataset.csv")
    input:
        viruses_checkv1 = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "all_HQ_viruses.fna"),
        proviruses_checkv1 = os.path.join(RESULTS_DIR, "viruses", "first_checkv", "all_HQ_proviruses.fna"),
        viruses_checkv2 = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "all_HQ_viruses.fna"),
        proviruses_checkv2 = os.path.join(RESULTS_DIR, "viruses", "second_checkv", "all_HQ_proviruses.fna"),
        HQ_viruses_with_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "HQ_viruses_with_duplicates.fna"),
        HQ_viruses_no_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "HQ_viruses_no_duplicates.fna"),
        all_viruses_with_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "all_viruses_with_duplicates.fna"),
        all_viruses_no_duplicates = os.path.join(RESULTS_DIR, "viruses", "final_set_virus", "all_viruses_no_duplicates.fna"),
    message: "Computing statistics of final viral dataset"
    shell:
        """
        quality1=$(dirname {input.viruses_checkv1})
        quality2=$(dirname {input.viruses_checkv2})
        echo "#input_checkv1;#complete_checkv1;#all_HQ_checkv1;#input_checkv2;#complete_checkv2;#all_HQ_checkv2;#HQ_checkv12;#HQ_nodup_checkv12;#all_checkv12;#all_nodup_checkv12" > {output}

        ncontigs1=$(find $quality1 -name "quality_summary.tsv" -exec tail -n +2 {{}} \; | wc -l);
        complete1=$(find $quality1 -name "quality_summary.tsv" -exec grep -c "Complete" {{}} \; | awk '{{sum+=$1}} END {{print sum}}');
        HQ1=$(find $quality1 -name "quality_summary.tsv" -exec grep -c "High-quality" {{}} \; | awk '{{sum+=$1}} END {{print sum}}');

        ncontigs2=$(find $quality2 -name "quality_summary.tsv" -exec tail -n +2 {{}} \; | wc -l);
        complete2=$(find $quality2 -name "quality_summary.tsv" -exec grep -c "Complete" {{}} \; | awk '{{sum+=$1}} END {{print sum}}');
        HQ2=$(find $quality2 -name "quality_summary.tsv" -exec grep -c "High-quality" {{}} \; | awk '{{sum+=$1}} END {{print sum}}');

        HQ_total12=$(grep ">" {input.HQ_viruses_with_duplicates} | wc -l);
        HQ_total12_dedup=$(grep ">" {input.HQ_viruses_no_duplicates} | wc -l);
        all_total12=$(grep ">" {input.all_viruses_with_duplicates} | wc -l);
        all_total12_dedup=$(grep ">" {input.all_viruses_no_duplicates} | wc -l);
        echo "$ncontigs1;$complete1;$HQ1;$ncontigs2;$complete2;$HQ2;$HQ_total12;$HQ_total12_dedup;$all_total12;$all_total12_dedup" >> {output}
        """
