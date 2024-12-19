# Mapping reads on assemblies with minimap2 for COBRA
rule minimap2_for_cobra:
    output: os.path.join(RESULTS_DIR, "viruses", "minimap2", "{assembly}_sorted_mapping_for_cobra.bam")
    input:
        assembly = os.path.join(RESULTS_DIR, "megahit", "{assembly}_megahit", "{assembly}.renamed.contigs.fa"),
        reads1 = os.path.join(RESULTS_DIR, "fastp", "{assembly}_R1_fastp.fastq.gz"),
        reads2 = os.path.join(RESULTS_DIR, "fastp", "{assembly}_R2_fastp.fastq.gz"),
    params:
        reference_bases_dealt_at_once = "16G"
    conda: os.path.join(ENV_DIR, "cobra.yaml")
    threads: config['minimap2']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_mapping_for_cobra.log")
    message: "Running minimap2 to prepare COBRA"
    shell:
        """(date && minimap2 -ax sr -I {params.reference_bases_dealt_at_once} -t {threads} {input.assembly} {input.reads1} {input.reads2} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output} && date) &> {log}"""

# Calculate coverages from bam files for COBRA
rule metabat2_summarize:
    output: os.path.join(RESULTS_DIR, "viruses", "metabat2", "{assembly}_contig_depths.txt")
    input: os.path.join(RESULTS_DIR, "viruses", "minimap2", "{assembly}_sorted_mapping_for_cobra.bam")
    conda: os.path.join(ENV_DIR, "cobra.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_metabat2_summarize.log")
    message: "Running metabat2 to prepare COBRA"
    shell:
        """(date && jgi_summarize_bam_contig_depths --outputDepth {output} {input} && date) &> {log}"""

rule metabat2_conversion:
    output: os.path.join(RESULTS_DIR, "viruses", "metabat2", "{assembly}_contig_depths_converted.txt")
    input: os.path.join(RESULTS_DIR, "viruses", "metabat2", "{assembly}_contig_depths.txt")
    conda: os.path.join(ENV_DIR, "cobra.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_metabat2_conversion.log")
    message: "Converting metabat2 coverages to prepare COBRA"
    shell:
        """(date && python scripts/coverage.transfer.py -i {input} -o {output} && date) &> {log}"""

# Launch COBRA using all prepared inputs
rule cobra:
    output: os.path.join(RESULTS_DIR, "viruses", "cobra", "{assembly}_cobra", "COBRA_joining_summary.txt"),
    input:
        all_contigs = os.path.join(RESULTS_DIR, "megahit", "{assembly}_megahit", "{assembly}.renamed.contigs.fa"),
        viral_contigs = os.path.join(RESULTS_DIR, "viruses", "prep_cobra", "{assembly}_for_COBRA.fna"),
        all_mapping = os.path.join(RESULTS_DIR, "viruses", "minimap2", "{assembly}_sorted_mapping_for_cobra.bam"),
        all_coverage = os.path.join(RESULTS_DIR, "viruses", "metabat2", "{assembly}_contig_depths_converted.txt"),
    conda: os.path.join(ENV_DIR, "cobra.yaml")
    threads: config['cobra']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_cobra.log")
    message: "Running COBRA"
    shell:
        """(date && mkdir -p $(dirname $(dirname {output})) && rm -r $(dirname {output}) &&
        cobra-meta -t {threads} -f {input.all_contigs} -q {input.viral_contigs} -c {input.all_coverage} -m {input.all_mapping} -a megahit -mink 21 -maxk 255 -o $(dirname {output}) && date) &> {log}"""

# Statistics of viruses obtained post COBRA
rule statistics_cobra:
    output: os.path.join(RESULTS_DIR, "statistics", "statistics_COBRA.csv")
    input: expand(os.path.join(RESULTS_DIR, "viruses", "cobra", "{assembly}_cobra", "COBRA_joining_summary.txt"), assembly=SAMPLES_VIR),
    message: "Computing statistics post COBRA extension"
    shell:
        """
        echo "sample;#self_circular;#extended_circular;#extended_circular_unique;#extended_partial;#extended_partial_unique;#extended_failed;#orphan_end" > {output}
        for cobra_file in {input}; do
            cobra_dir=$(dirname $cobra_file);
            self_circular=$(grep -q -Po '(?<=# Category i   - Self_circular: )[0-9]+' $(echo $cobra_dir)/log && grep -Po '(?<=# Category i   - Self_circular: )[0-9]+' $(echo $cobra_dir)/log || echo 0);
            extended_circular=$(grep -q -Po '(?<=# Category ii  - Extended_circular: )[0-9]+' $(echo $cobra_dir)/log && grep -Po '(?<=# Category ii  - Extended_circular: )[0-9]+' $(echo $cobra_dir)/log || echo 0);
            extended_circular_unique=$(grep -q ">" $(echo $cobra_dir)/COBRA_category_ii-a_extended_circular_unique.fasta && grep ">" $(echo $cobra_dir)/COBRA_category_ii-a_extended_circular_unique.fasta | wc -l || echo 0);
            extended_partial=$(grep -q -Po '(?<=# Category ii  - Extended_partial: )[0-9]+' $(echo $cobra_dir)/log && grep -Po '(?<=# Category ii  - Extended_partial: )[0-9]+' $(echo $cobra_dir)/log || echo 0);
            extended_partial_unique=$(grep -q ">" $(echo $cobra_dir)/COBRA_category_ii-b_extended_partial_unique.fasta && grep ">" $(echo $cobra_dir)/COBRA_category_ii-b_extended_partial_unique.fasta | wc -l || echo 0);
            extended_failed=$(grep -q -Po '(?<=# Category ii  - Extended_failed \(due to COBRA rules\): )[0-9]+' $(echo $cobra_dir)/log && grep -Po '(?<=# Category ii  - Extended_failed \(due to COBRA rules\): )[0-9]+' $(echo $cobra_dir)/log || echo 0);
            orphan_end=$(grep -q -Po '(?<=# Category iii - Orphan end: )[0-9]+' $(echo $cobra_dir)/log && grep -Po '(?<=# Category iii - Orphan end: )[0-9]+' $(echo $cobra_dir)/log || echo 0);
            echo "$(basename $cobra_dir);$self_circular;$extended_circular;$extended_circular_unique;$extended_partial;$extended_partial_unique;$extended_failed;$orphan_end"
        done >> {output}
        """
