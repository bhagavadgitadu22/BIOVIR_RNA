# Megahit rule
rule assembly_all:
    input: expand(os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.renamed.contigs.fa"), sample=SAMPLES)

rule megahit:
    output: os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.contigs.fa")
    input: 
        reads1 = os.path.join(RESULTS_DIR, "fastp", "{sample}_R1_fastp.fastq.gz"),
        reads2 = os.path.join(RESULTS_DIR, "fastp", "{sample}_R2_fastp.fastq.gz")
    conda: os.path.join(ENV_DIR, "assembly.yaml")
    threads: config['megahit']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_megahit.log")
    message: "Running MEGAHIT"
    shell:
        """(date && megahit --presets meta-sensitive --min-contig-len 500 -t {threads} -f -o $(dirname {output}) --out-prefix {wildcards.sample} -1 {input.reads1} -2 {input.reads2} && date) &> {log}"""

rule rename_contigs:
    output: os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.renamed.contigs.fa")
    input: os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.contigs.fa")
    message: "Running renaming of MEGAHIT"
    shell:
        """sed -E 's/>k.*?_/>{wildcards.sample}_/g' {input} | sed -E 's/ flag=.*$//g' > {output}"""

# Report
rule multiqc_assembly:
    output: os.path.join(RESULTS_DIR, "multiqc_assembly/multiqc_report.html")
    input: expand(os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.contigs.fa"), sample=SAMPLES)
    params: RESULTS_DIR
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "multiqc_assembly.log")
    message: "Running MultiQC"
    shell:
        """(date && multiqc {params}/megahit/*_megahit -o $(dirname {output}) && date) &> {log}"""
