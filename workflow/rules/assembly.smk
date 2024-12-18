# Megahit rule
rule assembly_all:
    input: expand(os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.contigs.fa"), assembly=SAMPLES)

rule megahit:
    output: os.path.join(RESULTS_DIR, "megahit", "{sample}_megahit", "{sample}.contigs.fa")
    input: 
        reads1 = os.path.join(READS_DIR, "{sample}_R1.fastq.gz"),,
        reads2 = os.path.join(READS_DIR, "{sample}_R2.fastq.gz")
    conda: os.path.join(ENV_DIR, "assembly.yaml")
    threads: config['megahit']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_megahit.log")
    message: "Running MEGAHIT"
    shell:
        """(date && megahit --presets meta-large --min-contig-len 300 -t {threads} -f -o $(dirname {output}) --out-prefix {wildcards.sample} -1 {input.reads1} -2 {input.reads2} && date) &> {log}"""
