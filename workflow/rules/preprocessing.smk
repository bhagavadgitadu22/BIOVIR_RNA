# Use cutadapt to remove adapters. Note that cutadapt works on the paired sequences
rule fastp:
    output:
        cut1 = os.path.join(RESULTS_DIR, "cutadapt", "{sample}_R1_cutadapt.fastq.gz"),
        cut2 = os.path.join(RESULTS_DIR, "cutadapt", "{sample}_R2_cutadapt.fastq.gz"),
        json = os.path.join(RESULTS_DIR, "cutadapt", "{sample}_cutadapt.json"),
    input:
        read1 = os.path.join(READS_DIR, "{sample}_R1.fastq.gz"),
        read2 = os.path.join(READS_DIR, "{sample}_R2.fastq.gz"),
    params:
        adapter_fwd = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adapter_rev = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        min_length = 50,
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    threads: config['cutadapt']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_cutadapt.log")
    message: "Running cutadapt" 
    shell:
        """(date && cutadapt -j {threads} -a {params.adapter_fwd} -A {params.adapter_rev} -m {params.min_length} --json {output.json} -o {output.cut1} -p {output.cut2} {input.read1} {input.read2} && date) &> {log}"""

# FastQC
rule fastqc:
    output: os.path.join(RESULTS_DIR, "fastqc", "{sample}_{read}_cutadapt_fastqc.html")
    input: os.path.join(READS_DIR, "{sample}_{read}.fastq.gz"),
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{read}_cutadapt_fastqc.log")
    message: "Running FastQC"
    shell:
        """(date && fastqc -o $(dirname {output}) {input} && date) &> {log}"""
        
# Combine all the samples. We still want to keep the read pairs separate
rule multiqc:
    output: os.path.join(RESULTS_DIR, "multiqc_preprocessing/multiqc_report.html")
    input: expand(os.path.join(RESULTS_DIR, "fastqc", "{sample}_{read}_cutadapt_fastqc.html"), sample=SAMPLES, read=["R1", "R2"])
    params: RESULTS_DIR
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "multiqc_preprocessing.log")
    message: "Running MultiQC"
    shell:
        """(date && multiqc {params}/cutadapt {params}/fastqc -o $(dirname {output}) && date) &> {log}"""
