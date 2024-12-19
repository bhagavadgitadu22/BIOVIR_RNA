# Use fastp to remove bad reads
rule fastp:
    output:
        out1 = os.path.join(RESULTS_DIR, "fastp", "{sample}_R1_fastp.fastq.gz"),
        out2 = os.path.join(RESULTS_DIR, "fastp", "{sample}_R2_fastp.fastq.gz"),
        json = os.path.join(RESULTS_DIR, "fastp", "{sample}_fastp.json"),
        html = os.path.join(RESULTS_DIR, "fastp", "{sample}_fastp.html"),
    input:
        read1 = os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        read2 = os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"),
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    threads: config['fastp']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_fastp.log")
    message: "Running fastp" 
    shell:
        """(date && fastp -t {threads} -i {input.read1} -I {input.read2} --cut_mean_quality 30 --json {output.json} --html {output.html} -o {output.out1} -O {output.out2} && date) &> {log}"""

# Combine all the samples. We still want to keep the read pairs separate
rule multiqc_preprocessing:
    output: os.path.join(RESULTS_DIR, "multiqc_preprocessing/multiqc_report.html")
    input: expand(os.path.join(RESULTS_DIR, "fastp", "{sample}_fastp.html"), sample=SAMPLES)
    params: RESULTS_DIR
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "multiqc_preprocessing.log")
    message: "Running MultiQC"
    shell:
        """(date && multiqc {params}/fastp -o $(dirname {output}) && date) &> {log}"""
