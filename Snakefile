import os
import glob

FASTQ_DIR = "data/samples"

FASTQ_FILES = glob.glob(os.path.join(FASTQ_DIR, "*.fastq"))

SAMPLES = [os.path.splitext(os.path.basename(f))[0] for f in FASTQ_FILES]

rule all:
    input:
        expand("results/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        "results/multiqc/multiqc_report.html"

rule fastqc:
    input:
        lambda wildcards: os.path.join(FASTQ_DIR, f"{wildcards.sample}.fastq")
    output:
        html = "results/fastqc/{sample}_fastqc.html",
        zip  = "results/fastqc/{sample}_fastqc.zip"
    conda:
        "envs/environment.yaml"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc {input} --outdir=results/fastqc > {log} 2>&1"
        
rule multiqc:
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "envs/environment.yaml"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc results/fastqc --outdir results/multiqc > {log} 2>&1"

