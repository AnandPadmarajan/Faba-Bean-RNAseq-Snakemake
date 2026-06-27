"""
Snakefile: faba-bean-rnaseq-snakemake
======================================

Snakemake re-implementation of an RNA-seq + WGCNA differential expression
pipeline originally developed as coursework analyzing the transcriptomic
response of Vicia faba (faba bean) to high night temperature (HNT) stress
(MSc Agrobioinformatics, JLU Giessen).

Pipeline stages:
  1. QC            - FastQC + Trim Galore
  2. Alignment     - HISAT2 (chosen over Bowtie2 based on benchmarked mapping
                      rate: ~91% vs ~70% on test data in the original project)
  3. Quantification- featureCounts (subread >=2.1.1, avoids a known segfault
                      present in earlier subread releases)
  4. DE analysis   - DESeq2, with design formula and contrast set in config.yaml
  5. Co-expression - WGCNA, with explicit numeric->color module label mapping

Run with:
    snakemake --cores 4 --use-conda

See README.md for full setup instructions.
"""

import os
configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())
RESULTS = config["results_dir"]


def get_fastq(wildcards, read):
    return config["samples"][wildcards.sample][f"fastq_{read}"]


rule all:
    input:
        f"{RESULTS}/05_deseq2/deseq2_significant_genes.tsv",
        f"{RESULTS}/06_wgcna/wgcna_module_sizes.tsv",
        f"{RESULTS}/01_qc/multiqc_report.html"


# -----------------------------------------------------------------------------
# 1. Quality control
# -----------------------------------------------------------------------------
rule fastqc:
    input:
        r1=lambda wc: get_fastq(wc, 1),
        r2=lambda wc: get_fastq(wc, 2)
    output:
        html1=f"{RESULTS}/01_qc/{{sample}}_R1_fastqc.html",
        html2=f"{RESULTS}/01_qc/{{sample}}_R2_fastqc.html"
    conda:
        "envs/qc.yaml"
    params:
        outdir=f"{RESULTS}/01_qc"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} -o {params.outdir}
        """

rule multiqc:
    input:
        expand(f"{RESULTS}/01_qc/{{sample}}_R1_fastqc.html", sample=SAMPLES)
    output:
        f"{RESULTS}/01_qc/multiqc_report.html"
    conda:
        "envs/qc.yaml"
    params:
        qcdir=f"{RESULTS}/01_qc"
    shell:
        """
        # multiqc ships with fastqc in many bioconda setups; if not present,
        # add 'multiqc' to envs/qc.yaml.
        multiqc {params.qcdir} -o {params.qcdir} -n multiqc_report.html || \
            touch {output}
        """


# -----------------------------------------------------------------------------
# 2. Alignment (HISAT2)
# -----------------------------------------------------------------------------
rule hisat2_build:
    input:
        fasta=config["genome_fasta"]
    output:
        f"{config['hisat2_index_prefix']}.1.ht2"
    conda:
        "envs/hisat2.yaml"
    params:
        prefix=config["hisat2_index_prefix"]
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        hisat2-build {input.fasta} {params.prefix}
        """

rule hisat2_align:
    input:
        r1=lambda wc: get_fastq(wc, 1),
        r2=lambda wc: get_fastq(wc, 2),
        index=f"{config['hisat2_index_prefix']}.1.ht2"
    output:
        bam=f"{RESULTS}/02_aligned/{{sample}}.sorted.bam"
    conda:
        "envs/hisat2.yaml"
    threads: 4
    params:
        prefix=config["hisat2_index_prefix"]
    shell:
        """
        mkdir -p {RESULTS}/02_aligned
        hisat2 -p {threads} -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} | \
            samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# 3. Quantification (featureCounts)
# -----------------------------------------------------------------------------
rule feature_counts:
    input:
        bams=expand(f"{RESULTS}/02_aligned/{{sample}}.sorted.bam", sample=SAMPLES),
        gtf=config["genome_gtf"]
    output:
        counts=f"{RESULTS}/03_counts/raw_counts.tsv"
    conda:
        "envs/featurecounts.yaml"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS}/03_counts
        featureCounts -T {threads} -p --countReadPairs \
            -a {input.gtf} -o {output.counts} {input.bams}
        """

rule clean_counts_matrix:
    """
    featureCounts output has extra metadata columns and full BAM paths as
    column headers. This rule trims it down to a clean gene_id x sample
    matrix for DESeq2.
    """
    input:
        f"{RESULTS}/03_counts/raw_counts.tsv"
    output:
        f"{RESULTS}/03_counts/clean_counts.tsv"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t", skiprows=1)
        meta_cols = ["Chr", "Start", "End", "Strand", "Length"]
        df = df.drop(columns=meta_cols)
        df.columns = ["gene_id"] + SAMPLES
        df.to_csv(output[0], sep="\t", index=False)

rule sample_table:
    output:
        f"{RESULTS}/03_counts/sample_table.tsv"
    run:
        import pandas as pd
        rows = [
            {"sample": s, "condition": config["samples"][s]["condition"]}
            for s in SAMPLES
        ]
        pd.DataFrame(rows).set_index("sample").to_csv(output[0], sep="\t")


# -----------------------------------------------------------------------------
# 4. Differential expression (DESeq2)
# -----------------------------------------------------------------------------
rule deseq2:
    input:
        counts=f"{RESULTS}/03_counts/clean_counts.tsv",
        sample_table=f"{RESULTS}/03_counts/sample_table.tsv"
    output:
        all_genes=f"{RESULTS}/05_deseq2/deseq2_all_genes.tsv",
        sig_genes=f"{RESULTS}/05_deseq2/deseq2_significant_genes.tsv",
        vst=f"{RESULTS}/05_deseq2/deseq2_vst_counts.tsv"
    conda:
        "envs/deseq2.yaml"
    params:
        formula=config["design_formula"],
        contrast_factor=config["contrast"][0],
        contrast_num=config["contrast"][1],
        contrast_denom=config["contrast"][2],
        padj=config["padj_threshold"],
        out_prefix=f"{RESULTS}/05_deseq2/deseq2"
    shell:
        """
        mkdir -p {RESULTS}/05_deseq2
        Rscript scripts/run_deseq2.R \
            {input.counts} {input.sample_table} "{params.formula}" \
            {params.contrast_factor} {params.contrast_num} {params.contrast_denom} \
            {params.padj} {params.out_prefix}
        """


# -----------------------------------------------------------------------------
# 5. Co-expression network (WGCNA)
# -----------------------------------------------------------------------------
rule wgcna:
    input:
        vst=f"{RESULTS}/05_deseq2/deseq2_vst_counts.tsv"
    output:
        modules=f"{RESULTS}/06_wgcna/wgcna_module_assignments.tsv",
        sizes=f"{RESULTS}/06_wgcna/wgcna_module_sizes.tsv"
    conda:
        "envs/wgcna.yaml"
    params:
        min_size=config["wgcna"]["min_module_size"],
        merge_height=config["wgcna"]["merge_cut_height"],
        power=config["wgcna"]["soft_power"],
        out_prefix=f"{RESULTS}/06_wgcna/wgcna"
    shell:
        """
        mkdir -p {RESULTS}/06_wgcna
        Rscript scripts/run_wgcna.R \
            {input.vst} {params.min_size} {params.merge_height} {params.power} \
            {params.out_prefix}
        """
