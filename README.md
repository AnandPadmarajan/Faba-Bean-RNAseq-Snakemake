# Faba Bean RNA-seq + WGCNA Pipeline (Snakemake)

A reproducible Snakemake re-implementation of an RNA-seq differential
expression and co-expression network analysis pipeline, originally developed
as part of MSc Agrobioinformatics coursework (Systems Biology II, JLU
Giessen) analyzing the transcriptomic response of faba bean (*Vicia faba*)
to high night temperature (HNT) stress.

## Why this repository exists

The original analysis was run on a university-provided VM using coursework
data that is **unpublished and not mine to redistribute**. Rather than leave
the pipeline locked away on infrastructure I won't have permanent access to,
this repository extracts the pipeline *logic* — tool choices, parameters,
and a couple of debugging fixes discovered along the way — into a properly
engineered Snakemake workflow, and demonstrates it end-to-end on a small
**public** RNA-seq test dataset so that anyone can clone this repo and run it
without needing access to the original data.

To run the pipeline on your own RNA-seq data (including a real *Vicia faba*
dataset, if you have rights to one), see [Running on your own data](#running-on-your-own-data)
below — no code changes are required, only `config.yaml`.

## Pipeline overview

| Stage | Tool | Notes |
|---|---|---|
| QC | FastQC + MultiQC | Standard pre-alignment quality check |
| Alignment | HISAT2 | Chosen over Bowtie2 after benchmarking on test reads (~91% vs ~70% mapping rate in the original project) |
| Quantification | featureCounts (subread ≥2.1.1) | Earlier subread versions segfaulted on this workload; pinned to a fixed version |
| Differential expression | DESeq2 | Design formula and contrast are configurable, not hardcoded |
| Co-expression network | WGCNA | Module labels are explicitly mapped from numeric IDs to color names — see [Module naming fix](#module-naming-fix) |

```
Raw FASTQ → FastQC/MultiQC → HISAT2 alignment → featureCounts → DESeq2 → WGCNA
```

## Module naming fix

During the original project, WGCNA's `blockwiseModules()` was used with
numeric module labels, which were then easy to silently confuse with row
indices or module rank in downstream tables. `scripts/run_wgcna.R` always
explicitly converts numeric labels to WGCNA's standard color names via
`labels2colors()` before writing any output — this is baked into the
pipeline by default, not left as a manual post-processing step.

## Repository structure

```
.
├── Snakefile              # Pipeline rules
├── config.yaml            # Samples, paths, and parameters (edit this to use your own data)
├── envs/                  # One Conda environment per tool, for reproducibility
├── scripts/
│   ├── fetch_test_data.sh # Downloads the small public demo dataset
│   ├── run_deseq2.R
│   └── run_wgcna.R
├── test_data/             # Demo dataset lives here after running fetch_test_data.sh
├── results/               # Pipeline outputs (generated, not committed)
└── .github/workflows/     # CI: runs the full pipeline on every push
```

## Quick start

Requires [Conda or Mamba](https://github.com/conda-forge/miniforge) and
Snakemake (tested on WSL2/Ubuntu and native Linux).

```bash
# 1. Create an environment for Snakemake itself
mamba create -n snakemake -c conda-forge -c bioconda snakemake-minimal -y
conda activate snakemake

# 2. Clone this repo
git clone https://github.com/<your-username>/faba-bean-rnaseq-snakemake.git
cd faba-bean-rnaseq-snakemake

# 3. Fetch the small public demo dataset
bash scripts/fetch_test_data.sh
# (after this completes, double check the downloaded filenames against
#  config.yaml's `samples:` section -- see the note printed by the script)

# 4. Dry-run to see the planned jobs without executing anything
snakemake -n -p

# 5. Run the full pipeline (each rule creates its own Conda env on first run)
snakemake --cores 4 --use-conda
```

Outputs land in `results/`, including:
- `results/01_qc/multiqc_report.html` — aggregated QC report
- `results/05_deseq2/deseq2_significant_genes.tsv` — DEGs at the configured padj threshold
- `results/06_wgcna/wgcna_module_assignments.tsv` — gene → module (color name) mapping

## Running on your own data

Edit `config.yaml`:
1. Point `genome_fasta` / `genome_gtf` at your reference files
2. List your samples under `samples:`, with `condition` and FASTQ paths
3. Adjust `design_formula` / `contrast` to match your experiment
4. Re-run `snakemake --cores N --use-conda`

No script changes are needed — the design formula, contrast, and WGCNA
parameters are all read from `config.yaml`.

## Continuous integration

Every push runs the full pipeline end-to-end on the public test dataset via
GitHub Actions (see `.github/workflows/ci.yaml`), so the badge below
reflects whether the pipeline is currently reproducible from a clean clone —
not just whether it worked once on a laptop.

## Future improvements

- Containerize each rule (Docker/Singularity) as an alternative to per-rule
  Conda environments, for environments where Conda solve times are a
  bottleneck
- Add automated unit tests for the Python helper logic in the Snakefile
  (count matrix cleaning, sample table generation)
- Extend WGCNA module-trait correlation and GO enrichment as additional
  downstream rules

## Background

This pipeline's design choices (HISAT2 vs Bowtie2, the subread version pin,
explicit WGCNA module color mapping, and the `~condition`-only-style design
formula over an interaction model) were developed and validated on a real
faba bean HNT-stress RNA-seq dataset as part of MSc Agrobioinformatics
coursework. That original analysis is documented in a full scientific report
(not included here due to data restrictions); this repository represents
the reusable engineering work extracted from that project.
