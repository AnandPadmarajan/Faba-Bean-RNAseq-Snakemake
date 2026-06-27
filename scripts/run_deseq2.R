#!/usr/bin/env Rscript
# =============================================================================
# run_deseq2.R
#
# Runs DESeq2 differential expression on a featureCounts count matrix.
# Design formula, contrast, and significance threshold are all read from
# config.yaml so this script is reusable across datasets without editing code.
#
# Usage (called from Snakemake):
#   Rscript run_deseq2.R <counts_matrix.tsv> <sample_table.tsv> <design_formula> \
#       <contrast_factor> <contrast_numerator> <contrast_denominator> \
#       <padj_threshold> <output_prefix>
# =============================================================================

suppressMessages({
  library(DESeq2)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
counts_path     <- args[1]
sample_tbl_path <- args[2]
design_formula  <- as.formula(args[3])
contrast_factor <- args[4]
contrast_num    <- args[5]
contrast_denom  <- args[6]
padj_threshold  <- as.numeric(args[7])
out_prefix      <- args[8]

# --- Load data ---
counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
sample_table <- read.delim(sample_tbl_path, row.names = 1)

# Sanity check: sample order must match between counts and sample table.
# This is a common, silent source of bugs in RNA-seq pipelines.
stopifnot(all(colnames(counts) == rownames(sample_table)))

# --- Build DESeq2 dataset ---
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_table,
  design    = design_formula
)

# Pre-filter very low count genes (standard DESeq2 recommendation)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(contrast_factor, contrast_num, contrast_denom)
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

sig_genes <- res_df %>% filter(!is.na(padj), padj < padj_threshold)

write.table(res_df, paste0(out_prefix, "_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_genes, paste0(out_prefix, "_significant_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Variance-stabilized counts, used downstream by WGCNA.
# vst() needs a reasonably large number of genes to fit a dispersion trend;
# on very small gene sets (e.g. toy/test data) it can error out, so we fall
# back to the slower but more robust varianceStabilizingTransformation().
vsd <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("vst() failed (likely too few genes for trend fitting); ",
            "falling back to varianceStabilizingTransformation()")
    varianceStabilizingTransformation(dds, blind = TRUE)
  }
)
vsd_mat <- assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")
write.table(vsd_mat, paste0(out_prefix, "_vst_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(
  "DESeq2 complete: %d genes tested, %d significant at padj < %s\n",
  nrow(res_df), nrow(sig_genes), padj_threshold
))
