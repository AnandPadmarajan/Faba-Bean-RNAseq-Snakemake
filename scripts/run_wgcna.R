#!/usr/bin/env Rscript
# =============================================================================
# run_wgcna.R
#
# Builds a weighted gene co-expression network from variance-stabilized counts
# and identifies co-expression modules.
#
# NOTE on module naming: WGCNA's blockwiseModules() can return purely numeric
# module labels (e.g. "1", "2") instead of color names (e.g. "turquoise",
# "violet") depending on settings. Numeric labels are easy to mix up with row
# indices and module sizes downstream, so this script always explicitly maps
# numeric labels to color names via labels2colors() -- this mirrors a fix
# (fix_module_names.R) developed during the original coursework project after
# a numeric/color mismatch bug surfaced in downstream analysis.
#
# Usage (called from Snakemake):
#   Rscript run_wgcna.R <vst_counts.tsv> <min_module_size> <merge_cut_height> \
#       <soft_power> <output_prefix>
# =============================================================================

suppressMessages({
  library(WGCNA)
  library(tidyverse)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
vst_path         <- args[1]
min_module_size  <- as.numeric(args[2])
merge_cut_height <- as.numeric(args[3])
soft_power       <- as.numeric(args[4])
out_prefix       <- args[5]

# --- Load and orient data: WGCNA expects samples as rows, genes as columns ---
expr <- read.delim(vst_path, row.names = 1, check.names = FALSE)
expr_t <- t(expr)

# --- Soft-thresholding power check ---
# soft_power is read from config as a starting point. In a real run, inspect
# the scale-free topology fit plot (written below) and update config.yaml's
# wgcna.soft_power if the chosen value doesn't show a clear elbow.
sft <- pickSoftThreshold(expr_t, powerVector = c(1:20), verbose = 0)
write.table(sft$fitIndices, paste0(out_prefix, "_soft_threshold_fit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Network construction + module detection ---
net <- blockwiseModules(
  expr_t,
  power = soft_power,
  minModuleSize = min_module_size,
  mergeCutHeight = merge_cut_height,
  TOMType = "unsigned",
  numericLabels = TRUE,   # deliberately TRUE: we map to colors explicitly below
  saveTOMs = FALSE,
  verbose = 0
)

# --- Explicit numeric -> color mapping (the fix_module_names.R logic) ---
module_colors <- labels2colors(net$colors)

module_assignments <- tibble(
  gene_id = colnames(expr_t),
  module_numeric = net$colors,
  module_color = module_colors
)

write.table(module_assignments, paste0(out_prefix, "_module_assignments.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

module_sizes <- module_assignments %>%
  count(module_color, name = "n_genes") %>%
  arrange(desc(n_genes))

write.table(module_sizes, paste0(out_prefix, "_module_sizes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(
  "WGCNA complete: %d genes assigned to %d modules (color-named, not numeric)\n",
  nrow(module_assignments), nrow(module_sizes)
))
