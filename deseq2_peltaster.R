#!/usr/bin/env Rscript
# ============================================================
# DESeq2 Differential Expression Analysis
# Peltaster fructicola: WT vs gh31del
# HCS 7004 — Genome Analytics
#
# Usage: Rscript deseq2_peltaster.R
# Expected input: featureCounts output at path specified below
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# --- User-configurable paths ---
COUNTS_FILE <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/counts/peltaster_counts.txt"
OUTDIR      <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/dge"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 1. Load count matrix
# ----------------------------------------------------------
raw <- read.table(COUNTS_FILE, header = TRUE, skip = 1, sep = "\t",
                  row.names = 1, check.names = FALSE)

# featureCounts: first 5 columns are Geneid, Chr, Start, End, Strand, Length
# Count columns start at column 6
count_mat <- raw[, 6:ncol(raw)]
colnames(count_mat) <- sub(".*/", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.sorted\\.bam$", "", colnames(count_mat))

cat("=== Count Matrix ===\n")
cat("Dimensions:", dim(count_mat), "\n")
cat("Sample names:", colnames(count_mat), "\n")
cat("\nTotal counts per sample:\n")
print(colSums(count_mat))

# ----------------------------------------------------------
# 2. Sample metadata
# ----------------------------------------------------------
sample_info <- data.frame(
  row.names = colnames(count_mat),
  condition = factor(
    c("WT", "WT", "WT", "gh31del", "gh31del", "gh31del"),
    levels = c("WT", "gh31del")
  )
)
cat("\nSample metadata:\n")
print(sample_info)

# ----------------------------------------------------------
# 3. DESeqDataSet construction and pre-filtering
# ----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ condition
)

# Pre-filter: keep genes with >= 10 total reads
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("\nGenes retained after filtering:", nrow(dds), "\n")

# ----------------------------------------------------------
# 4. Run DESeq2
# ----------------------------------------------------------
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

# Size factors (library size normalization)
cat("\nSize factors:\n")
print(sizeFactors(dds))

# ----------------------------------------------------------
# 5. Extract results
# Positive log2FC = higher expression in WT
# ----------------------------------------------------------
res <- results(dds,
               contrast = c("condition", "WT", "gh31del"),
               alpha    = 0.05)

cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

write.csv(res_df,
          file.path(OUTDIR, "peltaster_deseq2_results.csv"),
          row.names = FALSE)

# Significant DEGs
sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

cat("\n=== Significant DEGs (padj < 0.05, |log2FC| >= 1) ===\n")
cat("Total:", nrow(sig), "\n")
cat("Upregulated in WT (log2FC > 1):", sum(sig$log2FoldChange > 0), "\n")
cat("Downregulated in WT (log2FC < -1):", sum(sig$log2FoldChange < 0), "\n")

cat("\nTop 20 DEGs:\n")
print(head(sig[, c("gene_id","baseMean","log2FoldChange","padj")], 20))

write.csv(sig,
          file.path(OUTDIR, "peltaster_sig_DEGs.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 6. Variance-stabilizing transformation
# ----------------------------------------------------------
vst_data <- vst(dds, blind = FALSE)

# ----------------------------------------------------------
# 7. PCA plot
# ----------------------------------------------------------
pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data,
                   aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 5) +
  geom_text_repel(size = 3.5) +
  scale_color_manual(values = c("WT" = "#2196F3", "gh31del" = "#E53935")) +
  labs(title = "PCA — VST normalized counts",
       subtitle = "P. fructicola WT vs. gh31 deletion",
       x = paste0("PC1: ", pca_var[1], "% variance"),
       y = paste0("PC2: ", pca_var[2], "% variance"),
       color = "Condition") +
  theme_bw(base_size = 13)

ggsave(file.path(OUTDIR, "pca_plot.pdf"), pca_plot, width = 7, height = 5)
cat("\nPCA plot saved.\n")

# ----------------------------------------------------------
# 8. Volcano plot
# ----------------------------------------------------------
res_df <- res_df %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up in WT",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down in WT",
      TRUE ~ "Not significant"
    )
  )

# Genes to label
genes_label <- res_df %>%
  filter(grepl("AMS68_008039|AMS68_000995", gene_id))

volcano <- ggplot(
  res_df %>% filter(!is.na(padj)),
  aes(x = log2FoldChange, y = -log10(padj), color = significance)
) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c(
    "Not significant" = "grey60",
    "Up in WT"        = "#C62828",
    "Down in WT"      = "#1565C0"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", linewidth = 0.4) +
  geom_point(data = genes_label,
             size = 3.5, shape = 21, color = "black", fill = "#FFD600") +
  geom_label_repel(data = genes_label,
                   aes(label = gene_id),
                   size = 3, fill = "white", max.overlaps = 20,
                   box.padding = 0.5) +
  labs(
    title    = "Volcano Plot: WT vs. gh31del",
    subtitle = "Peltaster fructicola differential expression",
    x        = expression(log[2]~"Fold Change (WT / gh31del)"),
    y        = expression(-log[10]~"adjusted p-value"),
    color    = ""
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  theme_bw(base_size = 13)

ggsave(file.path(OUTDIR, "volcano_plot.pdf"), volcano, width = 8, height = 6)
cat("Volcano plot saved.\n")

# ----------------------------------------------------------
# 9. MA plot
# ----------------------------------------------------------
pdf(file.path(OUTDIR, "ma_plot.pdf"), width = 7, height = 5)
plotMA(res, main = "MA Plot: WT vs gh31del", alpha = 0.05, ylim = c(-6, 6))
dev.off()
cat("MA plot saved.\n")

# ----------------------------------------------------------
# 10. Dispersion estimates
# ----------------------------------------------------------
pdf(file.path(OUTDIR, "dispersion_plot.pdf"), width = 7, height = 5)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()
cat("Dispersion plot saved.\n")

# ----------------------------------------------------------
# 11. Sample distance heatmap
# ----------------------------------------------------------
sampleDists <- dist(t(assay(vst_data)))
sdm         <- as.matrix(sampleDists)

pdf(file.path(OUTDIR, "sample_distance_heatmap.pdf"), width = 6, height = 5)
pheatmap(sdm,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main  = "Sample-to-Sample Distances (VST)",
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()
cat("Sample distance heatmap saved.\n")

cat("\n=== Analysis complete ===\n")
cat("Output files in:", OUTDIR, "\n")
