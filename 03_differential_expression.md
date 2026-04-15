# Module 03 — Differential Expression Analysis with DESeq2

## Overview

With a raw count matrix in hand, you are ready to identify genes that are significantly differentially expressed between wild type and the *gh31* deletion mutant. In this module you will:

1. Load count data into R and construct a DESeq2 object
2. Run the DESeq2 statistical workflow
3. Visualize results with a volcano plot and MA plot
4. Identify and interpret significantly differentially expressed genes
5. Contextualize findings within the biology of GH31

---

## 3.1 — Prepare the DESeq2 R Script

All DESeq2 analysis runs inside the Bioconductor container. Create the R script first, then submit it as a SLURM job.

```bash
mkdir -p ${TUTORIAL}/dge
cat > ${TUTORIAL}/scripts/deseq2_peltaster.R << 'REOF'
# ============================================================
# DESeq2 Differential Expression Analysis
# P. fructicola: WT vs gh31del
# HCS 7004 — Genome Analytics
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

# ----------------------------------------------------------
# 1. Load count matrix
# ----------------------------------------------------------
counts_file <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/counts/peltaster_counts.txt"

raw <- read.table(counts_file, header = TRUE, skip = 1, sep = "\t",
                  row.names = 1, check.names = FALSE)

# featureCounts columns 2-6 are metadata; counts start at column 6
# Column names include full path — clean them up
count_mat <- raw[, 6:ncol(raw)]
colnames(count_mat) <- sub(".*/", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.sorted\\.bam$", "", colnames(count_mat))

cat("Count matrix dimensions:", dim(count_mat), "\n")
cat("Sample names:", colnames(count_mat), "\n")
cat("Total read counts per sample:\n")
print(colSums(count_mat))

# ----------------------------------------------------------
# 2. Build sample metadata (colData)
# ----------------------------------------------------------
sample_info <- data.frame(
  row.names = colnames(count_mat),
  condition = factor(c("WT", "WT", "WT",
                       "gh31del", "gh31del", "gh31del"),
                     levels = c("WT", "gh31del"))
)
cat("\nSample metadata:\n")
print(sample_info)

# ----------------------------------------------------------
# 3. Create DESeqDataSet and pre-filter
# ----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ condition
)

# Remove genes with very low counts (< 10 total reads)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("\nGenes after low-count filtering:", nrow(dds), "\n")

# ----------------------------------------------------------
# 4. Run DESeq2
# ----------------------------------------------------------
dds <- DESeq(dds)

# ----------------------------------------------------------
# 5. Extract results: WT vs gh31del
# (positive log2FC = higher in WT)
# ----------------------------------------------------------
res <- results(dds,
               contrast = c("condition", "WT", "gh31del"),
               alpha = 0.05)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

cat("\nSummary of DESeq2 results:\n")
summary(res)

outdir <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/dge"

write.csv(res_df,
          file.path(outdir, "peltaster_deseq2_results.csv"),
          row.names = FALSE)

# Significant DEGs: padj < 0.05 and |log2FC| >= 1
sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

cat("\nSignificant DEGs (padj < 0.05, |log2FC| >= 1):", nrow(sig), "\n")
cat("Up in WT:", sum(sig$log2FoldChange > 0), "\n")
cat("Down in WT:", sum(sig$log2FoldChange < 0), "\n")
cat("\nTop 20 significant DEGs:\n")
print(head(sig, 20))

write.csv(sig,
          file.path(outdir, "peltaster_sig_DEGs.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 6. Variance-stabilizing transformation for visualization
# ----------------------------------------------------------
vst_data <- vst(dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "PCA of VST-normalized counts",
       subtitle = "P. fructicola WT vs gh31del",
       x = paste0("PC1: ", pca_var[1], "% variance"),
       y = paste0("PC2: ", pca_var[2], "% variance"),
       color = "Condition") +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "pca_plot.pdf"), pca_plot, width = 7, height = 5)

# ----------------------------------------------------------
# 7. Volcano plot
# ----------------------------------------------------------
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange >= 1]  <- "Up in WT"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down in WT"

# Highlight specific genes of interest
genes_of_interest <- c("t1.AMS68_008039",   # GH31 itself
                       "t1.AMS68_000995")   # check if it appears

volcano <- ggplot(res_df %>% filter(!is.na(padj)),
                  aes(x = log2FoldChange, y = -log10(padj),
                      color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Not significant" = "grey60",
                                "Up in WT"        = "#D73027",
                                "Down in WT"      = "#4575B4")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  # Annotate top gene
  geom_point(data = res_df %>% filter(gene_id %in% genes_of_interest),
             size = 3, shape = 21, color = "black", fill = "gold") +
  ggrepel::geom_label_repel(
    data = res_df %>% filter(gene_id %in% genes_of_interest),
    aes(label = gene_id), size = 3, fill = "white", max.overlaps = 10) +
  labs(title = "Volcano Plot: WT vs gh31del",
       subtitle = "P. fructicola differential expression",
       x = expression(log[2]~"Fold Change (WT / gh31del)"),
       y = expression(-log[10]~"(adjusted p-value)"),
       color = "") +
  theme_bw(base_size = 12) +
  xlim(-6, 6)

ggsave(file.path(outdir, "volcano_plot.pdf"), volcano, width = 8, height = 6)

# ----------------------------------------------------------
# 8. MA plot
# ----------------------------------------------------------
pdf(file.path(outdir, "ma_plot.pdf"), width = 7, height = 5)
plotMA(res, main = "MA Plot: WT vs gh31del", alpha = 0.05, ylim = c(-6, 6))
dev.off()

# ----------------------------------------------------------
# 9. Sample distance heatmap
# ----------------------------------------------------------
library(pheatmap)
sampleDists <- dist(t(assay(vst_data)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf(file.path(outdir, "sample_distance_heatmap.pdf"), width = 6, height = 5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample Distances (VST)",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255))
dev.off()

cat("\nAll output files written to:", outdir, "\n")
cat("DESeq2 analysis complete.\n")
REOF
```

Note: The PCA plot uses `ggrepel` for label placement. If the container does not include `ggrepel` or `pheatmap`, install them locally inside the container invocation (see below).

---

## 3.2 — Run DESeq2 on OSC

```bash
cat > ${TUTORIAL}/03_deseq2.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=deseq2_peltaster
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/deseq2_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/deseq2.sif

echo "[$(date)] Running DESeq2..."

apptainer exec --bind ${TUTORIAL} ${SIF} \
    Rscript ${TUTORIAL}/scripts/deseq2_peltaster.R

echo "[$(date)] DESeq2 complete."
ls -lh ${TUTORIAL}/dge/
EOF

sbatch ${TUTORIAL}/03_deseq2.sh
```

Monitor the job log in real time while it runs:

```bash
tail -f ${TUTORIAL}/logs/deseq2_*.log
```

---

## 3.3 — Explore the Results

After the job completes, download the output files:

```bash
# On your LOCAL machine:
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/dge/*.pdf ~/Desktop/
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/dge/*.csv ~/Desktop/
```

Open `pca_plot.pdf`, `volcano_plot.pdf`, and `ma_plot.pdf`.

Then on OSC, examine the significant DEG list:

```bash
# How many significant DEGs?
awk -F',' 'NR>1 && $7!="NA" && $7<0.05 && ($6>1 || $6<-1) {print}' \
    ${TUTORIAL}/dge/peltaster_deseq2_results.csv | wc -l

# Top upregulated genes in WT
awk -F',' 'NR>1 && $7!="NA" && $7<0.05 && $6>0 {print}' \
    ${TUTORIAL}/dge/peltaster_deseq2_results.csv \
    | sort -t',' -k6 -rn | head -10

# Top downregulated genes in WT
awk -F',' 'NR>1 && $7!="NA" && $7<0.05 && $6<0 {print}' \
    ${TUTORIAL}/dge/peltaster_deseq2_results.csv \
    | sort -t',' -k6 -n | head -10
```

---

## 3.4 — Interpretation Questions

> **Q16:** Look at the PCA plot. Do the three WT replicates cluster together? Do the three *gh31del* replicates cluster together? What does good replicate clustering indicate about the quality of the experiment?

> **Q17:** How many genes are significantly upregulated in WT (padj < 0.05, log2FC > 1) and how many are significantly downregulated in WT (padj < 0.05, log2FC < −1)?

> **Q18:** Identify the most significantly upregulated gene in WT. Record its gene ID, log2 fold change, and adjusted p-value. Based on its gene ID and the annotation, what does this gene encode? (Hint: you may need to look it up in the GFF3 file or NCBI.)

```bash
# Look up the top upregulated gene in the annotation
TOPGENE="<insert top upregulated gene ID>"
grep "${TOPGENE}" ${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3 | head -5
```

> **Q19:** Look for gene `AMS68_008039` (GH31 itself) in your DESeq2 results. Is it present? If not, why might the GH31 gene be absent from or have unusual behavior in the differential expression output when comparing WT to a *gh31* deletion strain?

> **Q20:** GH31 is a glycoside hydrolase involved in cleaving α-glucosides. Given that GH31 deletion causes upregulation of certain genes and downregulation of others, propose a biological model explaining how a secreted enzyme could influence transcriptional regulation. What kind of signaling mechanism might link extracellular carbohydrate hydrolysis to changes in nuclear gene expression?

> **Q21:** The DESeq2 results likely contain some genes with moderate fold changes (|log2FC| ~1–2) that are statistically significant. Are all statistically significant DEGs biologically meaningful? How would you distinguish true regulatory targets of GH31 from statistical noise in a real experiment?

---

## 3.5 — Dispersion and Model Diagnostics

Review the DESeq2 dispersion estimates to assess whether the negative binomial model fits your data well:

```bash
cat >> ${TUTORIAL}/scripts/deseq2_peltaster.R << 'REOF'
# Dispersion plot
pdf(file.path(outdir, "dispersion_plot.pdf"), width = 7, height = 5)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()
REOF
```

> **Q22:** In the dispersion plot, do the final (shrunken) dispersions follow the fitted trend line reasonably well? What would extremely scattered dispersions suggest about the quality of the replicates or the appropriateness of the model?

---

## Summary

By the end of this module you should have:

- [x] DESeq2 results CSV with all genes tested
- [x] Significant DEGs CSV (padj < 0.05, |log2FC| ≥ 1)
- [x] Volcano plot, MA plot, PCA plot, sample distance heatmap
- [x] A biological hypothesis about GH31's regulatory role

**Proceed to [Module 04 — ATAC-seq Processing](04_atacseq_processing.md)**
