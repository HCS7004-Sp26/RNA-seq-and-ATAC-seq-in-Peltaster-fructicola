# Module 03 — Differential Expression Analysis with DESeq2

## Overview

With a raw count matrix in hand, you are ready to identify genes that are significantly differentially expressed between wild type and the *gh31* deletion mutant. In this module you will:

1. Load count data into R and construct a DESeq2 object
2. Run the DESeq2 statistical workflow
3. Visualize results with a volcano plot and MA plot
4. Identify and interpret significantly differentially expressed genes
5. Contextualize findings within the biology of GH31

---

## 3.0 — R Environment on OSC

This module uses the R installation available as an OSC module rather than an Apptainer container. OSC provides R 4.5.2 via the GCC-compiled module stack:

```bash
module load gcc/12.3.0
module load R/4.5.2
```

Required packages are installed **once** to your personal R library and then reused in all subsequent jobs. You only need to run the installation step once per OSC account.

### 3.0.1 — Set Up Your Personal R Library

By default, R on OSC will install user packages to `~/R/x86_64-pc-linux-gnu-library/4.5/`. Confirm this path and install all required packages in a single interactive session:

```bash
# Start an interactive session (do not install packages on the login node)
sinteractive -A PAS3260 -t 00:30:00 --mem=8G

# Load the modules
module load gcc/12.3.0
module load R/4.5.2

# Launch R
R
```

Inside R, run the following to install all packages needed for Modules 03 and 06:

```r
# Set your personal library path (R will prompt you the first time — say yes)
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

# CRAN packages
install.packages(
  c("ggplot2", "dplyr", "tidyr", "ggrepel", "RColorBrewer",
    "pheatmap", "scales"),
  lib = lib,
  repos = "https://cloud.r-project.org"
)

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")

BiocManager::install(
  c("DESeq2", "BiocGenerics", "S4Vectors", "IRanges",
    "GenomicRanges", "SummarizedExperiment"),
  lib = lib,
  update = FALSE,
  ask = FALSE
)

# Confirm all packages load correctly
pkgs <- c("DESeq2", "ggplot2", "dplyr", "tidyr",
          "ggrepel", "pheatmap", "RColorBrewer")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
cat("\nPackage availability:\n")
print(ok)
```

Exit R and the interactive session when done:

```r
q(save = "no")
```

```bash
exit
```

> **Tip:** If `BiocManager::install()` asks whether to update existing packages, entering `n` (no) is fine — you just need the packages present, not necessarily the latest versions.

> **Q0 (setup check):** After installation, which version of DESeq2 was installed? Run `packageVersion("DESeq2")` inside R to check. Why does it matter that the R version (4.5.2) and Bioconductor release are compatible?

---

## 3.1 — Prepare the DESeq2 R Script

Create the R script in your scripts directory. This script is called by the SLURM job in section 3.2.

```bash
cat > ${TUTORIAL}/scripts/deseq2_peltaster.R << 'REOF'
# ============================================================
# DESeq2 Differential Expression Analysis
# P. fructicola: WT vs gh31del
# HCS 7004 — Genome Analytics
# ============================================================

# Ensure the personal R library is on the search path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ----------------------------------------------------------
# 1. Load count matrix
# ----------------------------------------------------------
counts_file <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/counts/peltaster_counts.txt"

raw <- read.table(counts_file, header = TRUE, skip = 1, sep = "\t",
                  row.names = 1, check.names = FALSE)

# featureCounts: first 5 columns after Geneid are Chr/Start/End/Strand/Length
# Count columns start at column 6
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

cat("\nSize factors (library normalization):\n")
print(sizeFactors(dds))

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
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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
print(head(sig[, c("gene_id", "baseMean", "log2FoldChange", "padj")], 20))

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

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                  color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("WT" = "#2196F3", "gh31del" = "#E53935")) +
  labs(title = "PCA of VST-normalized counts",
       subtitle = "P. fructicola WT vs gh31del",
       x = paste0("PC1: ", pca_var[1], "% variance"),
       y = paste0("PC2: ", pca_var[2], "% variance"),
       color = "Condition") +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "pca_plot.pdf"), pca_plot, width = 7, height = 5)
cat("PCA plot saved.\n")

# ----------------------------------------------------------
# 7. Volcano plot
# ----------------------------------------------------------
res_df$significance <- "Not significant"
res_df$significance[!is.na(res_df$padj) &
                    res_df$padj < 0.05 &
                    res_df$log2FoldChange >= 1]  <- "Up in WT"
res_df$significance[!is.na(res_df$padj) &
                    res_df$padj < 0.05 &
                    res_df$log2FoldChange <= -1] <- "Down in WT"

# Genes to highlight
genes_of_interest <- c("t1.AMS68_008039",  # GH31 itself
                       "t1.AMS68_000995")  # GH31 master regulator target

volcano <- ggplot(res_df %>% filter(!is.na(padj)),
                  aes(x = log2FoldChange, y = -log10(padj),
                      color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Not significant" = "grey60",
                                "Up in WT"        = "#D73027",
                                "Down in WT"      = "#4575B4")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  geom_point(data = res_df %>% filter(gene_id %in% genes_of_interest),
             size = 3, shape = 21, color = "black", fill = "gold") +
  geom_label_repel(
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
cat("Volcano plot saved.\n")

# ----------------------------------------------------------
# 8. MA plot
# ----------------------------------------------------------
pdf(file.path(outdir, "ma_plot.pdf"), width = 7, height = 5)
plotMA(res, main = "MA Plot: WT vs gh31del", alpha = 0.05, ylim = c(-6, 6))
dev.off()
cat("MA plot saved.\n")

# ----------------------------------------------------------
# 9. Dispersion estimates
# ----------------------------------------------------------
pdf(file.path(outdir, "dispersion_plot.pdf"), width = 7, height = 5)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()
cat("Dispersion plot saved.\n")

# ----------------------------------------------------------
# 10. Sample distance heatmap
# ----------------------------------------------------------
sampleDists       <- dist(t(assay(vst_data)))
sampleDistMatrix  <- as.matrix(sampleDists)

pdf(file.path(outdir, "sample_distance_heatmap.pdf"), width = 6, height = 5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main  = "Sample-to-Sample Distances (VST)",
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()
cat("Sample distance heatmap saved.\n")

cat("\nAll output files written to:", outdir, "\n")
cat("DESeq2 analysis complete.\n")
REOF
```

---

## 3.2 — Run DESeq2 on OSC

The SLURM script loads the OSC R modules and calls the script directly — no container needed.

```bash
cat > ${TUTORIAL}/scripts/03_deseq2.sh << 'EOF'
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

# Load OSC R modules
module load gcc/12.3.0
module load R/4.5.2

echo "[$(date)] R version: $(R --version | head -1)"
echo "[$(date)] Running DESeq2..."

Rscript ${TUTORIAL}/scripts/deseq2_peltaster.R

echo "[$(date)] DESeq2 complete."
ls -lh ${TUTORIAL}/dge/
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/03_deseq2.sh
```

Monitor the job log in real time while it runs:

```bash
tail -f ${TUTORIAL}/logs/deseq2_*.log
```

A successful run ends with lines like:

```
Sample distance heatmap saved.

All output files written to: /fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/dge
DESeq2 analysis complete.
```

If you see `Error in library(DESeq2) : there is no package called 'DESeq2'`, the package installation in section 3.0 did not complete — return to that step.

---

## 3.3 — Explore the Results

After the job completes, download the output files to your local machine:

```bash
# Run on your LOCAL terminal (not on OSC)
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/hcs7004/<your_username>/rnaseq_atacseq/dge/*.pdf ~/Desktop/
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/hcs7004/<your_username>/rnaseq_atacseq/dge/*.csv ~/Desktop/
```

Open `pca_plot.pdf`, `volcano_plot.pdf`, `ma_plot.pdf`, and `dispersion_plot.pdf`.

Then on OSC, examine the significant DEG list:

```bash
# How many significant DEGs?
awk -F',' 'NR>1 && $7!="NA" && $7<0.05 && ($6>1 || $6<-1) {print}' \
    ${TUTORIAL}/dge/peltaster_deseq2_results.csv | wc -l

# Top upregulated genes in WT (sorted by log2FC)
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

The dispersion plot was already generated in section 3.1 (step 9 of the R script) and saved as `dispersion_plot.pdf`. Open it alongside your other outputs.

> **Q22:** In the dispersion plot, do the final (shrunken) dispersions follow the fitted trend line reasonably well? What would extremely scattered dispersions suggest about the quality of the replicates or the appropriateness of the model?

---

## Summary

By the end of this module you should have:

- [x] DESeq2 results CSV with all genes tested
- [x] Significant DEGs CSV (padj < 0.05, |log2FC| ≥ 1)
- [x] Volcano plot, MA plot, PCA plot, dispersion plot, sample distance heatmap
- [x] A biological hypothesis about GH31's regulatory role

**Proceed to [Module 04 — ATAC-seq Processing](04_atacseq_processing.md)**
