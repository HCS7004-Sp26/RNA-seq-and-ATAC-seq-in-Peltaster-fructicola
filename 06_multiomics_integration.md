# Module 06 — Multi-Omics Integration: Linking Chromatin Accessibility to Gene Expression

## Overview

This is the culminating module of the tutorial. You now have two complementary datasets:

- **RNA-seq DEGs**: Genes whose expression changes when GH31 is deleted
- **ATAC-seq peaks**: Regions of open chromatin in wild-type *P. fructicola*

The central question is: **Are differentially expressed genes preferentially located near accessible chromatin?** And more specifically: **Does the GH31 locus itself show exceptional chromatin accessibility, and what does its regulatory landscape look like?**

In this module you will:

1. Intersect DEGs with ATAC-seq peaks
2. Ask whether DEG promoters are enriched for open chromatin
3. Examine the GH31 genomic locus in detail
4. Synthesize biological findings into a regulatory model

---

## 6.1 — Prepare Input Files

You need two key files from previous modules:

```bash
# Significant DEGs from Module 03
DEG_FILE=${TUTORIAL}/dge/peltaster_sig_DEGs.csv

# Consensus ATAC peaks from Module 05
PEAKS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed

# Reference annotation
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
```

First, create a BED file of promoter regions for all DEGs. You need to look up the genomic coordinates of each significant DEG in the GFF3:

```bash
cat > ${TUTORIAL}/scripts/deg_promoters.sh << 'EOF'
#!/bin/bash
set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
DEG_FILE=${TUTORIAL}/dge/peltaster_sig_DEGs.csv
OUTDIR=${TUTORIAL}/integration

mkdir -p ${OUTDIR}

# Extract gene IDs from DEG CSV (column 1, skip header)
tail -n +2 ${DEG_FILE} | cut -d',' -f1 | sed 's/"//g' \
    > ${OUTDIR}/deg_gene_ids.txt

echo "Total significant DEGs: $(wc -l < ${OUTDIR}/deg_gene_ids.txt)"

# Extract promoter coordinates for DEGs from GFF3
# Match gene IDs (accounting for the "t1." prefix in GFF3 gene_id attribute)
while IFS= read -r gene; do
    grep -P "\tgene\t" ${GFF} | grep -F "${gene}" 
done < ${OUTDIR}/deg_gene_ids.txt \
| awk 'BEGIN{OFS="\t"} {
    match($9, /gene_id=([^;]+)/, arr)
    if ($7=="+") print $1, ($4>2000 ? $4-2001 : 0), $4, arr[1], ".", $7
    else          print $1, $5-1, $5+2000, arr[1], ".", $7
}' | sort -k1,1 -k2,2n \
> ${OUTDIR}/deg_promoters.bed

echo "DEG promoter regions: $(wc -l < ${OUTDIR}/deg_promoters.bed)"
EOF

bash ${TUTORIAL}/scripts/deg_promoters.sh
```

---

## 6.2 — Test for Enrichment of ATAC Peaks at DEG Promoters

Do DEG promoters have significantly more ATAC-seq peaks than expected by chance? We test this by comparing peak overlap at DEG promoters to a random set of background gene promoters.

```bash
cat > ${TUTORIAL}/06_deg_peak_overlap.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=deg_peak_overlap
#SBATCH --account=PAS3260
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/deg_peak_overlap_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
BEDTOOLS_SIF=${CLASSDATA}/containers/bedtools.sif
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
PEAKS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed
OUTDIR=${TUTORIAL}/integration

# ---- DEG promoters vs. ATAC peaks ----
echo "=== DEG Promoters with ATAC Peaks ==="
apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${OUTDIR}/deg_promoters.bed \
        -b ${PEAKS} \
        -u \
> ${OUTDIR}/deg_promoters_with_peaks.bed

DEG_TOTAL=$(wc -l < ${OUTDIR}/deg_promoters.bed)
DEG_WITH_PEAKS=$(wc -l < ${OUTDIR}/deg_promoters_with_peaks.bed)
DEG_PCT=$(echo "scale=1; 100*${DEG_WITH_PEAKS}/${DEG_TOTAL}" | bc)
echo "DEG promoters: ${DEG_TOTAL}"
echo "DEG promoters with ATAC peak: ${DEG_WITH_PEAKS} (${DEG_PCT}%)"

# ---- All gene promoters vs. ATAC peaks (background) ----
echo ""
echo "=== All Gene Promoters with ATAC Peaks (background) ==="

grep -P "\tgene\t" ${GFF} \
    | awk 'BEGIN{OFS="\t"} {
        match($9, /gene_id=([^;]+)/, arr)
        if ($7=="+") print $1, ($4>2000 ? $4-2001 : 0), $4, arr[1], ".", $7
        else          print $1, $5-1, $5+2000, arr[1], ".", $7
      }' | sort -k1,1 -k2,2n \
> ${OUTDIR}/all_gene_promoters.bed

apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${OUTDIR}/all_gene_promoters.bed \
        -b ${PEAKS} \
        -u \
> ${OUTDIR}/all_promoters_with_peaks.bed

ALL_TOTAL=$(wc -l < ${OUTDIR}/all_gene_promoters.bed)
ALL_WITH_PEAKS=$(wc -l < ${OUTDIR}/all_promoters_with_peaks.bed)
ALL_PCT=$(echo "scale=1; 100*${ALL_WITH_PEAKS}/${ALL_TOTAL}" | bc)
echo "All gene promoters: ${ALL_TOTAL}"
echo "All promoters with ATAC peak: ${ALL_WITH_PEAKS} (${ALL_PCT}%)"
echo ""
echo "Enrichment (DEG% / All%): $(echo "scale=2; ${DEG_PCT}/${ALL_PCT}" | bc)"
EOF

sbatch ${TUTORIAL}/06_deg_peak_overlap.sh
```

---

## 6.3 — Examine the GH31 Locus and Its Regulatory Region

Now zoom in on GH31 (AMS68_008039) and its immediate genomic neighborhood. GH31 is located on chromosome CP051143.1 at coordinates 2,682,570–2,685,565 (+ or − strand — check the GFF3).

```bash
# Find GH31 in the annotation
grep "AMS68_008039" ${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3 | head -5

# Check for ATAC peaks near the GH31 locus
# Define a 5 kb window around GH31
echo -e "CP051143.1\t2677570\t2690565\tGH31_region" \
    > ${TUTORIAL}/integration/gh31_region.bed

apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${CLASSDATA}/containers/bedtools.sif \
    bedtools intersect \
        -a ${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed \
        -b ${TUTORIAL}/integration/gh31_region.bed \
        -wa \
> ${TUTORIAL}/integration/gh31_nearby_peaks.bed

echo "ATAC peaks within 5 kb of GH31:"
cat ${TUTORIAL}/integration/gh31_nearby_peaks.bed
```

Also check the most accessible peak in the entire dataset — the peak with the highest MACS3 score:

```bash
# Combine all three replicates' narrowPeak scores and find the maximum
cat ${TUTORIAL}/atacseq/peaks/atac_wt_rep*/atac_wt_rep*_peaks.narrowPeak \
    | sort -k9,9rn | head -10
```

The `-log10(q-value)` score is in column 9 of the narrowPeak format. Higher scores indicate higher confidence peaks.

> **Q31:** Does the GH31 promoter region have an ATAC-seq peak? Based on the MACS3 score, how does this peak rank compared to all peaks in the genome? What does exceptional chromatin accessibility at the GH31 locus suggest about its regulation?

---

## 6.4 — Integrate DEG and ATAC Data in R

Now perform a formal integration analysis in R: for every gene that has both expression data and an ATAC peak in its promoter, examine whether accessibility correlates with expression change.

```bash
cat > ${TUTORIAL}/scripts/multiomics_integration.R << 'REOF'
# ============================================================
# Multi-Omics Integration: RNA-seq + ATAC-seq
# P. fructicola — GH31 regulatory analysis
# HCS 7004 — Genome Analytics
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUTDIR <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/integration"
DGEDIR <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/dge"

# ----------------------------------------------------------
# 1. Load DESeq2 results
# ----------------------------------------------------------
deseq2 <- read.csv(file.path(DGEDIR, "peltaster_deseq2_results.csv"),
                   stringsAsFactors = FALSE)

cat("Total genes tested:", nrow(deseq2), "\n")

# ----------------------------------------------------------
# 2. Load ATAC peak annotation
# (which gene promoters have peaks)
# ----------------------------------------------------------
# promoters_with_peaks.bed has: chr, start, end, gene_id, ., strand
peak_genes <- read.table(
    file.path("/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq/peaks/consensus",
              "peaks_in_promoters.bed"),
    header = FALSE, sep = "\t",
    col.names = c("chr","start","end","gene_id","score","strand",
                  "peak_chr","peak_start","peak_end")
  )[, "gene_id", drop = FALSE]

peak_genes <- peak_genes %>%
  mutate(has_atac_peak = TRUE) %>%
  distinct()

cat("Genes with promoter ATAC peak:", nrow(peak_genes), "\n")

# ----------------------------------------------------------
# 3. Merge RNA-seq and ATAC data
# ----------------------------------------------------------
merged <- deseq2 %>%
  left_join(peak_genes, by = "gene_id") %>%
  mutate(has_atac_peak = ifelse(is.na(has_atac_peak), FALSE, TRUE))

cat("Genes with both RNA-seq and ATAC data:",
    sum(!is.na(merged$padj) & merged$has_atac_peak), "\n")

# ----------------------------------------------------------
# 4. Compare expression changes: genes with vs. without peaks
# ----------------------------------------------------------
merged_filtered <- merged %>%
  filter(!is.na(log2FoldChange), !is.na(padj))

p_boxplot <- ggplot(merged_filtered,
       aes(x = has_atac_peak, y = log2FoldChange,
           fill = has_atac_peak)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("FALSE" = "#AAAAAA", "TRUE" = "#2196F3")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Expression change vs. ATAC-seq accessibility",
       subtitle = "P. fructicola WT vs gh31del",
       x = "Promoter has ATAC peak",
       y = expression(log[2]~"Fold Change (WT / gh31del)"),
       fill = "ATAC peak") +
  theme_bw(base_size = 12)

ggsave(file.path(OUTDIR, "boxplot_atac_vs_expression.pdf"),
       p_boxplot, width = 6, height = 5)

# Wilcoxon rank-sum test
wt <- wilcox.test(
  log2FoldChange ~ has_atac_peak,
  data = merged_filtered,
  alternative = "two.sided"
)
cat("\nWilcoxon test (log2FC: peak vs. no peak):\n")
cat(sprintf("  W = %.0f, p = %.4g\n", wt$statistic, wt$p.value))

# ----------------------------------------------------------
# 5. Focus on significant DEGs: are they enriched for ATAC peaks?
# ----------------------------------------------------------
sig_threshold <- 0.05
fc_threshold  <- 1.0

sig_degs <- merged_filtered %>%
  filter(padj < sig_threshold, abs(log2FoldChange) >= fc_threshold)

up_in_wt   <- sig_degs %>% filter(log2FoldChange > 0)
down_in_wt <- sig_degs %>% filter(log2FoldChange < 0)

cat("\n--- Significant DEGs with ATAC peaks ---\n")
cat(sprintf("Upregulated in WT: %d total, %d with ATAC peak (%.0f%%)\n",
            nrow(up_in_wt),
            sum(up_in_wt$has_atac_peak),
            100*mean(up_in_wt$has_atac_peak)))
cat(sprintf("Downregulated in WT: %d total, %d with ATAC peak (%.0f%%)\n",
            nrow(down_in_wt),
            sum(down_in_wt$has_atac_peak),
            100*mean(down_in_wt$has_atac_peak)))
cat(sprintf("Background (not DEG): %d total, %d with ATAC peak (%.0f%%)\n",
            nrow(merged_filtered) - nrow(sig_degs),
            sum(!sig_degs$has_atac_peak & merged_filtered$padj >= sig_threshold),
            100*mean(merged_filtered$has_atac_peak[merged_filtered$padj >= sig_threshold |
                                                    is.na(merged_filtered$padj)])))

# ----------------------------------------------------------
# 6. The GH31 locus: expression and accessibility
# ----------------------------------------------------------
gh31_gene   <- "AMS68_008039"
gh31_master <- "AMS68_000995"   # The key downstream regulator

for (goi in c(gh31_gene, gh31_master)) {
  row <- merged %>% filter(grepl(goi, gene_id))
  if (nrow(row) > 0) {
    cat(sprintf("\n%s:\n", goi))
    cat(sprintf("  log2FC (WT/del) = %.2f\n", row$log2FoldChange))
    cat(sprintf("  padj = %.2e\n", row$padj))
    cat(sprintf("  Has ATAC peak = %s\n", row$has_atac_peak))
  }
}

# ----------------------------------------------------------
# 7. Summary table of key DEGs
# ----------------------------------------------------------
key_genes <- sig_degs %>%
  arrange(desc(abs(log2FoldChange))) %>%
  select(gene_id, log2FoldChange, padj, has_atac_peak) %>%
  head(20)

cat("\n--- Top 20 DEGs by |log2FC| ---\n")
print(key_genes)

write.csv(key_genes,
          file.path(OUTDIR, "top_degs_with_atac.csv"),
          row.names = FALSE)

cat("\nIntegration analysis complete.\n")
REOF
```

Submit the integration analysis:

```bash
cat > ${TUTORIAL}/06_multiomics.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=multiomics
#SBATCH --account=PAS3260
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/multiomics_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/deseq2.sif

apptainer exec --bind ${TUTORIAL} ${SIF} \
    Rscript ${TUTORIAL}/scripts/multiomics_integration.R

echo "[$(date)] Integration analysis complete."
ls -lh ${TUTORIAL}/integration/
EOF

sbatch ${TUTORIAL}/06_multiomics.sh
```

---

## 6.5 — Synthesize Your Findings

By this point, you have:

- A list of genes that change expression in the *gh31* deletion
- A map of open chromatin across the WT *P. fructicola* genome
- An understanding of which DEGs have accessible promoters

> **Q32:** Is the percentage of DEG promoters with ATAC peaks significantly higher than the background rate for all genes? Report the enrichment ratio and interpret what this means for the relationship between chromatin accessibility and transcriptional response to GH31 deletion.

> **Q33:** Look at the top upregulated gene in WT (from Q18). Does it have an ATAC-seq peak in its promoter? Based on your understanding of chromatin accessibility, what does this suggest about how this gene is regulated?

> **Q34:** Below is a simplified model of what you might observe. Fill in the expected relationships and evaluate whether your data support or contradict each prediction:

| Prediction | Your result | Supported? |
|---|---|---|
| GH31 promoter has strong ATAC peak (highly accessible) | | |
| Top upregulated gene in WT has promoter ATAC peak | | |
| Genes downregulated in WT have promoter ATAC peaks | | |
| Most background (non-DEG) genes lack promoter ATAC peaks | | |

> **Q35:** The GH31 gene encodes a secreted glycoside hydrolase — it functions outside the cell, degrading host plant cell wall polysaccharides. Yet deleting it causes transcriptional changes inside the fungus. Propose a signal transduction model explaining how loss of extracellular enzyme activity could feed back to change gene expression. Consider: nutrient sensing, product inhibition, host-derived signals, or two-component signaling pathways.

> **Q36 (Synthesis):** Based on your combined RNA-seq and ATAC-seq results, write a 200–300 word paragraph summarizing the regulatory role of GH31 in *P. fructicola*. Your paragraph should address: (1) which genes change in expression when GH31 is deleted, (2) the chromatin accessibility context of those genes, (3) what this tells us about GH31's indirect role in transcriptional regulation, and (4) what experiment you would perform next to test your model.

---

## 6.6 — Bonus: Visualize the GH31 Region in IGV

For a direct visual check of your data, load your results into the Integrative Genomics Viewer (IGV). Download to your local machine:

```bash
# On your LOCAL machine:
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/reference/peltaster_fructicola_genome.fa ~/Desktop/
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/reference/peltaster_fructicola_annotation.gff3 ~/Desktop/
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/aligned/rnaseq_wt_rep1.sorted.bam ~/Desktop/
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq/aligned/atac_wt_rep1.filtered.bam ~/Desktop/
```

In IGV:
1. Load the *P. fructicola* genome (FASTA + index)
2. Load the GFF3 annotation
3. Load the RNA-seq BAM (WT rep1)
4. Load the ATAC-seq BAM (WT rep1)
5. Navigate to `CP051143.1:2,680,000-2,690,000` (the GH31 locus region)

> **Q37:** Describe what you see in the IGV browser at the GH31 locus. Is the chromatin accessibility (ATAC coverage) highest at the promoter, gene body, or 3′ end? How does the RNA-seq read density compare between the WT and *gh31del* tracks?

---

## Summary

By the end of this module you should have:

- [x] Intersected DEG promoters with ATAC peaks
- [x] Computed enrichment of ATAC peaks at DEG vs. background promoters
- [x] Integration R analysis output and box plots
- [x] A synthesized biological model for GH31's regulatory role

---

## Overall Tutorial Summary

Congratulations on completing the RNA-seq and ATAC-seq analysis workflow! Here is what you accomplished across all seven modules:

| Module | Key Output |
|---|---|
| 00 Setup | Computing environment, containers, data |
| 01 RNA-seq QC | Trimmed reads, QC reports |
| 02 Alignment & Quantification | BAMs, count matrix, TPM values |
| 03 Differential Expression | DESeq2 results, volcano plot, significant DEGs |
| 04 ATAC-seq Processing | Filtered BAMs, fragment size distributions |
| 05 Peak Calling | Consensus ATAC peaks, FRiP scores, peak annotation |
| 06 Multi-Omics Integration | DEG–ATAC overlap, enrichment analysis, regulatory model |

The complete workflow demonstrates how combining transcriptomic (RNA-seq) and epigenomic (ATAC-seq) data can generate richer biological hypotheses than either dataset alone — a cornerstone of modern multi-omics approaches in plant pathology and beyond.
