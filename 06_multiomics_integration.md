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
cat > ${TUTORIAL}/scripts/06_deg_promoters.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=deg_promoters
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/deg_promoters_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
DEG_FILE=${TUTORIAL}/dge/peltaster_sig_DEGs.csv
OUTDIR=${TUTORIAL}/integration

mkdir -p ${OUTDIR}

# --- Step 1: Extract DEG gene IDs from the CSV ---
# Column 1, skip header, strip quotes and the "t1." prefix
# DESeq2 inherits the "t1." prefix from featureCounts/GTF gene IDs,
# but the GFF3-derived promoter BED uses bare IDs (e.g. AMS68_000995).
# Stripping the prefix here ensures the grep in Step 3 finds matches.
tail -n +2 ${DEG_FILE} \
    | cut -d',' -f1 \
    | sed 's/"//g' \
    | sed 's/^t1\.//' \
    > ${OUTDIR}/deg_gene_ids.txt

echo "Significant DEGs: $(wc -l < ${OUTDIR}/deg_gene_ids.txt)"

# --- Step 2: Build promoter BED for all genes in GFF3 ---
# Single pass through the GFF3 — much faster than one grep per gene
grep -P "\tgene\t" ${GFF} \
    | awk 'BEGIN{OFS="\t"} {
        match($9, /gene_id=([^;]+)/, arr)
        gid = arr[1]
        if ($7=="+") {
            s = ($4 > 2001 ? $4 - 2001 : 0)
            print $1, s, $4, gid, ".", $7
        } else {
            print $1, $5 - 1, $5 + 2000, gid, ".", $7
        }
    }' | sort -k1,1 -k2,2n \
    > ${OUTDIR}/all_gene_promoters_2kb.bed

echo "All gene promoters: $(wc -l < ${OUTDIR}/all_gene_promoters_2kb.bed)"

# --- Step 3: Subset to DEG promoters ---
# grep -F matches any line containing one of the DEG gene IDs
# The IDs may carry a "t1." prefix in the GFF3 — match flexibly
grep -F -f ${OUTDIR}/deg_gene_ids.txt \
    ${OUTDIR}/all_gene_promoters_2kb.bed \
    | sort -k1,1 -k2,2n \
    > ${OUTDIR}/deg_promoters.bed

echo "DEG promoter regions: $(wc -l < ${OUTDIR}/deg_promoters.bed)"

# --- Sanity check: warn if the counts differ noticeably ---
N_IDS=$(wc -l < ${OUTDIR}/deg_gene_ids.txt)
N_BED=$(wc -l < ${OUTDIR}/deg_promoters.bed)
if [ "${N_BED}" -lt "${N_IDS}" ]; then
    echo "WARNING: ${N_BED} promoter regions found for ${N_IDS} DEG IDs."
    echo "Some gene IDs may not match the GFF3 — check ID format consistency."
fi

echo "[$(date)] Done."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/06_deg_promoters.sh
```

After it runs, check the log and verify the counts look right before moving on to the bedtools intersection step:
```bash
cat ${TUTORIAL}/logs/deg_promoters_*.log
head -3 ${TUTORIAL}/integration/deg_promoters.bed
```
---

## 6.2 — Test for Enrichment of ATAC Peaks at DEG Promoters

Do DEG promoters have significantly more ATAC-seq peaks than expected by chance? We test this by comparing peak overlap at DEG promoters to a random set of background gene promoters.

```bash
cat > ${TUTORIAL}/scripts/06_deg_peak_overlap.sh << 'EOF'
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
BEDTOOLS_SIF=${TUTORIAL}/containers/bedtools.sif
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
PEAKS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed
OUTDIR=${TUTORIAL}/integration

# ---- DEG promoters vs. ATAC peaks ----
echo "=== DEG Promoters with ATAC Peaks ==="

apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${OUTDIR}/deg_promoters.bed \
        -b ${PEAKS} \
        -u \
> ${OUTDIR}/deg_promoters_with_peaks.bed

DEG_TOTAL=$(wc -l < ${OUTDIR}/deg_promoters.bed)
DEG_WITH_PEAKS=$(wc -l < ${OUTDIR}/deg_promoters_with_peaks.bed)
DEG_PCT=$(echo "scale=1; 100*${DEG_WITH_PEAKS}/${DEG_TOTAL}" | bc)
echo "DEG promoters:                ${DEG_TOTAL}"
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

apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${OUTDIR}/all_gene_promoters.bed \
        -b ${PEAKS} \
        -u \
> ${OUTDIR}/all_promoters_with_peaks.bed

ALL_TOTAL=$(wc -l < ${OUTDIR}/all_gene_promoters.bed)
ALL_WITH_PEAKS=$(wc -l < ${OUTDIR}/all_promoters_with_peaks.bed)
ALL_PCT=$(echo "scale=1; 100*${ALL_WITH_PEAKS}/${ALL_TOTAL}" | bc)
echo "All gene promoters:                ${ALL_TOTAL}"
echo "All promoters with ATAC peak:      ${ALL_WITH_PEAKS} (${ALL_PCT}%)"
echo ""
echo "Enrichment (DEG% / All%): $(echo "scale=2; ${DEG_PCT}/${ALL_PCT}" | bc)"

# ---- Write gene ID list for the R integration step ----
awk '{print $4}' ${OUTDIR}/all_promoters_with_peaks.bed \
    | sort -u > ${OUTDIR}/genes_with_promoter_peaks.txt

echo ""
echo "Gene IDs with promoter ATAC peak: $(wc -l < ${OUTDIR}/genes_with_promoter_peaks.txt)"
echo "[$(date)] Done."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/06_deg_peak_overlap.sh
```

---

## 6.3 — Examine the GH31 Locus and Its Regulatory Region

Now zoom in on GH31 (AMS68_008039) and its immediate genomic neighborhood. GH31 is located on chromosome CP051143.1 at coordinates 2,682,570–2,685,565 (+ or − strand — check the GFF3).

```bash
# Find GH31 in the annotation
grep "AMS68_008039" ${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3 | head -5
```

Check for ATAC peaks near the GH31 locus by defining a 5 kb window around the gene and intersecting it with the consensus peak set:

```bash
BEDTOOLS_SIF=${TUTORIAL}/containers/bedtools.sif

# Define a 5 kb window spanning the GH31 locus
echo -e "CP051143.1\t2677570\t2690565\tGH31_region" \
    > ${TUTORIAL}/integration/gh31_region.bed

# Find consensus peaks overlapping that window
apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed \
        -b ${TUTORIAL}/integration/gh31_region.bed \
        -wa \
> ${TUTORIAL}/integration/gh31_nearby_peaks.bed

echo "ATAC peaks within 5 kb of GH31:"
cat ${TUTORIAL}/integration/gh31_nearby_peaks.bed
```

Also check the most accessible peaks in the entire dataset — the peaks with the highest MACS3 score:

```bash
# Combine all three replicates' narrowPeak files and sort by -log10(q-value) (column 9)
cat ${TUTORIAL}/atacseq/peaks/atac_wt_rep1/atac_wt_rep1_peaks.narrowPeak \
    ${TUTORIAL}/atacseq/peaks/atac_wt_rep2/atac_wt_rep2_peaks.narrowPeak \
    ${TUTORIAL}/atacseq/peaks/atac_wt_rep3/atac_wt_rep3_peaks.narrowPeak \
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

# Ensure the personal R library is on the search path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# --- Paths ---
BASE   <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq"
DGEDIR <- file.path(BASE, "dge")
OUTDIR <- file.path(BASE, "integration")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 1. Load DESeq2 results
# ----------------------------------------------------------
deseq2 <- read.csv(file.path(DGEDIR, "peltaster_deseq2_results.csv"),
                   stringsAsFactors = FALSE)
cat(sprintf("DESeq2 results loaded: %d genes\n", nrow(deseq2)))

# Strip the "t1." prefix from gene IDs to match the GFF3-derived IDs
deseq2$gene_id_clean <- sub("^t1\\.", "", deseq2$gene_id)

# ----------------------------------------------------------
# 2. Load genes with promoter ATAC peaks
#    (produced by 06a_promoter_peak_overlap.sh)
# ----------------------------------------------------------
peak_genes_file <- file.path(OUTDIR, "genes_with_promoter_peaks.txt")

if (!file.exists(peak_genes_file)) {
  stop("genes_with_promoter_peaks.txt not found. Run 06_deg_peak_overlap.sh first.")
}

peak_gene_ids <- readLines(peak_genes_file)
peak_gene_ids <- sub("^t1\\.", "", trimws(peak_gene_ids))
cat(sprintf("Genes with promoter ATAC peak: %d\n", length(peak_gene_ids)))

# ----------------------------------------------------------
# 3. Merge and categorize
# ----------------------------------------------------------
merged <- deseq2 %>%
  mutate(
    has_atac_peak = gene_id_clean %in% peak_gene_ids,
    deg_category  = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up in WT",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down in WT",
      TRUE ~ "Background"
    )
  )

cat(sprintf("Genes with ATAC peak + expression data: %d\n",
            sum(merged$has_atac_peak & !is.na(merged$padj))))

# ----------------------------------------------------------
# 4. ATAC peak rate by DEG category
# ----------------------------------------------------------
atac_by_deg <- merged %>%
  filter(!is.na(padj)) %>%
  group_by(deg_category) %>%
  summarise(
    n_genes     = n(),
    n_with_peak = sum(has_atac_peak),
    pct_peak    = round(100 * mean(has_atac_peak), 1),
    .groups     = "drop"
  )

cat("\n=== ATAC Peak Rate by Expression Category ===\n")
print(atac_by_deg)
write.csv(atac_by_deg,
          file.path(OUTDIR, "atac_rate_by_deg_category.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 5. Fisher's exact test: Up-in-WT vs. Background
# ----------------------------------------------------------
up_peak   <- sum(merged$deg_category == "Up in WT"   & merged$has_atac_peak,  na.rm = TRUE)
up_nopeak <- sum(merged$deg_category == "Up in WT"   & !merged$has_atac_peak, na.rm = TRUE)
bg_peak   <- sum(merged$deg_category == "Background" & merged$has_atac_peak,  na.rm = TRUE)
bg_nopeak <- sum(merged$deg_category == "Background" & !merged$has_atac_peak, na.rm = TRUE)

ct <- matrix(c(up_peak, up_nopeak, bg_peak, bg_nopeak),
             nrow = 2,
             dimnames = list(
               c("Has ATAC peak", "No ATAC peak"),
               c("Up in WT", "Background")
             ))

cat("\n=== Fisher's Exact Test: Up-in-WT vs. Background ===\n")
print(ct)
ft <- fisher.test(ct)
cat(sprintf("Odds ratio: %.2f\n", ft$estimate))
cat(sprintf("p-value:    %.4g\n", ft$p.value))

# ----------------------------------------------------------
# 6. Stacked bar chart
# ----------------------------------------------------------
bar_df <- atac_by_deg %>%
  mutate(
    pct_no_peak  = 100 - pct_peak,
    deg_category = factor(deg_category,
                          levels = c("Up in WT", "Down in WT", "Background"))
  ) %>%
  pivot_longer(cols = c("pct_peak", "pct_no_peak"),
               names_to  = "peak_status",
               values_to = "pct") %>%
  mutate(peak_status = ifelse(peak_status == "pct_peak",
                              "Has ATAC peak", "No ATAC peak"))

p_bar <- ggplot(bar_df,
                aes(x = deg_category, y = pct, fill = peak_status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Has ATAC peak" = "#1976D2",
                               "No ATAC peak"  = "#E0E0E0")) +
  geom_text(data = atac_by_deg %>%
              mutate(deg_category = factor(deg_category,
                                           levels = c("Up in WT","Down in WT","Background"))),
            aes(x = deg_category, y = pct_peak + 3,
                label = paste0(pct_peak, "%"), fill = NULL),
            size = 4) +
  labs(title    = "Chromatin accessibility at DEG promoters",
       subtitle = "Fraction of genes with promoter ATAC-seq peak",
       x        = "Expression category (WT vs. gh31del)",
       y        = "Percentage of genes (%)",
       fill     = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "bar_atac_by_deg_category.pdf"),
       p_bar, width = 7, height = 5)
cat("Bar chart saved.\n")

# ----------------------------------------------------------
# 7. Violin + boxplot: log2FC vs. ATAC peak presence
# ----------------------------------------------------------
p_violin <- ggplot(
  merged %>% filter(!is.na(log2FoldChange)),
  aes(x = has_atac_peak, y = log2FoldChange, fill = has_atac_peak)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.08, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red", linewidth = 0.6) +
  scale_x_discrete(labels = c("FALSE" = "No peak", "TRUE" = "Has peak")) +
  scale_fill_manual(values = c("FALSE" = "#B0BEC5", "TRUE" = "#1565C0"),
                    guide  = "none") +
  labs(title    = "Expression change vs. chromatin accessibility",
       subtitle = "P. fructicola WT vs. gh31del",
       x        = "Promoter ATAC-seq peak",
       y        = expression(log[2]~"Fold Change (WT / gh31del)")) +
  theme_bw(base_size = 12)

ggsave(file.path(OUTDIR, "violin_atac_vs_lfc.pdf"),
       p_violin, width = 6, height = 5)
cat("Violin plot saved.\n")

# ----------------------------------------------------------
# 8. Genes of interest: GH31 and its key target
# ----------------------------------------------------------
genes_of_interest <- c("AMS68_008039", "AMS68_000995")
cat("\n=== Key Genes of Interest ===\n")
for (g in genes_of_interest) {
  row <- merged %>% filter(grepl(g, gene_id))
  if (nrow(row) > 0) {
    cat(sprintf("\nGene: %s\n", g))
    cat(sprintf("  log2FC (WT/del): %s\n",
                ifelse(is.na(row$log2FoldChange), "NA",
                       sprintf("%.3f", row$log2FoldChange))))
    cat(sprintf("  padj:            %s\n",
                ifelse(is.na(row$padj), "NA", sprintf("%.3e", row$padj))))
    cat(sprintf("  Has ATAC peak:   %s\n",  row$has_atac_peak))
    cat(sprintf("  DEG category:    %s\n",  row$deg_category))
  } else {
    cat(sprintf("\n%s: not found in DESeq2 results\n", g))
  }
}

# ----------------------------------------------------------
# 9. Write final integrated table
# ----------------------------------------------------------
final_table <- merged %>%
  filter(!is.na(padj)) %>%
  select(gene_id, baseMean, log2FoldChange, lfcSE,
         stat, pvalue, padj, has_atac_peak, deg_category) %>%
  arrange(padj)

write.csv(final_table,
          file.path(OUTDIR, "integrated_rnaseq_atacseq.csv"),
          row.names = FALSE)

cat(sprintf("\nFull integrated table: %d genes — written to integrated_rnaseq_atacseq.csv\n",
            nrow(final_table)))
cat("Integration analysis complete.\n")
REOF
```

Submit the integration analysis:

```bash
cat > ${TUTORIAL}/scripts/06b_multiomics_r.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=multiomics_r
#SBATCH --account=PAS3260
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/multiomics_r_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq

# Load OSC R modules
module load gcc/12.3.0
module load R/4.5.2

echo "[$(date)] Running multi-omics integration..."

Rscript ${TUTORIAL}/scripts/multiomics_integration.R

echo "[$(date)] Done."
ls -lh ${TUTORIAL}/integration/
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/06b_multiomics_r.sh
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

## 6.6 — Visualize the GH31 Region in IGV

If you have not yet loaded your data into IGV, follow the complete step-by-step guide in **[Module 05 — Section 5.6](05_peak_calling.md)**. That section covers downloading all required files, loading the custom genome, bigWig tracks, consensus peaks BED, and RNA-seq BAMs, and navigating to any locus of interest.

Once your session is set up, navigate to the GH31 locus:

```
CP051143.1:2,680,000-2,690,000
```

> **Q39:** Describe what you see in the IGV browser at the GH31 locus. Is the chromatin accessibility (ATAC coverage) highest at the promoter, gene body, or 3′ end? How does the RNA-seq read density compare between the WT and *gh31del* tracks?

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
