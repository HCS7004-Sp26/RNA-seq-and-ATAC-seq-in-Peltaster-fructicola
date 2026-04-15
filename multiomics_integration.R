#!/usr/bin/env Rscript
# ============================================================
# Multi-Omics Integration: RNA-seq + ATAC-seq
# Peltaster fructicola — GH31 regulatory analysis
# HCS 7004 — Genome Analytics
#
# Usage: Rscript multiomics_integration.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# --- Paths ---
BASE      <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq"
DGEDIR    <- file.path(BASE, "dge")
ATACDIR   <- file.path(BASE, "atacseq/peaks/consensus")
OUTDIR    <- file.path(BASE, "integration")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 1. Load DESeq2 results
# ----------------------------------------------------------
deseq2 <- read.csv(file.path(DGEDIR, "peltaster_deseq2_results.csv"),
                   stringsAsFactors = FALSE)
cat(sprintf("DESeq2 results: %d genes tested\n", nrow(deseq2)))

# Clean gene IDs for joining (strip "t1." prefix if needed)
deseq2$gene_id_clean <- sub("^t1\\.", "", deseq2$gene_id)

# ----------------------------------------------------------
# 2. Load genes with promoter ATAC peaks
# ----------------------------------------------------------
peaks_in_promoters_file <- file.path(ATACDIR, "peaks_in_promoters.bed")

if (file.exists(peaks_in_promoters_file)) {
  peak_genes_raw <- read.table(peaks_in_promoters_file,
                                header = FALSE, sep = "\t")
  # Column 4 is gene_id from the promoter BED
  peak_gene_ids <- unique(peak_genes_raw[, 4])
  peak_gene_ids_clean <- sub("^t1\\.", "", peak_gene_ids)
  cat(sprintf("Genes with promoter ATAC peak: %d\n", length(peak_gene_ids_clean)))
} else {
  cat("WARNING: peaks_in_promoters.bed not found. Run Module 05 first.\n")
  peak_gene_ids_clean <- character(0)
}

# ----------------------------------------------------------
# 3. Merge datasets
# ----------------------------------------------------------
merged <- deseq2 %>%
  mutate(has_atac_peak = gene_id_clean %in% peak_gene_ids_clean)

cat(sprintf("Genes with ATAC peak AND expression data: %d\n",
            sum(merged$has_atac_peak & !is.na(merged$padj))))

# ----------------------------------------------------------
# 4. Categorize genes
# ----------------------------------------------------------
merged <- merged %>%
  mutate(
    deg_category = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up in WT",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down in WT",
      TRUE ~ "Background"
    )
  )

# ATAC peak rate by DEG category
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
write.csv(atac_by_deg, file.path(OUTDIR, "atac_rate_by_deg_category.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 5. Statistical test: Fisher's exact test
#    Up-in-WT vs. Background — are peaks enriched at up-regulated genes?
# ----------------------------------------------------------
up_peak     <- sum(merged$deg_category == "Up in WT"  & merged$has_atac_peak,  na.rm = TRUE)
up_nopeak   <- sum(merged$deg_category == "Up in WT"  & !merged$has_atac_peak, na.rm = TRUE)
bg_peak     <- sum(merged$deg_category == "Background" & merged$has_atac_peak, na.rm = TRUE)
bg_nopeak   <- sum(merged$deg_category == "Background" & !merged$has_atac_peak, na.rm = TRUE)

contingency <- matrix(c(up_peak, up_nopeak, bg_peak, bg_nopeak),
                      nrow = 2,
                      dimnames = list(
                        c("Has ATAC peak", "No ATAC peak"),
                        c("Up in WT", "Background")
                      ))

cat("\n=== Fisher's Exact Test: Up-in-WT vs Background ===\n")
print(contingency)
ft <- fisher.test(contingency)
cat(sprintf("Odds ratio: %.2f\n", ft$estimate))
cat(sprintf("p-value: %.4g\n", ft$p.value))

# ----------------------------------------------------------
# 6. Visualization: stacked bar chart
# ----------------------------------------------------------
bar_df <- atac_by_deg %>%
  mutate(pct_no_peak = 100 - pct_peak) %>%
  tidyr::pivot_longer(cols = c("pct_peak", "pct_no_peak"),
                      names_to = "peak_status",
                      values_to = "pct") %>%
  mutate(peak_status = ifelse(peak_status == "pct_peak",
                              "Has ATAC peak", "No ATAC peak"))

p_bar <- ggplot(bar_df,
       aes(x = deg_category, y = pct, fill = peak_status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Has ATAC peak" = "#1976D2",
                               "No ATAC peak"  = "#E0E0E0")) +
  geom_text(data = atac_by_deg,
            aes(x = deg_category, y = pct_peak + 3,
                label = paste0(pct_peak, "%"),
                fill = NULL),
            size = 4) +
  labs(title = "Chromatin accessibility at DEG promoters",
       subtitle = "Fraction of genes with promoter ATAC-seq peak",
       x = "Expression category (WT vs. gh31del)",
       y = "Percentage of genes (%)",
       fill = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "bar_atac_by_deg_category.pdf"),
       p_bar, width = 7, height = 5)
cat("\nBar chart saved.\n")

# ----------------------------------------------------------
# 7. Scatter: log2FC vs. presence/absence of ATAC peak
# ----------------------------------------------------------
p_violin <- ggplot(
  merged %>% filter(!is.na(log2FoldChange)),
  aes(x = has_atac_peak, y = log2FoldChange,
      fill = has_atac_peak)
) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.08, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_x_discrete(labels = c("FALSE" = "No peak", "TRUE" = "Has peak")) +
  scale_fill_manual(values = c("FALSE" = "#B0BEC5", "TRUE" = "#1565C0"),
                    guide = "none") +
  labs(
    title    = "Expression change vs. chromatin accessibility",
    subtitle = "P. fructicola WT vs. gh31del",
    x        = "Promoter ATAC-seq peak",
    y        = expression(log[2]~"Fold Change (WT / gh31del)")
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUTDIR, "violin_atac_vs_lfc.pdf"),
       p_violin, width = 6, height = 5)
cat("Violin plot saved.\n")

# ----------------------------------------------------------
# 8. GH31 and its key targets
# ----------------------------------------------------------
genes_of_interest <- c("AMS68_008039", "AMS68_000995")
cat("\n=== Key Genes of Interest ===\n")
for (g in genes_of_interest) {
  row <- merged %>% filter(grepl(g, gene_id))
  if (nrow(row) > 0) {
    cat(sprintf("\nGene: %s\n", g))
    cat(sprintf("  log2FC (WT/del): %.3f\n", row$log2FoldChange))
    cat(sprintf("  padj:            %.3e\n", row$padj))
    cat(sprintf("  Has ATAC peak:   %s\n",  row$has_atac_peak))
    cat(sprintf("  DEG category:    %s\n",  row$deg_category))
  } else {
    cat(sprintf("\nGene: %s — not found in DESeq2 results\n", g))
  }
}

# ----------------------------------------------------------
# 9. Write final integrated table
# ----------------------------------------------------------
final_table <- merged %>%
  filter(!is.na(padj)) %>%
  select(gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         has_atac_peak, deg_category) %>%
  arrange(padj)

write.csv(final_table,
          file.path(OUTDIR, "integrated_rnaseq_atacseq.csv"),
          row.names = FALSE)

cat(sprintf("\nFull integrated table written: %d genes\n", nrow(final_table)))
cat("Integration analysis complete.\n")
