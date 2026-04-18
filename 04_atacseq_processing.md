# Module 04 — ATAC-seq Processing

## Overview

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) uses the bacterial Tn5 transposase to insert sequencing adapters into open chromatin regions. Because Tn5 cannot access nucleosome-bound DNA, the resulting reads reflect nucleosome-free, accessible chromatin. In this module you will:

1. Run FastQC and Trimmomatic on raw ATAC-seq reads
2. Align reads to the *P. fructicola* genome with Bowtie2
3. Apply post-alignment filters specific to ATAC-seq
4. Assess library quality (fragment size distribution, FRiP score)

---

## Background: ATAC-seq vs. RNA-seq Read Processing

ATAC-seq alignment differs from RNA-seq in several important ways:

| Aspect | RNA-seq (HISAT2) | ATAC-seq (Bowtie2) |
|---|---|---|
| Aligner | Splice-aware | Not splice-aware needed |
| Fragment size | Uniform (RNA size) | Bimodal (sub-nucleosomal + mono-nucleosomal) |
| Key filter | Keep uniquely mapped | Remove duplicates, keep properly paired |
| Mitochondrial | Not typically a concern | Must remove (high signal) |
| Insert size | Not diagnostic | Informative (~200 bp ladder) |

In fungal ATAC-seq, mitochondrial DNA (and rDNA) is often highly accessible and will dominate a library unless explicitly removed.

> **Q23:** Tn5 inserts adapters as a dimer, cutting DNA at two sites separated by ~9 bp. Why does this mean ATAC-seq reads need a +4/−5 bp offset correction for precise footprinting analyses? Is this offset correction necessary for calling broad accessibility peaks?

---

## 4.1 — FastQC and Trimmomatic for ATAC-seq

The same QC and trimming approach used for RNA-seq applies here, but ATAC-seq reads commonly show higher adapter content because short fragments from sub-nucleosomal regions may be entirely adapter sequence.

```bash
cat > ${TUTORIAL}/scripts/04_atacseq_qc_trim.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=atac_qc_trim
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-3
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/atac_qc_trim_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
INDIR=${TUTORIAL}/atacseq
TRIMDIR=${TUTORIAL}/atacseq/trimmed
FQCDIR=${TUTORIAL}/atacseq/fastqc_raw

mkdir -p ${TRIMDIR} ${FQCDIR}

SAMPLES=(atac_wt_rep1 atac_wt_rep2 atac_wt_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Processing ATAC-seq sample: ${SAMPLE}"

# FastQC on raw reads
apptainer exec --bind ${TUTORIAL} ${TUTORIAL}/containers/fastqc.sif \
    fastqc \
        --outdir ${FQCDIR} \
        --threads 4 \
        ${INDIR}/${SAMPLE}_R1.fastq.gz \
        ${INDIR}/${SAMPLE}_R2.fastq.gz

# Trimmomatic — ATAC-seq uses Nextera adapters
apptainer exec --bind ${TUTORIAL} ${TUTORIAL}/containers/trimmomatic.sif \
    trimmomatic PE \
        -threads 8 \
        ${INDIR}/${SAMPLE}_R1.fastq.gz \
        ${INDIR}/${SAMPLE}_R2.fastq.gz \
        ${TRIMDIR}/${SAMPLE}_R1_trimmed.fastq.gz \
        ${TRIMDIR}/${SAMPLE}_R1_unpaired.fastq.gz \
        ${TRIMDIR}/${SAMPLE}_R2_trimmed.fastq.gz \
        ${TRIMDIR}/${SAMPLE}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:/opt/conda/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa:2:30:10 \
        -phred33 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

echo "[$(date)] Done: ${SAMPLE}"
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/04_atacseq_qc_trim.sh
```

> **Q24:** ATAC-seq uses **Nextera** adapters rather than the TruSeq adapters used for RNA-seq. Why does the library preparation chemistry differ, and what property of the Tn5 transposase drives adapter choice?

---

## 4.2 — Build a Bowtie2 Index and Align ATAC-seq Reads

For ATAC-seq, we use Bowtie2 in **very sensitive** mode with a liberal insert size range to capture both sub-nucleosomal and di-nucleosomal fragments:

```bash
cat > ${TUTORIAL}/scripts/04_bowtie2_index.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=bowtie2_index
#SBATCH --account=PAS3260
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/bowtie2_index_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/bowtie2.sif

mkdir -p ${TUTORIAL}/reference/bowtie2_index

apptainer exec --bind ${TUTORIAL} ${SIF} \
    bowtie2-build \
        --threads 4 \
        ${TUTORIAL}/reference/peltaster_fructicola_genome.fa \
        ${TUTORIAL}/reference/bowtie2_index/peltaster

echo "[$(date)] Bowtie2 index complete."
ls -lh ${TUTORIAL}/reference/bowtie2_index/
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/04_bowtie2_index.sh
```

Now align all three ATAC-seq replicates:

```bash
cat > ${TUTORIAL}/scripts/04_bowtie2_align_atac.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=bowtie2_atac
#SBATCH --account=PAS3260
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=1-3
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/bowtie2_atac_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
BT2_SIF=${TUTORIAL}/containers/bowtie2.sif
SAM_SIF=${TUTORIAL}/containers/samtools.sif
INDEX=${TUTORIAL}/reference/bowtie2_index/peltaster
TRIMDIR=${TUTORIAL}/atacseq/trimmed
OUTDIR=${TUTORIAL}/atacseq/aligned

mkdir -p ${OUTDIR}

SAMPLES=(atac_wt_rep1 atac_wt_rep2 atac_wt_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Aligning ATAC-seq: ${SAMPLE}"

# Align — very-sensitive, no mixed/discordant, max insert 2000 bp
apptainer exec --bind ${TUTORIAL} ${BT2_SIF} \
    bowtie2 \
        -p 8 \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        -X 2000 \
        -x ${INDEX} \
        -1 ${TRIMDIR}/${SAMPLE}_R1_trimmed.fastq.gz \
        -2 ${TRIMDIR}/${SAMPLE}_R2_trimmed.fastq.gz \
        2>${TUTORIAL}/logs/bowtie2_${SAMPLE}.summary \
| apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
    samtools sort \
        -@ 8 -m 2G \
        -o ${OUTDIR}/${SAMPLE}.sorted.bam

apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
    samtools index ${OUTDIR}/${SAMPLE}.sorted.bam

echo "[$(date)] Done: ${SAMPLE}"
cat ${TUTORIAL}/logs/bowtie2_${SAMPLE}.summary
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/04_bowtie2_align_atac.sh
```

---

## 4.3 — Post-Alignment Filtering

ATAC-seq requires several filtering steps after alignment. Perform these in sequence:

**Step 1: Remove duplicates** (PCR and optical duplicates inflate accessibility signal)  
**Step 2: Filter low-quality and multi-mapped reads** (MAPQ < 30)  
**Step 3: Keep only properly paired reads**  
**Step 4: Remove mitochondrial reads** (if a mitochondrial contig is present)

```bash
cat > ${TUTORIAL}/scripts/04_atac_filter.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=atac_filter
#SBATCH --account=PAS3260
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-3
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/atac_filter_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/samtools.sif
INDIR=${TUTORIAL}/atacseq/aligned
TMPDIR=${TUTORIAL}/atacseq/tmp

mkdir -p ${TMPDIR}

SAMPLES=(atac_wt_rep1 atac_wt_rep2 atac_wt_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Processing: ${SAMPLE}"

# Step 1: Name-sort the coordinate-sorted BAM
# fixmate requires mates to be adjacent (name order)
echo "[$(date)] Step 1: Name-sort..."
apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools sort \
        -n \
        -@ 8 \
        -m 2G \
        -T ${TMPDIR}/${SAMPLE}_namesort \
        -o ${INDIR}/${SAMPLE}.namesorted.bam \
        ${INDIR}/${SAMPLE}.sorted.bam

# Step 2: fixmate — adds MC tag and mate score tags needed by markdup
echo "[$(date)] Step 2: fixmate..."
apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools fixmate \
        -m \
        -@ 8 \
        ${INDIR}/${SAMPLE}.namesorted.bam \
        ${INDIR}/${SAMPLE}.fixmate.bam

# Step 3: Re-sort by coordinate for markdup
echo "[$(date)] Step 3: Coordinate-sort after fixmate..."
apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools sort \
        -@ 8 \
        -m 2G \
        -T ${TMPDIR}/${SAMPLE}_coordsort \
        -o ${INDIR}/${SAMPLE}.fixmate.sorted.bam \
        ${INDIR}/${SAMPLE}.fixmate.bam

# Step 4: Mark and remove duplicates
echo "[$(date)] Step 4: markdup..."
apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools markdup \
        -r \
        -@ 8 \
        -s \
        -f ${TUTORIAL}/logs/markdup_${SAMPLE}.stats \
        ${INDIR}/${SAMPLE}.fixmate.sorted.bam \
        ${INDIR}/${SAMPLE}.dedup.bam

# Step 5: Filter — properly paired (-f 2), high MAPQ (-q 30),
#          remove unmapped/secondary/supplementary/QC-fail (-F 1804)
echo "[$(date)] Step 5: Filter..."
apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools view \
        -@ 8 \
        -b \
        -f 2 \
        -F 1804 \
        -q 30 \
        ${INDIR}/${SAMPLE}.dedup.bam \
| apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools sort \
        -@ 8 \
        -m 2G \
        -T ${TMPDIR}/${SAMPLE}_final \
        -o ${INDIR}/${SAMPLE}.filtered.bam

apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools index ${INDIR}/${SAMPLE}.filtered.bam

# Step 6: Read count summary at each stage
echo ""
echo "=== ${SAMPLE} read counts ==="
echo -n "Raw aligned:      "
apptainer exec --bind ${TUTORIAL} ${SIF} samtools view -c    ${INDIR}/${SAMPLE}.sorted.bam
echo -n "After fixmate:    "
apptainer exec --bind ${TUTORIAL} ${SIF} samtools view -c    ${INDIR}/${SAMPLE}.fixmate.sorted.bam
echo -n "After dedup:      "
apptainer exec --bind ${TUTORIAL} ${SIF} samtools view -c    ${INDIR}/${SAMPLE}.dedup.bam
echo -n "After filter:     "
apptainer exec --bind ${TUTORIAL} ${SIF} samtools view -c -f 2 ${INDIR}/${SAMPLE}.filtered.bam

echo ""
echo "Duplicate stats:"
cat ${TUTORIAL}/logs/markdup_${SAMPLE}.stats

# Step 7: Clean up large intermediates
echo "[$(date)] Cleaning up intermediates..."
rm -f ${INDIR}/${SAMPLE}.namesorted.bam \
      ${INDIR}/${SAMPLE}.fixmate.bam \
      ${INDIR}/${SAMPLE}.fixmate.sorted.bam \
      ${INDIR}/${SAMPLE}.dedup.bam

echo "[$(date)] Done: ${SAMPLE}"
echo "Final filtered BAM: ${INDIR}/${SAMPLE}.filtered.bam"
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/04_atac_filter.sh

# After the jobs are done, confirm all three filtered BAMs exist and are indexed
ls -lh ${TUTORIAL}/atacseq/aligned/*.filtered.bam
ls -lh ${TUTORIAL}/atacseq/aligned/*.filtered.bam.bai

# Check duplication rates
cat ${TUTORIAL}/logs/markdup_atac_wt_rep*.stats
```

---

## 4.4 — Fragment Size Distribution

A key quality metric for ATAC-seq is the fragment size distribution, which should show a characteristic nucleosomal ladder: a prominent sub-nucleosomal peak (~100–150 bp), a mono-nucleosomal peak (~200 bp), and decreasing di- and tri-nucleosomal peaks.

```bash
cat > ${TUTORIAL}/scripts/04_extract_insert_sizes.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=insert_sizes
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-3
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/insert_sizes_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/samtools.sif

SAMPLES=(atac_wt_rep1 atac_wt_rep2 atac_wt_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Extracting insert sizes: ${SAMPLE}"

# -f 2: properly paired reads only
# $9 is the TLEN (template length) field — positive for the leftmost read of each pair
apptainer exec --bind ${TUTORIAL} ${SIF} \
    samtools view -f 2 \
        ${TUTORIAL}/atacseq/aligned/${SAMPLE}.filtered.bam \
| awk '$9 > 0 && $9 < 1000 {print $9}' \
| sort -n \
| uniq -c \
> ${TUTORIAL}/atacseq/${SAMPLE}_insert_sizes.txt

echo "Lines written: $(wc -l < ${TUTORIAL}/atacseq/${SAMPLE}_insert_sizes.txt)"
echo "[$(date)] Done: ${SAMPLE}"
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/04_extract_insert_sizes.sh
```

After it finishes, verify the format of one output file:
```bash
head -5 ${TUTORIAL}/atacseq/atac_wt_rep1_insert_sizes.txt
```

### Write the R script

```bash
cat > ${TUTORIAL}/scripts/plot_fragment_sizes.R << 'REOF'
# ============================================================
# ATAC-seq Fragment Size Distribution
# P. fructicola wild type
# HCS 7004 — Genome Analytics
# ============================================================

# Ensure personal R library is on the search path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages(library(ggplot2))

# --- Paths ---
atacdir <- "/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq"
samples <- c("atac_wt_rep1", "atac_wt_rep2", "atac_wt_rep3")

# --- Load insert size tables ---
# uniq -c output: col1 = count, col2 = insert size
all_frags <- lapply(samples, function(s) {
  f <- file.path(atacdir, paste0(s, "_insert_sizes.txt"))
  if (!file.exists(f)) {
    message("WARNING: file not found — ", f)
    return(NULL)
  }
  d <- read.table(f, header = FALSE,
                  col.names = c("count", "size"))
  d <- d[d$size > 0 & d$size < 1000, ]
  d$sample <- s
  d
})

# Check that at least one file loaded
if (all(sapply(all_frags, is.null))) {
  stop("No insert size files found. Run 04_extract_insert_sizes.sh first.")
}

frag_df <- do.call(rbind, Filter(Negate(is.null), all_frags))

cat(sprintf("Total size classes loaded: %d across %d samples\n",
            nrow(frag_df), length(unique(frag_df$sample))))

# --- Plot ---
p <- ggplot(frag_df, aes(x = size, y = count, color = sample)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  scale_y_log10(
    labels = scales::label_comma()
  ) +
  scale_x_continuous(
    breaks = seq(0, 1000, 100),
    limits = c(0, 1000)
  ) +
  # Nucleosomal size reference lines
  geom_vline(xintercept = 150, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 200, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 400, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = 155, y = 10, label = "NFR",
           hjust = 0, size = 3, color = "grey30") +
  annotate("text", x = 205, y = 10, label = "Mono-nuc",
           hjust = 0, size = 3, color = "grey30") +
  annotate("text", x = 405, y = 10, label = "Di-nuc",
           hjust = 0, size = 3, color = "grey30") +
  labs(
    title    = "ATAC-seq Fragment Size Distribution",
    subtitle = "P. fructicola wild type (3 replicates)",
    x        = "Insert size (bp)",
    y        = "Read count (log10 scale)",
    color    = "Replicate"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

outfile <- file.path(atacdir, "fragment_size_distribution.pdf")
ggsave(outfile, p, width = 8, height = 5)
cat("Fragment size plot saved to:", outfile, "\n")
REOF
```

### Submit the SLURM job

```bash
cat > ${TUTORIAL}/scripts/04_plot_fragment_sizes.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=frag_size_plot
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/frag_size_plot_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq

module load gcc/12.3.0
module load R/4.5.2

echo "[$(date)] Plotting fragment size distributions..."

Rscript ${TUTORIAL}/scripts/plot_fragment_sizes.R

echo "[$(date)] Done."
ls -lh ${TUTORIAL}/atacseq/fragment_size_distribution.pdf
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/04_plot_fragment_sizes.sh
```

> **Q25:** Describe the fragment size distribution you observe. Do you see the expected nucleosomal ladder pattern? The sub-nucleosomal peak (< 150 bp) represents DNA fragments from nucleosome-free regions. The mono-nucleosomal peak (~200 bp) represents DNA wrapped around a single nucleosome that was not accessible. Why are both populations present in an ATAC-seq library?

> **Q26:** The *P. fructicola* genome is highly compact (~17 Mb). How might the compact genome organization of this species affect the ATAC-seq fragment size distribution compared to a typical mammalian ATAC-seq experiment?

---

## Summary

By the end of this module you should have:

- [x] FastQC/Trimmomatic output for all three ATAC-seq samples
- [x] Bowtie2 index in `reference/bowtie2_index/`
- [x] Sorted, indexed, filtered BAM files in `atacseq/aligned/`
- [x] Fragment size distribution plots

**Proceed to [Module 05 — Peak Calling](05_peak_calling.md)**
