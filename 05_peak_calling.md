# Module 05 — ATAC-seq Peak Calling with MACS3

## Overview

Peak calling identifies genomic regions where ATAC-seq read pileup is significantly higher than background — i.e., regions of open chromatin. In this module you will:

1. Call peaks on individual replicates using **MACS3**
2. Create a consensus peak set across replicates
3. Assess peak quality with FRiP scores
4. Annotate peaks relative to gene features

---

## Background: MACS3 and ATAC-seq Peak Calling

MACS3 (Model-based Analysis of ChIP-Seq, version 3) models the expected background signal and uses a Poisson-based approach to identify regions of significant enrichment. For ATAC-seq without a matched input/control library:

- Use `--nomodel` to skip the shift-model step (or let MACS3 estimate it)
- Use `--format BAMPE` for paired-end data (MACS3 uses insert sizes directly)
- The `--nolambda` flag can be used if local lambda estimation is not desired, but the default local lambda is generally recommended
- Expected peaks: ~1,000–5,000 for a compact fungal genome

---

## 5.1 — Per-Replicate Peak Calling

Call peaks on each of the three filtered BAM files:

```bash
cat > ${TUTORIAL}/scripts/05_macs3_individual.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=macs3_atac
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-3
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/macs3_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/macs3.sif
BAMDIR=${TUTORIAL}/atacseq/aligned
OUTDIR=${TUTORIAL}/atacseq/peaks

mkdir -p ${OUTDIR}

SAMPLES=(atac_wt_rep1 atac_wt_rep2 atac_wt_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Calling peaks for: ${SAMPLE}"

apptainer exec --bind ${TUTORIAL} ${SIF} \
    macs3 callpeak \
        --treatment ${BAMDIR}/${SAMPLE}.filtered.bam \
        --format BAMPE \
        --name ${SAMPLE} \
        --outdir ${OUTDIR}/${SAMPLE} \
        --gsize 1.7e7 \
        --qvalue 0.05 \
        --keep-dup all \
        --call-summits \
        2>&1 | tee ${TUTORIAL}/logs/macs3_${SAMPLE}.log

echo "[$(date)] Done: ${SAMPLE}"
echo "Peak count:"
wc -l ${OUTDIR}/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/05_macs3_individual.sh
```

After all three jobs finish, check peak counts:

```bash
for SAMPLE in atac_wt_rep1 atac_wt_rep2 atac_wt_rep3; do
    echo -n "${SAMPLE}: "
    wc -l < ${TUTORIAL}/atacseq/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
done
```

> **Q27:** How many peaks were called in each replicate? Given a genome size of ~17 Mb and ~8,100 genes, roughly what percentage of the genome falls within called peaks? Does this seem biologically reasonable for a compact fungal genome with tightly packed genes?

---

## 5.2 — Create a Consensus Peak Set

Individual replicates may call slightly different sets of peaks due to read sampling variation. A consensus peak set retains peaks supported by at least two of three replicates. We use **bedtools** for this step (already available as a container):

```bash
cat > ${TUTORIAL}/scripts/05_consensus_peaks.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=consensus_peaks
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/consensus_peaks_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
BEDTOOLS_SIF=${TUTORIAL}/containers/bedtools.sif
PEAKDIR=${TUTORIAL}/atacseq/peaks
OUTDIR=${TUTORIAL}/atacseq/peaks/consensus

mkdir -p ${OUTDIR}

REP1=${PEAKDIR}/atac_wt_rep1/atac_wt_rep1_peaks.narrowPeak
REP2=${PEAKDIR}/atac_wt_rep2/atac_wt_rep2_peaks.narrowPeak
REP3=${PEAKDIR}/atac_wt_rep3/atac_wt_rep3_peaks.narrowPeak

# Step 1: Pool all peaks and merge overlapping ones into a union set
echo "[$(date)] Building union peak set..."
cat ${REP1} ${REP2} ${REP3} \
    | sort -k1,1 -k2,2n \
    | apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
        bedtools merge -i stdin \
    > ${OUTDIR}/union_peaks.bed

echo "Union peak count: $(wc -l < ${OUTDIR}/union_peaks.bed)"

# Step 2: Identify which union peaks are supported by each replicate individually
# A union peak is "supported" by a replicate if it overlaps at least one peak
# called in that replicate's narrowPeak file
echo "[$(date)] Counting per-replicate support..."

for REP in 1 2 3; do
    REPFILE=${PEAKDIR}/atac_wt_rep${REP}/atac_wt_rep${REP}_peaks.narrowPeak
    apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
        bedtools intersect \
            -a ${OUTDIR}/union_peaks.bed \
            -b ${REPFILE} \
            -u \
    > ${OUTDIR}/union_supported_by_rep${REP}.bed
    echo "  Rep${REP} support: $(wc -l < ${OUTDIR}/union_supported_by_rep${REP}.bed) peaks"
done

# Step 3: Keep peaks found in at least 2 of 3 replicates
# bedtools multiinter counts how many input files overlap each base pair.
# Filtering for $4 >= 2 retains positions covered by ≥2 replicates.
# bedtools merge then reassembles adjacent bases back into intervals.
echo "[$(date)] Applying ≥2/3 reproducibility filter..."

apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
    bedtools multiinter \
        -i ${OUTDIR}/union_supported_by_rep1.bed \
           ${OUTDIR}/union_supported_by_rep2.bed \
           ${OUTDIR}/union_supported_by_rep3.bed \
| awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
| apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
    bedtools merge -i stdin \
> ${OUTDIR}/consensus_peaks_2of3.bed

echo "Consensus peak count (≥2/3 replicates): $(wc -l < ${OUTDIR}/consensus_peaks_2of3.bed)"
echo "[$(date)] Done creating consensus peak set."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/05_consensus_peaks.sh
```

> **Q28:** How many peaks are in the consensus set (present in ≥2 of 3 replicates) compared to the union set? What does this tell you about the reproducibility of your ATAC-seq experiment?

---

## 5.3 — FRiP Score (Fraction of Reads in Peaks)

The FRiP score is a widely used quality metric for ATAC-seq. It measures what fraction of all mapped reads fall within called peaks. FRiP > 0.2 is generally considered acceptable for ATAC-seq; higher is better.

```bash
cat > ${TUTORIAL}/scripts/05_frip_score.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=frip_score
#SBATCH --account=PAS3260
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/frip_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SAM_SIF=${TUTORIAL}/containers/samtools.sif
BEDTOOLS_SIF=${TUTORIAL}/containers/bedtools.sif
CONSENSUS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed
OUTFILE=${TUTORIAL}/atacseq/peaks/consensus/frip_scores.txt

echo "Sample | Total_reads | Reads_in_peaks | FRiP" | tee ${OUTFILE}
echo "-------|-------------|----------------|------" | tee -a ${OUTFILE}

for SAMPLE in atac_wt_rep1 atac_wt_rep2 atac_wt_rep3; do
    BAM=${TUTORIAL}/atacseq/aligned/${SAMPLE}.filtered.bam
    echo "[$(date)] Processing: ${SAMPLE}"

    # Total properly paired reads in the filtered BAM
    # -f 2: properly paired
    # -F 4: exclude unmapped
    TOTAL=$(apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
        samtools view -c -f 2 -F 4 ${BAM})

    # Reads overlapping consensus peaks by any amount (≥1 bp)
    # bedtools intersect -u: report each -a read once if it overlaps any peak
    # samtools view -c applies the same -f 2 -F 4 filter for consistency with TOTAL
    IN_PEAKS=$(apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
        bedtools intersect \
            -a ${BAM} \
            -b ${CONSENSUS} \
            -u \
    | apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
        samtools view -c -f 2 -F 4)

    FRIP=$(echo "scale=4; ${IN_PEAKS} / ${TOTAL}" | bc)

    echo "${SAMPLE} | ${TOTAL} | ${IN_PEAKS} | ${FRIP}" | tee -a ${OUTFILE}
done

echo ""
echo "FRiP scores written to: ${OUTFILE}"
echo "Note: FRiP > 0.2 is generally considered acceptable for ATAC-seq."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/05_frip_score.sh
```

> **Q29:** Report the FRiP score for each of your three replicates. Do they meet the common quality threshold of FRiP > 0.2? What biological and technical factors affect the FRiP score in a fungal ATAC-seq experiment?

---

## 5.4 — Annotate Peaks Relative to Gene Features

To understand where open chromatin occurs in the genome, annotate peaks with respect to annotated features. We use bedtools to intersect peaks with the GFF3 annotation:

```bash
cat > ${TUTORIAL}/scripts/05_annotate_peaks.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=annotate_peaks
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/annotate_peaks_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
BEDTOOLS_SIF=${TUTORIAL}/containers/bedtools.sif
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
PEAKS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed
OUTDIR=${TUTORIAL}/atacseq/peaks/consensus

# Extract promoter regions (2 kb upstream of each gene TSS)
# GFF3 format uses ID= (not gene_id=) in the attributes column.
# Using plain ERE ID=([^;]+) for gawk compatibility — no PCRE groups.
# This extracts gene-AMS68_XXXXXX IDs consistent with the rest of the pipeline.
echo "[$(date)] Building gene promoter BED..."

grep -P "\tgene\t" ${GFF} \
    | awk 'BEGIN{OFS="\t"} {
        match($9, /ID=([^;]+)/, arr)
        gid = arr[1]
        if ($7=="+") print $1, ($4>2000 ? $4-2001 : 0), $4, gid, ".", $7
        else          print $1, $5-1, $5+2000, gid, ".", $7
    }' \
    | sort -k1,1 -k2,2n \
    > ${OUTDIR}/gene_promoters_2kb.bed

echo "Gene promoter regions: $(wc -l < ${OUTDIR}/gene_promoters_2kb.bed)"
echo "Sample gene IDs:"
head -3 ${OUTDIR}/gene_promoters_2kb.bed | cut -f4

# Intersect consensus peaks with promoters
# -wa -wb: output all columns from both -a (peaks) and -b (promoters)
# Output columns: chr_peak, start_peak, end_peak,
#                 chr_prom, start_prom, end_prom, gene_id, ".", strand
echo "[$(date)] Intersecting peaks with promoters..."

apptainer exec --bind ${TUTORIAL} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${PEAKS} \
        -b ${OUTDIR}/gene_promoters_2kb.bed \
        -wa -wb \
    > ${OUTDIR}/peaks_in_promoters.bed

echo "Peaks overlapping promoters: $(wc -l < ${OUTDIR}/peaks_in_promoters.bed)"
echo "Total consensus peaks: $(wc -l < ${PEAKS})"
echo ""
echo "Promoter peaks (first 10):"
head -10 ${OUTDIR}/peaks_in_promoters.bed

# Extract gene IDs near peaks (column 7 of the intersect output)
# Using -F'\t' to split on tabs, ensuring empty fields are not collapsed
awk -F'\t' '{print $7}' ${OUTDIR}/peaks_in_promoters.bed \
    | grep -v '^$' \
    | sort -u \
    > ${OUTDIR}/genes_with_promoter_peaks.txt

echo ""
echo "Unique genes with promoter ATAC peaks: $(wc -l < ${OUTDIR}/genes_with_promoter_peaks.txt)"
echo "Sample gene IDs:"
head -3 ${OUTDIR}/genes_with_promoter_peaks.txt

echo "[$(date)] Done."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/05_annotate_peaks.sh
```

> **Q30:** What percentage of your consensus ATAC-seq peaks overlap with annotated promoter regions (within 2 kb of a TSS)? Where are the remaining peaks located (e.g., within gene bodies, intergenic regions)? Does this distribution make sense given that ATAC-seq primarily captures regulatory chromatin?

---

## 5.5 — Generate Normalized Coverage Tracks with deepTools

Generate bigWig files with bamCoverage
bamCoverage converts each filtered BAM into a normalized bigWig file. We use RPGC normalization (reads per genomic content), which scales each sample to an equivalent sequencing depth — making replicates directly comparable in IGV.
Key parameters for ATAC-seq:

--effectiveGenomeSize 17000000 — approximate mappable size of the P. fructicola genome
--normalizeUsing RPGC — depth normalization
--binSize 10 — 10 bp resolution, appropriate for a compact ~17 Mb genome
--extendReads — extends fragments to their actual insert size rather than just plotting read ends (important for ATAC-seq where the signal comes from the full fragment)
--ignoreDuplicates — skip any remaining duplicates (belt-and-suspenders after markdup)

```bash
cat > ${TUTORIAL}/scripts/05_bamcoverage.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-3
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/bamcoverage_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/deeptools.sif

SAMPLES=(atac_wt_rep1 atac_wt_rep2 atac_wt_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Generating bigWig: ${SAMPLE}"

apptainer exec --bind ${TUTORIAL} ${SIF} \
    bamCoverage \
        --bam ${TUTORIAL}/atacseq/aligned/${SAMPLE}.filtered.bam \
        --outFileName ${TUTORIAL}/atacseq/aligned/${SAMPLE}.rpgc.bw \
        --outFileFormat bigwig \
        --effectiveGenomeSize 17000000 \
        --normalizeUsing RPGC \
        --binSize 10 \
        --extendReads \
        --ignoreDuplicates \
        --numberOfProcessors 8 \
        2>&1

echo "[$(date)] Done: ${SAMPLE}"
ls -lh ${TUTORIAL}/atacseq/aligned/${SAMPLE}.rpgc.bw
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/05_bamcoverage.sh
```

After all three array tasks finish, confirm the output:
```bash
ls -lh ${TUTORIAL}/atacseq/aligned/*.rpgc.bw
# Expect three files, each typically 2–10 MB for a compact fungal genome
```

## 5.6 — Visualizing Data in IGV: Step-by-Step Guide

Download IGV from https://igv.org/doc/desktop/#DownloadPage/ if not already installed. The current version is 2.19.7; these instructions apply to any 2.x release.

### Download files for IGV

You need five things on your local machine. Run these commands in your **local terminal** (not on OSC):

```bash
# Create a local directory for this
mkdir -p ~/Desktop/peltaster_igv
cd ~/Desktop/peltaster_igv

# 1. Reference genome FASTA + index
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/reference/peltaster_fructicola_genome.fa .
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/reference/peltaster_fructicola_genome.fa.fai .

# 2. Gene annotation
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/reference/peltaster_fructicola_annotation.gff3 .

# 3. ATAC-seq bigWig tracks (all three replicates)
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq/aligned/atac_wt_rep1.rpgc.bw \
    <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq/aligned/atac_wt_rep2.rpgc.bw \
    <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq/aligned/atac_wt_rep3.rpgc.bw \
    .

# 4. Consensus ATAC peaks
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/atacseq/peaks/consensus/consensus_peaks_2of3.bed .

# 5. RNA-seq BAMs + indexes (WT rep1 and gh31del rep1 for comparison)
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/aligned/rnaseq_wt_rep1.sorted.bam .
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/aligned/rnaseq_wt_rep1.sorted.bam.bai .
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/aligned/rnaseq_gh31del_rep1.sorted.bam .
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/aligned/rnaseq_gh31del_rep1.sorted.bam.bai .
```

---

### Step 1 — Load the reference genome

IGV needs the *P. fructicola* genome loaded as a custom reference because it is not in IGV's built-in genome list.

1. Open IGV
2. Click **Genomes** → **Load Genome from File...**
3. Navigate to `~/Desktop/peltaster_igv/` and select `peltaster_fructicola_genome.fa`
4. IGV will detect the `.fai` index automatically and load all five chromosomes (CP051139.1 through CP051143.1)
5. The chromosome selector dropdown (top left) will now show the five contigs

### Step 2 — Load the gene annotation

1. Click **File** → **Load from File...**
2. Select `peltaster_fructicola_annotation.gff3`
3. A gene track appears at the bottom of the panel — genes are shown as arrows indicating strand direction
4. Right-click the track name → **Expanded** to see individual gene models rather than collapsed features

### Step 3 — Load the ATAC-seq bigWig tracks

1. Click **File** → **Load from File...**
2. Select all three `*.rpgc.bw` files at once (Shift+click or Cmd/Ctrl+click)
3. Three signal tracks appear — each shows RPGC-normalized read depth across the genome
4. To make replicates visually comparable, right-click each bigWig track → **Set Data Range** → set minimum to `0` and maximum to the same value (e.g., `50`) on all three tracks

> **Tip:** Right-click a bigWig track → **Change Track Color** to assign each replicate a distinct color (e.g., shades of blue).

### Step 4 — Load the consensus peaks BED

1. Click **File** → **Load from File...**
2. Select `consensus_peaks_2of3.bed`
3. This appears as a feature track showing called peak intervals — peaks appear as colored blocks directly below the bigWig signal, letting you visually confirm that peaks coincide with signal enrichment

### Step 5 — Load the RNA-seq BAM files

1. Click **File** → **Load from File...**
2. Select `rnaseq_wt_rep1.sorted.bam` (IGV auto-detects the `.bai` index if it is in the same folder)
3. Repeat for `rnaseq_gh31del_rep1.sorted.bam`
4. Two RNA-seq coverage tracks appear — right-click each → **Set Data Range** → set both to the same maximum (e.g., `200`) so they are directly comparable

### Step 6 — Navigate to the GH31 locus

1. In the search box at the top, type:
CP051143.1:2,680,000-2,690,000
and press Enter
2. This shows a 10 kb window spanning GH31 (AMS68_008039, coordinates 2,682,570–2,685,565) with flanking sequence for context

Your panel should now show (top to bottom):

- Chromosome ruler with scale bar
- Three ATAC-seq bigWig tracks (accessibility signal)
- Consensus peaks BED track (called peak intervals)
- WT RNA-seq BAM track (read coverage)
- *gh31del* RNA-seq BAM track (read coverage)
- Gene annotation track (GH31 gene model)

### Step 7 — Useful navigation tips

| Action | How |
|---|---|
| Zoom in | Click and drag on the ruler, or scroll wheel |
| Zoom out | `Ctrl+Shift+scroll` or the `-` button |
| Jump to a gene | Type gene ID (e.g. `AMS68_000995`) in the search box |
| Expand track height | Click and drag the bottom edge of any track |
| Collapse all tracks | Right-click any track → **Squished** |
| Save a screenshot | **File** → **Save PNG Image...** |
| Save the session | **File** → **Save Session...** — saves track layout as an `.xml` file that can be shared with classmates (requires the same local data files to reload) |

---

> **Q37 (IGV):** At the GH31 locus, describe the spatial relationship between the ATAC-seq signal peak, the called peak interval in the BED track, and the gene model. Is the accessibility signal centered on the promoter, the gene body, or elsewhere? Does the RNA-seq WT track show uniform coverage across the gene body, consistent with active transcription?

> **Q38 (IGV):** Navigate to gene `AMS68_000995` (the top upregulated DEG from Module 03). Compare the ATAC-seq signal at its promoter to what you observed at GH31. Then compare the WT vs. *gh31del* RNA-seq tracks at this locus. Does the visual evidence in IGV match what DESeq2 reported numerically?
---

## Summary

By the end of this module you should have:

- [x] Per-replicate narrowPeak files in `atacseq/peaks/`
- [x] Consensus peak set in `atacseq/peaks/consensus/consensus_peaks_2of3.bed`
- [x] FRiP scores for all three replicates
- [x] Peak annotation overlapping gene promoters

**Proceed to [Module 06 — Multi-Omics Integration](06_multiomics_integration.md)**
