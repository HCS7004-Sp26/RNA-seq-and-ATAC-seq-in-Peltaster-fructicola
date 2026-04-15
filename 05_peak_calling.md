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
cat > ${TUTORIAL}/05_macs3_individual.sh << 'EOF'
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

sbatch ${TUTORIAL}/05_macs3_individual.sh
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
cat > ${TUTORIAL}/05_consensus_peaks.sh << 'EOF'
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
BEDTOOLS_SIF=${CLASSDATA}/containers/bedtools.sif
PEAKDIR=${TUTORIAL}/atacseq/peaks
OUTDIR=${TUTORIAL}/atacseq/peaks/consensus
SIZES=${TUTORIAL}/atacseq/genome.sizes

mkdir -p ${OUTDIR}

REP1=${PEAKDIR}/atac_wt_rep1/atac_wt_rep1_peaks.narrowPeak
REP2=${PEAKDIR}/atac_wt_rep2/atac_wt_rep2_peaks.narrowPeak
REP3=${PEAKDIR}/atac_wt_rep3/atac_wt_rep3_peaks.narrowPeak

# Step 1: Pool all peaks and merge overlapping ones into a union set
cat ${REP1} ${REP2} ${REP3} \
    | sort -k1,1 -k2,2n \
    | apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
        bedtools merge -i stdin \
    > ${OUTDIR}/union_peaks.bed

echo "Union peak count: $(wc -l < ${OUTDIR}/union_peaks.bed)"

# Step 2: Count how many replicates support each union peak
# A union peak is "reproducible" if ≥2 replicates overlap it

for REP in 1 2 3; do
    REPFILE=${PEAKDIR}/atac_wt_rep${REP}/atac_wt_rep${REP}_peaks.narrowPeak
    apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
        bedtools intersect \
            -a ${OUTDIR}/union_peaks.bed \
            -b ${REPFILE} \
            -u \
    > ${OUTDIR}/union_supported_by_rep${REP}.bed
done

# Keep peaks found in at least 2 of 3 replicates
apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
    bedtools intersect \
        -a ${OUTDIR}/union_peaks.bed \
        -b ${OUTDIR}/union_supported_by_rep1.bed \
           ${OUTDIR}/union_supported_by_rep2.bed \
           ${OUTDIR}/union_supported_by_rep3.bed \
        -u \
        -filenames \
| sort -k1,1 -k2,2n \
| apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
    bedtools merge -i stdin \
> ${OUTDIR}/consensus_peaks_2of3.bed

echo "Consensus peak count (≥2/3 replicates): $(wc -l < ${OUTDIR}/consensus_peaks_2of3.bed)"

# Also create a merged all-replicates BAM for visualization
echo "[$(date)] Done creating consensus peak set."
EOF

sbatch ${TUTORIAL}/05_consensus_peaks.sh
```

> **Q28:** How many peaks are in the consensus set (present in ≥2 of 3 replicates) compared to the union set? What does this tell you about the reproducibility of your ATAC-seq experiment?

---

## 5.3 — FRiP Score (Fraction of Reads in Peaks)

The FRiP score is a widely used quality metric for ATAC-seq. It measures what fraction of all mapped reads fall within called peaks. FRiP > 0.2 is generally considered acceptable for ATAC-seq; higher is better.

```bash
cat > ${TUTORIAL}/05_frip_score.sh << 'EOF'
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
BEDTOOLS_SIF=${CLASSDATA}/containers/bedtools.sif
CONSENSUS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed

echo "Sample | Total_reads | Reads_in_peaks | FRiP"
echo "-------|-------------|----------------|------"

for SAMPLE in atac_wt_rep1 atac_wt_rep2 atac_wt_rep3; do
    BAM=${TUTORIAL}/atacseq/aligned/${SAMPLE}.filtered.bam

    TOTAL=$(apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
        samtools view -c -f 2 ${BAM})

    IN_PEAKS=$(apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
        bedtools intersect -a ${BAM} -b ${CONSENSUS} -u -f 0.5 \
    | apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
        samtools view -c)

    FRIP=$(echo "scale=4; ${IN_PEAKS} / ${TOTAL}" | bc)
    echo "${SAMPLE} | ${TOTAL} | ${IN_PEAKS} | ${FRIP}"
done
EOF

sbatch ${TUTORIAL}/05_frip_score.sh
```

> **Q29:** Report the FRiP score for each of your three replicates. Do they meet the common quality threshold of FRiP > 0.2? What biological and technical factors affect the FRiP score in a fungal ATAC-seq experiment?

---

## 5.4 — Annotate Peaks Relative to Gene Features

To understand where open chromatin occurs in the genome, annotate peaks with respect to annotated features. We use bedtools to intersect peaks with the GFF3 annotation:

```bash
cat > ${TUTORIAL}/05_annotate_peaks.sh << 'EOF'
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
BEDTOOLS_SIF=${CLASSDATA}/containers/bedtools.sif
GFF=${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3
PEAKS=${TUTORIAL}/atacseq/peaks/consensus/consensus_peaks_2of3.bed
OUTDIR=${TUTORIAL}/atacseq/peaks/consensus

# Extract promoter regions (2 kb upstream of each gene TSS)
# First get gene features from GFF3
grep -P "\tgene\t" ${GFF} \
    | awk 'BEGIN{OFS="\t"} {
        match($9, /gene_id=([^;]+)/, arr)
        # For + strand: promoter is [start-2000, start]
        # For - strand: promoter is [end, end+2000]
        if ($7=="+") print $1, ($4>2000 ? $4-2001 : 0), $4, arr[1], ".", $7
        else print $1, $5-1, $5+2000, arr[1], ".", $7
    }' \
    | sort -k1,1 -k2,2n \
    > ${OUTDIR}/gene_promoters_2kb.bed

# Intersect consensus peaks with promoters
apptainer exec --bind ${TUTORIAL},${CLASSDATA} ${BEDTOOLS_SIF} \
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

# Extract just gene IDs near peaks
awk '{print $7}' ${OUTDIR}/peaks_in_promoters.bed | sort -u \
    > ${OUTDIR}/genes_with_promoter_peaks.txt
echo ""
echo "Unique genes with promoter ATAC peaks: $(wc -l < ${OUTDIR}/genes_with_promoter_peaks.txt)"
EOF

sbatch ${TUTORIAL}/05_annotate_peaks.sh
```

> **Q30:** What percentage of your consensus ATAC-seq peaks overlap with annotated promoter regions (within 2 kb of a TSS)? Where are the remaining peaks located (e.g., within gene bodies, intergenic regions)? Does this distribution make sense given that ATAC-seq primarily captures regulatory chromatin?

---

## 5.5 — Visualize Peak Tracks (Optional)

For visual inspection of peaks in a genome browser, create BigWig files from the filtered BAMs. These can be loaded into IGV or the UCSC Genome Browser.

```bash
# This step requires deepTools — check if available in a container,
# or use samtools to create a per-base coverage BED instead:
for SAMPLE in atac_wt_rep1 atac_wt_rep2 atac_wt_rep3; do
    apptainer exec --bind ${TUTORIAL} ${TUTORIAL}/containers/samtools.sif \
        samtools depth \
            -a \
            ${TUTORIAL}/atacseq/aligned/${SAMPLE}.filtered.bam \
    | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' \
    | gzip \
    > ${TUTORIAL}/atacseq/aligned/${SAMPLE}.coverage.bedgraph.gz &
done
wait
echo "Coverage tracks created."
```

---

## Summary

By the end of this module you should have:

- [x] Per-replicate narrowPeak files in `atacseq/peaks/`
- [x] Consensus peak set in `atacseq/peaks/consensus/consensus_peaks_2of3.bed`
- [x] FRiP scores for all three replicates
- [x] Peak annotation overlapping gene promoters

**Proceed to [Module 06 — Multi-Omics Integration](06_multiomics_integration.md)**
