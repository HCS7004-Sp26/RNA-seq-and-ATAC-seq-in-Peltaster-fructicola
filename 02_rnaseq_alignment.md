# Module 02 — RNA-seq Alignment and Read Quantification

## Overview

With clean trimmed reads, you are ready to map them to the *P. fructicola* reference genome and count how many reads support each annotated gene. In this module you will:

1. Build a **HISAT2** splice-aware genome index
2. Align all six RNA-seq samples
3. Sort and index BAM files with **Samtools**
4. Count aligned reads per gene with **featureCounts**
5. Compute TPM values for exploratory analysis

---

## 2.1 — Build the HISAT2 Genome Index

HISAT2 requires a pre-built index of the reference genome. Because *P. fructicola* is a fungus, splicing is limited but present; HISAT2 handles this better than a simple short-read aligner like Bowtie2.

Build the index in a dedicated directory:

```bash
cat > ${TUTORIAL}/scripts/02_hisat2_index.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/hisat2_index_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/hisat2.sif

mkdir -p ${TUTORIAL}/reference/hisat2_index

echo "[$(date)] Building HISAT2 index..."

apptainer exec --bind ${TUTORIAL} ${SIF} \
    hisat2-build \
        -p 8 \
        ${TUTORIAL}/reference/peltaster_fructicola_genome.fa \
        ${TUTORIAL}/reference/hisat2_index/peltaster

echo "[$(date)] Index complete."
ls -lh ${TUTORIAL}/reference/hisat2_index/
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/02_hisat2_index.sh
```

When the job finishes, you should see eight `.ht2` index files in `reference/hisat2_index/`.

> **Q10:** HISAT2 is described as a "splice-aware" aligner. What does this mean, and why is it important even for a compact fungal genome? How does HISAT2 differ from an aligner like BWA-MEM in its treatment of spliced alignments?

---

## 2.2 — Align RNA-seq Reads

Align all six trimmed paired-end samples using a SLURM job array. Each job also immediately sorts and indexes the output BAM to minimize intermediate storage.

```bash
cat > ${TUTORIAL}/scripts/02_hisat2_align_rnaseq.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH --account=PAS3260
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=1-6
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/hisat2_align_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
HISAT2_SIF=${TUTORIAL}/containers/hisat2.sif
SAM_SIF=${TUTORIAL}/containers/samtools.sif
INDEX=${TUTORIAL}/reference/hisat2_index/peltaster
TRIMDIR=${TUTORIAL}/rnaseq/trimmed
OUTDIR=${TUTORIAL}/rnaseq/aligned

SAMPLES=(rnaseq_wt_rep1 rnaseq_wt_rep2 rnaseq_wt_rep3 \
          rnaseq_gh31del_rep1 rnaseq_gh31del_rep2 rnaseq_gh31del_rep3)
SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "[$(date)] Aligning: ${SAMPLE}"

# Align with HISAT2, pipe to samtools sort
apptainer exec --bind ${TUTORIAL} ${HISAT2_SIF} \
    hisat2 \
        -p 8 \
        --dta \
        -x ${INDEX} \
        -1 ${TRIMDIR}/${SAMPLE}_R1_trimmed.fastq.gz \
        -2 ${TRIMDIR}/${SAMPLE}_R2_trimmed.fastq.gz \
        --rg-id ${SAMPLE} \
        --rg "SM:${SAMPLE}" \
        --rg "PL:ILLUMINA" \
        2>${TUTORIAL}/logs/hisat2_${SAMPLE}.summary \
| apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
    samtools sort \
        -@ 8 \
        -m 2G \
        -o ${OUTDIR}/${SAMPLE}.sorted.bam

# Index the sorted BAM
apptainer exec --bind ${TUTORIAL} ${SAM_SIF} \
    samtools index ${OUTDIR}/${SAMPLE}.sorted.bam

echo "[$(date)] Done: ${SAMPLE}"
echo "Alignment summary for ${SAMPLE}:"
cat ${TUTORIAL}/logs/hisat2_${SAMPLE}.summary
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/02_hisat2_align_rnaseq.sh
```

After all jobs complete, check alignment rates from the summary files:

```bash
grep "overall alignment rate" ${TUTORIAL}/logs/hisat2_*.summary
```

> **Q11:** What overall alignment rate did you observe across samples? A rate below ~70% for RNA-seq data mapped to its own genome would be concerning — what possible causes would you investigate?

> **Q12:** The `--dta` flag in HISAT2 stands for "downstream transcriptome assembly." What does it do differently, and is it necessary here given that we already have a reference annotation?

---

## 2.3 — Assess Alignment Quality with Samtools Flagstat

```bash
cat > ${TUTORIAL}/scripts/02_flagstat_rnaseq.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=flagstat_rna
#SBATCH --account=PAS3260
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/flagstat_rna_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/samtools.sif

for BAM in ${TUTORIAL}/rnaseq/aligned/*.sorted.bam; do
    SAMPLE=$(basename ${BAM} .sorted.bam)
    echo "=== ${SAMPLE} ==="
    apptainer exec --bind ${TUTORIAL} ${SIF} \
        samtools flagstat ${BAM}
    echo ""
done > ${TUTORIAL}/rnaseq/aligned/flagstat_summary.txt

cat ${TUTORIAL}/rnaseq/aligned/flagstat_summary.txt
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/02_flagstat_rnaseq.sh
```

---

## 2.4 — Count Reads per Gene with featureCounts

featureCounts (part of the Subread package) assigns aligned reads to annotated gene features. It counts only reads that unambiguously map to a single gene.

Key parameters:
- `-p` paired-end mode
- `-s 2` reverse-strand library (standard for TruSeq stranded RNA-seq)
- `-T 8` number of threads
- `-Q 10` minimum mapping quality
- `-a` annotation file (GFF3)
- `-t exon` feature type to count
- `-g gene_id` attribute to group by

```bash
cat > ${TUTORIAL}/scripts/02_featurecounts.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/featurecounts_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/subread.sif

mkdir -p ${TUTORIAL}/rnaseq/counts

echo "[$(date)] Running featureCounts..."

apptainer exec --bind ${TUTORIAL} ${SIF} \
    featureCounts \
        -p \
        --countReadPairs \
        -s 2 \
        -T 8 \
        -Q 10 \
        -a ${TUTORIAL}/reference/peltaster_fructicola_transcripts.gtf \
        -t exon \
        -g gene_id \
        -o ${TUTORIAL}/rnaseq/counts/peltaster_counts.txt \
        ${TUTORIAL}/rnaseq/aligned/rnaseq_wt_rep1.sorted.bam \
        ${TUTORIAL}/rnaseq/aligned/rnaseq_wt_rep2.sorted.bam \
        ${TUTORIAL}/rnaseq/aligned/rnaseq_wt_rep3.sorted.bam \
        ${TUTORIAL}/rnaseq/aligned/rnaseq_gh31del_rep1.sorted.bam \
        ${TUTORIAL}/rnaseq/aligned/rnaseq_gh31del_rep2.sorted.bam \
        ${TUTORIAL}/rnaseq/aligned/rnaseq_gh31del_rep3.sorted.bam

echo "[$(date)] featureCounts complete."
head -3 ${TUTORIAL}/rnaseq/counts/peltaster_counts.txt
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/02_featurecounts.sh
```

After the job completes, inspect the count matrix:

```bash
# Skip the header comment line, look at structure
head -5 ${TUTORIAL}/rnaseq/counts/peltaster_counts.txt
wc -l ${TUTORIAL}/rnaseq/counts/peltaster_counts.txt   # number of genes + 2 header lines
```

The output file has columns: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, then one column per BAM file.

Also review the summary file:

```bash
cat ${TUTORIAL}/rnaseq/counts/peltaster_counts.txt.summary
```

> **Q13:** What fraction of reads were counted as "Assigned" vs. "Unassigned_NoFeatures"? What does "Unassigned_Ambiguity" mean, and why might a compact fungal genome have higher ambiguity rates than a plant genome with more intergenic space?

---

## 2.5 — Compute TPM Values

Counts are not directly comparable across samples because library sizes differ. TPM (Transcripts Per Million) normalizes by both gene length and total library size. We compute TPM here for exploratory visualization only — **DESeq2 uses raw counts** for statistical testing.

Create a simple TPM script and run it in an interactive session:

```bash
# Start an interactive session — do not run python scripts on the login node
sinteractive -A PAS3260 -t 00:15:00 --mem=4G
```

```bash
cat > ${TUTORIAL}/scripts/compute_tpm.py << 'PYEOF'
#!/usr/bin/env python3
"""
Compute TPM values from featureCounts output.
Usage: python3 compute_tpm.py counts.txt > tpm_matrix.tsv
"""
import sys
import csv

infile = sys.argv[1]
rows = []
samples = []

with open(infile) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        if line.startswith('Geneid'):
            header = line.strip().split('\t')
            # Sample columns start at index 6
            samples = [h.split('/')[-1].replace('.sorted.bam','') for h in header[6:]]
            continue
        parts = line.strip().split('\t')
        gene_id = parts[0]
        length = int(parts[5])
        counts = [float(x) for x in parts[6:]]
        rows.append({'gene': gene_id, 'length': length, 'counts': counts})

# Step 1: RPK (reads per kilobase)
for r in rows:
    r['rpk'] = [c / (r['length'] / 1000.0) for c in r['counts']]

# Step 2: scaling factor per sample
n_samples = len(samples)
scaling = [sum(r['rpk'][i] for r in rows) / 1e6 for i in range(n_samples)]

# Step 3: TPM
print('gene_id\t' + '\t'.join(samples))
for r in rows:
    tpm = [r['rpk'][i] / scaling[i] for i in range(n_samples)]
    print(r['gene'] + '\t' + '\t'.join(f'{v:.4f}' for v in tpm))
PYEOF

python3 ${TUTORIAL}/scripts/compute_tpm.py \
    ${TUTORIAL}/rnaseq/counts/peltaster_counts.txt \
    > ${TUTORIAL}/rnaseq/counts/peltaster_tpm.tsv

head -5 ${TUTORIAL}/rnaseq/counts/peltaster_tpm.tsv
wc -l ${TUTORIAL}/rnaseq/counts/peltaster_tpm.tsv
```

> **Q14:** Look up the TPM value for gene `AMS68_008039` (GH31) in your TPM matrix. Compare WT vs. *gh31del* replicates. What do you observe, and does this make biological sense given the experimental design?

> **Q15:** TPM is commonly used for comparing gene expression levels within a single sample. Why is it **not** appropriate to use TPM values as input to DESeq2 for differential expression testing? What statistical property of raw counts does DESeq2 model?

---

## Summary

By the end of this module you should have:

- [x] HISAT2 genome index in `reference/hisat2_index/`
- [x] Six sorted, indexed BAM files in `rnaseq/aligned/`
- [x] featureCounts output matrix in `rnaseq/counts/peltaster_counts.txt`
- [x] TPM matrix in `rnaseq/counts/peltaster_tpm.tsv`

**Proceed to [Module 03 — Differential Expression Analysis](03_differential_expression.md)**
