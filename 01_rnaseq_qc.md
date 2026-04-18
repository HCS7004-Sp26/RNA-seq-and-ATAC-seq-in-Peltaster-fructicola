# Module 01 — RNA-seq Quality Control and Adapter Trimming

## Overview

Before any alignment, you must assess the quality of your raw reads and remove adapter sequences and low-quality bases. In this module you will:

1. Run **FastQC** on all 12 RNA-seq FASTQ files
2. Summarize results with **MultiQC**
3. Trim adapters and low-quality bases with **Trimmomatic**
4. Re-run FastQC and MultiQC to confirm improvement

---

## 1.1 — FastQC on Raw Reads

FastQC computes per-read quality scores (Phred), GC content distribution, adapter contamination, and sequence duplication levels. Run it on all 12 RNA-seq files in a single SLURM job using a task array.

```bash
cat > ${TUTORIAL}/scripts/01_fastqc_raw_rnaseq.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=fastqc_raw_rna
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/fastqc_raw_rna_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/fastqc.sif

echo "[$(date)] Running FastQC on raw RNA-seq reads..."

apptainer exec --bind ${TUTORIAL} ${SIF} \
    fastqc \
        --outdir ${TUTORIAL}/rnaseq/fastqc_raw \
        --threads 4 \
        ${TUTORIAL}/rnaseq/rnaseq_wt_rep1_R1.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_wt_rep1_R2.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_wt_rep2_R1.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_wt_rep2_R2.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_wt_rep3_R1.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_wt_rep3_R2.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_gh31del_rep1_R1.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_gh31del_rep1_R2.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_gh31del_rep2_R1.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_gh31del_rep2_R2.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_gh31del_rep3_R1.fastq.gz \
        ${TUTORIAL}/rnaseq/rnaseq_gh31del_rep3_R2.fastq.gz

echo "[$(date)] FastQC complete."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/01_fastqc_raw_rnaseq.sh
```

After the job completes, inspect the output directory:

```bash
ls ${TUTORIAL}/rnaseq/fastqc_raw/
```

You will see one `_fastqc.html` and one `_fastqc.zip` per input file.

---

## 1.2 — Summarize with MultiQC

Rather than opening 12 individual FastQC reports, MultiQC aggregates them into a single interactive HTML:

```bash
cat > ${TUTORIAL}/scripts/01_multiqc_raw_rnaseq.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=multiqc_raw_rna
#SBATCH --account=PAS3260
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/multiqc_raw_rna_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/multiqc.sif

apptainer exec --bind ${TUTORIAL} ${SIF} \
    multiqc \
        ${TUTORIAL}/rnaseq/fastqc_raw/ \
        --outdir ${TUTORIAL}/rnaseq/fastqc_raw/ \
        --filename multiqc_rnaseq_raw.html \
        --title "RNA-seq Raw Reads QC"

echo "[$(date)] MultiQC complete."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/01_multiqc_raw_rnaseq.sh
```

Download the resulting `multiqc_rnaseq_raw.html` to your local machine and open it in a browser:

```bash
# Run this on your LOCAL terminal (not OSC)
scp <username>@pitzer.osc.edu:/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/rnaseq/fastqc_raw/multiqc_rnaseq_raw.html ~/Desktop/
```

> **Q5:** In the MultiQC report, what is the average per-base sequence quality (Phred score) across all samples at the 3' end of reads? Are there any samples that look noticeably worse than the others?

> **Q6:** Does the "Adapter Content" module flag any adapter contamination? What adapter sequences are being detected, and why must they be removed before alignment?

> **Q7:** The "Sequence Duplication Levels" module shows the percentage of duplicate reads. What causes high duplication in RNA-seq data, and is a moderate level of duplication necessarily a problem for differential expression analysis?

---

## 1.3 — Adapter Trimming with Trimmomatic

Trimmomatic performs four operations in this workflow:

- **ILLUMINACLIP**: Remove Illumina TruSeq adapter sequences
- **LEADING/TRAILING**: Remove low-quality bases (Phred < 3) from read ends
- **SLIDINGWINDOW**: Trim when a 4-base sliding window drops below Phred 15
- **MINLEN**: Discard reads shorter than 36 bp after trimming

Submit trimming for all six samples using a SLURM job array (array index 1–6 maps to each sample):

```bash
cat > ${TUTORIAL}/scripts/01_trimmomatic_rnaseq.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=trim_rnaseq
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-6
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/trim_rnaseq_%A_%a.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
SIF=${TUTORIAL}/containers/trimmomatic.sif
INDIR=${TUTORIAL}/rnaseq
OUTDIR=${TUTORIAL}/rnaseq/trimmed

# Map SLURM array index to sample name
SAMPLES=(rnaseq_wt_rep1 rnaseq_wt_rep2 rnaseq_wt_rep3 \
          rnaseq_gh31del_rep1 rnaseq_gh31del_rep2 rnaseq_gh31del_rep3)

SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "[$(date)] Trimming sample: ${SAMPLE}"

apptainer exec --bind ${TUTORIAL} ${SIF} \
    trimmomatic PE \
        -threads 8 \
        ${INDIR}/${SAMPLE}_R1.fastq.gz \
        ${INDIR}/${SAMPLE}_R2.fastq.gz \
        ${OUTDIR}/${SAMPLE}_R1_trimmed.fastq.gz \
        ${OUTDIR}/${SAMPLE}_R1_unpaired.fastq.gz \
        ${OUTDIR}/${SAMPLE}_R2_trimmed.fastq.gz \
        ${OUTDIR}/${SAMPLE}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:/opt/conda/share/trimmomatic-0.40-0/adapters/TruSeq3-PE-2.fa:2:30:10 \
        -phred33 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

echo "[$(date)] Done: ${SAMPLE}"
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/01_trimmomatic_rnaseq.sh
```

After all six array jobs complete, check the trimmed output:

```bash
ls ${TUTORIAL}/rnaseq/trimmed/*_trimmed.fastq.gz | wc -l   # expect: 12
```

The `_unpaired.fastq.gz` files contain reads whose mate was discarded. These are not used in subsequent steps.

---

## 1.4 — Post-Trimming QC

Re-run FastQC and MultiQC on the trimmed reads to confirm that adapter sequences are gone and overall quality has improved:

```bash
cat > ${TUTORIAL}/scripts/01_fastqc_trimmed_rnaseq.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=fastqc_trim_rna
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/fastqc_trim_rna_%j.log

set -euo pipefail

TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
mkdir -p ${TUTORIAL}/rnaseq/fastqc_trimmed

apptainer exec --bind ${TUTORIAL} ${TUTORIAL}/containers/fastqc.sif \
    fastqc \
        --outdir ${TUTORIAL}/rnaseq/fastqc_trimmed \
        --threads 4 \
        ${TUTORIAL}/rnaseq/trimmed/*_trimmed.fastq.gz

apptainer exec --bind ${TUTORIAL} ${TUTORIAL}/containers/multiqc.sif \
    multiqc \
        ${TUTORIAL}/rnaseq/fastqc_trimmed/ \
        --outdir ${TUTORIAL}/rnaseq/fastqc_trimmed/ \
        --filename multiqc_rnaseq_trimmed.html \
        --title "RNA-seq Trimmed Reads QC"

echo "[$(date)] Post-trimming QC complete."
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/01_fastqc_trimmed_rnaseq.sh
```

> **Q8:** Compare the pre- and post-trimming MultiQC reports. What changed in the Adapter Content and Per Base Sequence Quality panels? Approximately what percentage of reads were discarded or quality-trimmed across your samples?

> **Q9:** Why is it important to use the same trimming parameters across all samples in a differential expression experiment?

---

## Summary

By the end of this module you should have:

- [x] Raw FastQC/MultiQC reports for all 12 RNA-seq FASTQ files
- [x] 12 trimmed FASTQ files (paired) in `rnaseq/trimmed/`
- [x] Post-trimming MultiQC report confirming adapter removal

**Proceed to [Module 02 — RNA-seq Alignment & Quantification](02_rnaseq_alignment.md)**