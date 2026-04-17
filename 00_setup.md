# Module 00 — Setup and Computing Environment

## Overview

This module covers everything you need before running any analyses: setting up your directory structure, pulling software containers, and confirming that your OSC environment is ready.

All compute jobs in this tutorial run on the **Ohio Supercomputer Center (OSC)** under project allocation **PAS3260**, using **Apptainer** containers to ensure reproducible software environments. Because OSC does not provide most bioinformatics modules centrally, all tools are delivered as pre-built Apptainer images pulled from the Seqera Wave Community Registry.

---

## 0.1 — Log In and Set Up Your Working Directory

Log in to OSC via the OnDemand terminal at https://ondemand.osc.edu, or via SSH:

```bash
ssh <username>@pitzer.osc.edu
```

Create your working directory for this tutorial:

```bash
export TUTORIAL=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq
mkdir -p ${TUTORIAL}/{containers,reference,rnaseq/{fastqc_raw,trimmed,aligned,counts},atacseq/{fastqc_raw,trimmed,aligned,peaks},dge,integration,scripts,logs}
cd ${TUTORIAL}
```

Set a convenience variable pointing to the shared class data:

```bash
export CLASSDATA=/fs/scratch/PAS3260/RNAseq_ATACseq
```

> **Tip:** Add both `export` lines to your `~/.bashrc` so they persist across sessions.

---

## 0.2 — Symlink the Reference and Input Data

Rather than copying the large FASTQ files, create symbolic links:

```bash
# Reference genome and annotation
ln -s ${CLASSDATA}/reference/peltaster_fructicola_genome.fa     ${TUTORIAL}/reference/
ln -s ${CLASSDATA}/reference/peltaster_fructicola_genome.fa.fai ${TUTORIAL}/reference/
ln -s ${CLASSDATA}/reference/peltaster_fructicola_annotation.gff3 ${TUTORIAL}/reference/
ln -s ${CLASSDATA}/reference/peltaster_fructicola_transcripts.gtf ${TUTORIAL}/reference/

# RNA-seq reads
for f in ${CLASSDATA}/rnaseq_illumina/*.fastq.gz; do
    ln -s "$f" ${TUTORIAL}/rnaseq/
done

# ATAC-seq reads
for f in ${CLASSDATA}/atacseq/*.fastq.gz; do
    ln -s "$f" ${TUTORIAL}/atacseq/
done

# Genome sizes file (needed for ATAC-seq)
ln -s ${CLASSDATA}/05_atacseq/genome.sizes ${TUTORIAL}/atacseq/
```

Verify that you have 12 RNA-seq files (6 samples × 2 reads) and 6 ATAC-seq files (3 samples × 2 reads):

```bash
ls ${TUTORIAL}/rnaseq/*.fastq.gz | wc -l   # expect: 12
ls ${TUTORIAL}/atacseq/*.fastq.gz | wc -l  # expect: 6
```

> **Q1:** List the six RNA-seq sample names. What are the two experimental conditions being compared, and how many biological replicates are there per condition?

---

## 0.3 — Pull Software Containers

All tools are run through Apptainer container images. Pull them now so they are ready when your jobs run. This step requires internet access and can take several minutes — submit it as a SLURM job.

```bash
cat > ${TUTORIAL}/scripts/00_pull_containers.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=pull_containers
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/logs/pull_containers_%j.log

set -euo pipefail

CONTAINER_DIR=/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/containers

cd ${CONTAINER_DIR}

echo "[$(date)] Pulling FastQC..."
apptainer pull fastqc.sif \
    oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960

echo "[$(date)] Pulling MultiQC..."
apptainer pull multiqc.sif \
    oras://community.wave.seqera.io/library/multiqc:1.33--e3576ddf588fa00d

echo "[$(date)] Pulling Trimmomatic..."
apptainer pull trimmomatic.sif \
    oras://community.wave.seqera.io/library/trimmomatic:0.40--7b5b7590373e6fc4

echo "[$(date)] Pulling HISAT2..."
apptainer pull hisat2.sif \
    oras://community.wave.seqera.io/library/hisat2:2.2.2--0a747e741adc8dcc

echo "[$(date)] Pulling Samtools..."
apptainer pull samtools.sif \
    oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a

echo "[$(date)] Pulling Subread (featureCounts)..."
apptainer pull subread.sif \
    oras://community.wave.seqera.io/library/subread:2.1.1--bae420bffb4edf16

echo "[$(date)] Pulling Bowtie2..."
apptainer pull bowtie2.sif \
    oras://community.wave.seqera.io/library/bowtie2:2.5.5--21b0835eb76ba3c0

echo "[$(date)] Pulling MACS3..."
apptainer pull macs3.sif \
    oras://community.wave.seqera.io/library/macs3:3.0.4--1d41c250736f1138

echo "[$(date)] All containers pulled successfully."
ls -lh ${CONTAINER_DIR}/*.sif
EOF
cd ${TUTORIAL}
sbatch ${TUTORIAL}/scripts/00_pull_containers.sh
```

> **Tip:** Monitor job progress with `squeue -u OSCuser`. Container SIF files are large (~400 MB–1 GB each); the full pull job may take 30–45 minutes.

After the job completes, confirm that all `.sif` files are present:

```bash
ls -lh ${TUTORIAL}/containers/*.sif
```

You should see nine container images: `fastqc.sif`, `multiqc.sif`, `trimmomatic.sif`, `hisat2.sif`, `samtools.sif`, `subread.sif`, `bowtie2.sif`, and `macs3.sif`.

> **Q2:** What is the purpose of using Apptainer containers instead of loading modules directly with `module load`? Name one advantage this approach offers for reproducibility.

---

## 0.4 — Inspect the Reference Genome

Before any alignment, it is good practice to understand the structure of the reference:

```bash
# How many chromosomes/scaffolds?
grep "^>" ${TUTORIAL}/reference/peltaster_fructicola_genome.fa | wc -l

# Scaffold names and approximate sizes
grep "^>" ${TUTORIAL}/reference/peltaster_fructicola_genome.fa

# Check that the annotation matches the genome
cut -f1 ${TUTORIAL}/reference/peltaster_fructicola_annotation.gff3 | sort -u | grep -v "^#"
```

> **Q3:** How many chromosomes does the *P. fructicola* genome have? Look up the Wang et al. 2020 paper: what makes the genome organization of this species unusual compared to other Dothideomycetes?

---

## 0.5 — Experimental Design Summary

Before running any analysis, you should understand what you are comparing:

| Sample | Condition | Replicate |
|---|---|---|
| rnaseq_wt_rep1 | Wild type | 1 |
| rnaseq_wt_rep2 | Wild type | 2 |
| rnaseq_wt_rep3 | Wild type | 3 |
| rnaseq_gh31del_rep1 | *gh31* deletion | 1 |
| rnaseq_gh31del_rep2 | *gh31* deletion | 2 |
| rnaseq_gh31del_rep3 | *gh31* deletion | 3 |
| atac_wt_rep1 | Wild type | 1 |
| atac_wt_rep2 | Wild type | 2 |
| atac_wt_rep3 | Wild type | 3 |

The ATAC-seq data are from the **wild type only**. This represents the chromatin accessibility landscape of *P. fructicola* under normal conditions — and we will later ask whether genes that change expression in the *gh31* deletion mutant are located near accessible chromatin regions.

> **Q4:** What biological hypothesis does comparing RNA-seq between WT and *gh31del* allow you to test? State this as a specific scientific question.

---

## Summary

By the end of this module you should have:

- [x] A working directory with the full subdirectory structure
- [x] Symbolic links to all FASTQ files and reference data
- [x] All nine Apptainer container images pulled and confirmed present
- [x] Familiarity with the experimental design and reference genome structure

**Proceed to [Module 01 — RNA-seq QC & Trimming](01_rnaseq_qc.md)**
