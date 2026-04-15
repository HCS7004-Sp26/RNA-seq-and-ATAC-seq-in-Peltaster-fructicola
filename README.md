# Multi-Omics Analysis: RNA-seq and ATAC-seq in *Peltaster fructicola*

## HCS 7004 — Genome Analytics | Ohio State University

---

## Biological Context

*Peltaster fructicola* is a Dothideomycete fungal pathogen responsible for sooty blotch disease on apple fruit. Despite having one of the smallest known filamentous fungal genomes (~17 Mb, ~8,100 genes), it encodes a suite of carbohydrate-active enzymes (CAZymes) critical for colonizing host plant tissue. Among these, **GH31** (gene AMS68_008039; Family 31 Glycoside Hydrolase) has been proposed as a key enzyme involved in degrading host cell wall glycans during infection.

In this tutorial you will analyze two complementary genomic datasets generated from *P. fructicola* cultures:

| Experiment | Conditions | Replicates |
|---|---|---|
| RNA-seq (Illumina paired-end) | Wild type (WT) vs. *gh31* deletion mutant | 3 per condition |
| ATAC-seq (Illumina paired-end) | Wild type (WT) only | 3 replicates |

Together, these data allow you to ask: **When GH31 is deleted, which genes change in expression, and are those genes located in regions of open chromatin?**

---

## Learning Objectives

By the end of this tutorial, you will be able to:

1. Perform quality control and adapter trimming on short-read RNA-seq and ATAC-seq data
2. Build a genome index and align short reads with HISAT2 (RNA-seq) and Bowtie2 (ATAC-seq)
3. Quantify gene expression using featureCounts and compute TPM values
4. Identify differentially expressed genes (DEGs) using DESeq2
5. Call chromatin accessibility peaks from ATAC-seq data using MACS3
6. Integrate RNA-seq and ATAC-seq results to identify candidate regulatory loci
7. Interpret multi-omics findings in the context of fungal pathogen biology

---

## Tutorial Structure

| Module | Topic |
|---|---|
| [00 — Setup](00_setup.md) | Computing environment, data organization, software containers |
| [01 — RNA-seq QC & Trimming](01_rnaseq_qc.md) | FastQC, MultiQC, Trimmomatic |
| [02 — RNA-seq Alignment & Quantification](02_rnaseq_alignment.md) | HISAT2 index, alignment, featureCounts |
| [03 — Differential Expression Analysis](03_differential_expression.md) | DESeq2, volcano plots, DEG interpretation |
| [04 — ATAC-seq Processing](04_atacseq_processing.md) | FastQC, Trimmomatic, Bowtie2, filtering |
| [05 — Peak Calling](05_peak_calling.md) | MACS3 peak calling, peak annotation |
| [06 — Multi-Omics Integration](06_multiomics_integration.md) | bedtools, correlation analysis, biological synthesis |

---

## Data Organization

All input data are located under the shared class directory on OSC:

```
/fs/scratch/PAS3260/RNAseq_ATACseq/
├── reference/
│   ├── peltaster_fructicola_genome.fa
│   ├── peltaster_fructicola_genome.fa.fai
│   ├── peltaster_fructicola_annotation.gff3
│   └── peltaster_fructicola_proteins.faa
├── rnaseq_illumina/
│   ├── rnaseq_wt_rep[1-3]_R[1-2].fastq.gz
│   └── rnaseq_gh31del_rep[1-3]_R[1-2].fastq.gz
└── atacseq/
    ├── atac_wt_rep[1-3]_R[1-2].fastq.gz
    └── genome.sizes
```

Your working directory for all output will be:

```
/fs/scratch/PAS3260/<your_username>/rnaseq_atacseq/
```

---

## Reference

Wang et al. 2020. A chromosome-scale assembly of the smallest Dothideomycete genome reveals a unique genome compaction mechanism in filamentous fungi. *BMC Genomics* **21**:321. https://doi.org/10.1186/s12864-020-6732-8

---

> **Note for instructors:** INSTRUCTOR_ prefixed files (truth tables and expected peak BED) are located in the class data directory but are not exposed to students. These are for grading and validation only.
