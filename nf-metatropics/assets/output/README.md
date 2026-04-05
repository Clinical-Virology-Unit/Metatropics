# Metatropics result folders (detailed)

Paths below are relative to your pipeline **`--outdir`**. The pipeline also creates a Nextflow **`work/`** directory (scratch); what follows is only what is **published** into `outdir`.

---

## Overview

| Group | Role |
|--------|------|
| **`Basecalling/`** | POD5: Dorado basecalling and demultiplexing (optional). |
| **`Reads/`** | Read QC, trimming, human / optional host depletion. |
| **`Classification/`** | MetaMaps taxonomic outputs. |
| **`Viral_reads/`** | Per-virus extracted FASTQs. |
| **`Variant_calling/`** | Viral references, Medaka alignments/variants, merged depth tables. |
| **`Consensus/`** | iVar draft and Homopolish polished genomes. |
| **`Summary/`** | Final tables, read-count figure, coverage PDFs, software versions and Nextflow trace/reports. |

---

## `Basecalling/` (optional)

| Path | Contents |
|------|----------|
| `Basecalling/basecalling` | Intermediate FASTQ from Dorado before demultiplexing. |
| `Basecalling/demultiplexing` | Per-barcode FASTQ after demultiplexing. |

---

## `Reads/`

This stage removes low-quality reads and reads that match host background (e.g. human). The host-depleted read set feeds the rest of the pipeline. Alignments to host genomes are **not** copied to `outdir` (they stay in Nextflow `work/`); the table below lists the **published** FASTQs and QC reports.

| Path | Contents |
|------|----------|
| `Reads/fix` | Per-sample FASTQ after naming / format fixes (compressed). |
| `Reads/rarefaction` | Per-sample FASTQ after optional rarefaction subsampling. |
| `Reads/fastp` | Trimmed reads and FASTP reports. |
| `Reads/nanoplot` | Read-length and quality summaries (NanoPlot) on the **fixed** per-sample FASTQs (same inputs as `Reads/fix`), **before** fastp trimming. |
| `Reads/nohuman` | Human-depleted reads or reads not mapping to the human reference (`*_other.fastq.gz`). |
| `Reads/nohost` | Optional: FASTQ after extra host depletion (`*_other.fastq.gz`). |

---

## `Classification/`

Host-depleted reads are mapped to the metagenomic database and summarized so you can see **which organisms (taxa) are present** in each sample.

| Path | Contents |
|------|----------|
| `Classification/metamaps` | MetaMaps mapping and classification (`mapDirectly`, `Classify`). |

---

## `Viral_reads/`

Reads assigned to a given **virus** are extracted into FASTQs **per sample and per virus**.

| Path | Contents |
|------|----------|
| `Viral_reads/seqtk` | FASTQ per sample per virus. |

---

## `Variant_calling/`

For each candidate virus, reads are aligned to the **viral reference** (**BAMs**). **Variant calls** are differences from that reference.

**Medaka** ([Oxford Nanopore documentation](https://nanoporetech.github.io/medaka/index.html)) uses **neural-network inference on pileups** of reads aligned to the reference (not raw signal). This pipeline runs **haploid variant calling** (`medaka_haploid_variant`), producing **BAM** and **VCF**. 

| Path | Contents |
|------|----------|
| `Variant_calling/reffix` | Reference FASTA with **cleaned headers** for each virus. |
| `Variant_calling/medaka` | **BAMs**, **VCF** variant calls, and Medaka coverage text. |
| `Variant_calling/addingDepth` | Per-virus **`*.sdepth.tsv`** depth tables (coverage + consensus + classification). |

---

## `Consensus/`

A **consensus genome** is called per virus (iVar draft, **Homopolish** final genome).

| Path | Contents |
|------|----------|
| `Consensus/ivar` | Consensus sequences from iVar (input to Homopolish). |
| `Consensus/homopolish` | **Polished** consensus FASTA (typical final genome per virus per sample). |

---

## `Summary/`

Run-wide summaries and provenance.

| Path | Contents |
|------|----------|
| `Summary/final` | Combined **final TSV** across the whole run. |
| `Summary/rcoverage` | Coverage **PDFs** for identified viruses (if enabled). |
| `Summary/read_count` | Read-count **CSV/PDF** and staged inputs for the read-distribution figure (copies from **`Reads/`** and scans **`Classification/metamaps`**). |
| `Summary/pipeline_info` | Nextflow trace/reports and software versions. |

**Typical reporting priorities:** **`Consensus/homopolish`**, **`Summary/final`**, **`Summary/rcoverage`**, and **`Summary/read_count`**.
