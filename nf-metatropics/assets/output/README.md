# Metatropics result folders (detailed)

Paths below are relative to your pipeline **`--outdir`**. The pipeline also creates a Nextflow **`work/`** directory; what follows is only what is published into `outdir`.

---

## Overview

| Group | Role |
|--------|------|
| **`Basecalling/`** | Dorado basecalling and demultiplexing (optional). |
| **`Reads/`** | Read QC, trimming, human / optional host depletion. |
| **`Classification/`** | Virasign viral classification outputs and reports |
| **`Variant_calling/`** | Medaka read alignments and variant calls (VCFs). |
| **`Consensus/`** | iVar draft and Homopolish polished genomes. |
| **`Summary/`** | Final Metatropics report listing all viruses identified (`Summary/virasign/*.html`), plus read-count summaries and pipeline provenance (`Summary/pipeline_info`: software versions + Nextflow reports). |

---

## `Basecalling/` (optional)

| Path | Contents |
|------|----------|
| `Basecalling/basecalling` | Intermediate FASTQ from Dorado before demultiplexing. |
| `Basecalling/demultiplexing` | Per-barcode FASTQ after demultiplexing (trims barcodes and adapters). |

---

## `Reads/`

This stage optionally performs rarefaction, then removes low-quality reads, short reads, and reads that match host background (e.g. human, pan). The host-depleted read set is fed to the rest of the pipeline.

| Path | Contents |
|------|----------|
| `Reads/raw` | Per-sample raw reads after naming / format fixes (`*_fixed.fastq.gz`). |
| `Reads/rarefaction` | Optional: Per-sample FASTQ after optional rarefaction subsampling. |
| `Reads/nanoplot` | Read-length and quality summaries (NanoPlot) on the raw reads. |
| `Reads/fastplong` | Trimmed reads and fastplong QC reports. |
| `Reads/nohuman` | Human-depleted reads (`*_human_depleted.fastq.gz`); mapped reads (`*_human.fastq.gz`). |
| `Reads/nohost` | Host-depleted reads (`*_host_depleted.fastq.gz`); mapped reads (`*_host.fastq.gz`). |

---

## `Classification/`

Host-depleted reads are mapped to the metagenomic database (e.g., Refseq, RVDB) and summarized so you can see which organisms (taxa) are present in each sample.

| Path   | Contents |
|--------|----------|
| `Classification/virasign/` | Virasign per-sample viral classification results. |

| Output | Meaning  |
|--------|----------|
| `*_unfiltered_all_references.json` | All candidate viral hits. |
| `*_final_selected_references.json` | Final confident viral hits. |
| `*.fasta` | Best reference sequence(s) selected by Virasign. |
| `mreads.fastq.gz` | Reads mapped to the reference(s) for each virus. |
| `*.bam` | Read alignments against the best reference(s). |
| `*.bai` | BAM index for the alignments. |
| `*coverage*.pdf` | Coverage plots per virus/reference. |

For more Virasign details, see [`DaanJansen94/virasign`](https://github.com/DaanJansen94/virasign).

---

## `Variant_calling/`

For each candidate virus, **[Medaka](https://github.com/nanoporetech/medaka)** first writes a BAM of reads aligned to the viral reference. It then uses that BAM in haploid variant calling, applying neural-network inference on pileups of those aligned reads (networks trained for Oxford Nanopore basecalled data), and emits a VCF describing differences from that reference.

| Path | Contents |
|------|----------|
| `Variant_calling/reffix` | Reference FASTA with cleaned headers for each virus. |
| `Variant_calling/medaka` | BAMs, VCF variant calls, and Medaka coverage text. |
| `Variant_calling/addingDepth` | Per-virus **`*.sdepth.tsv`** depth tables (coverage + consensus + classification). |

---

## `Consensus/`

Consensus building takes the Medaka-produced BAM and feeds it to iVar, which calls a draft consensus sequence; Homopolish then polishes that draft. 

| Path | Contents |
|------|----------|
| `Consensus/ivar` | Consensus sequences from iVar (input to Homopolish). |
| `Consensus/homopolish` | Polished consensus FASTA (typical final genome per virus per sample). |

---

## `Summary/`

Run-wide summaries and provenance.

| Path | Contents |
|------|----------|
| `Summary/readcount` | Read distribution (`read_distribution.html`). |
| `Summary/virasign` | Final Metatropics HTML report for quick inspection across samples. |
| `Summary/pipeline_info` | Nextflow trace/reports and software versions. |

For a quick overview of viruses found, open **`Summary/virasign/results_summary_*.html`**. 
