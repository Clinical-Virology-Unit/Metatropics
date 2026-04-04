[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metatropics-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/metatropics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# Metatropics: A Viral Metagenomics ONT Pipeline

**Metatropics** is a Nextflow-driven bioinformatics pipeline for metagenomic Oxford Nanopore sequencing data. It is built to detect viral pathogens in complex samples (e.g. blood and swabs from a range of body sites) and, where coverage allows, to generate high-quality viral consensus genomes and to perform variant analysis.

**Metatropics** is an abbreviation of **Metagenomics for Tropical Fevers**, and reflects how the project began, with an emphasis on finding human viral pathogens in patients presenting with tropical fevers. The same pipeline has since been validated and applied outside that first setting, including for other febrile syndromes, for genomic surveillance, and for research and diagnostic questions around viral pathogens relevant to human health.

## Pipeline summary

![Figure](./nf-metatropics/assets/logo/Metatropics.jpg)

## 1. Clone the repository
```bash
git clone https://github.com/DaanJansen94/Metatropics.git
cd Metatropics
```

## 2. Java and Nextflow
You need **Java 17+** and **[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) ≥ 22.10.1**. On Debian/Ubuntu you can do:
```bash
sudo apt update && sudo apt install -y openjdk-17-jdk curl
curl -sSL https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/
nextflow -version
```

## 3. Containers
Use **[Docker](https://docs.docker.com/engine/install/)** on a typical Linux workstation (example below), or **[Singularity](https://sylabs.io/docs/) / [Apptainer](https://apptainer.org/docs/)** on many HPC clusters - then run with the matching Nextflow profile (e.g. `-profile docker`, `-profile singularity`).

**Docker** (example):
```bash
curl -fsSL https://get.docker.com/ | sudo sh
sudo usermod -aG docker "$USER"   # then log out and back in (or `newgrp docker`)
docker run --rm hello-world
```

## 4. Download databases
Download and unpack the required database (Viral RefSeq and human host genomes):

```
mkdir -p Databases && cd Databases
wget -c https://zenodo.org/records/13132915/files/combined_databases.tar.gz
tar -xzvf combined_databases.tar.gz
rm combined_databases.tar.gz
```

## 5. Configure paths (samplesheet, output, databases)

At the **repository root**, choose the params file to match how you start: **[`params_fastq.yaml`](params_fastq.yaml)** when you already have basecalled **reads (FASTQ)**, or **[`params_POD5.yaml`](params_POD5.yaml)** when you start from raw **POD5** signal (“squiggle”) data and need basecalling inside the pipeline. Edit that file using **absolute paths**.

| Setting | Purpose |
|-----|-------------------|
| `input` | Samplesheet CSV: copy **[`fastq.csv`](nf-metatropics/assets/submission/fastq.csv)** (FASTQ) or **[`POD5.csv`](nf-metatropics/assets/submission/POD5.csv)** (POD5) from [`nf-metatropics/assets/submission/`](nf-metatropics/assets/submission/), edit it, then set `input` to that file’s absolute path. 
| `outdir` | Where results are written. |
| `dbmeta` | Set database path (from step 4). |
| `Human_host_fasta` | Set human host FASTA path (from step 4). |
| `input_dir` | POD5 only: directory containing POD5 files (with **`params_POD5.yaml`**) |

Additional options: **[`nf-metatropics/assets/submission/all_options.md`](nf-metatropics/assets/submission/all_options.md)**.

## 6. Running Metatropics

With your `params_fastq.yaml` or `params_POD5.yaml` in place, run from the **repository root** (swap `-profile docker` for e.g. `-profile singularity`, if needed):

```
nextflow run nf-metatropics/ -profile docker -params-file params_fastq.yaml -resume
```

## 7. Output

Results are written under your chosen `--outdir` and summarized below:

### POD5 basecalling and demultiplexing (optional)

| Folder | Contents |
|--------|----------|
| `basecalling` | Intermediate FASTQ from Dorado before demultiplexing. |
| `demultiplexing` | Per-barcode FASTQ after demultiplexing. |

### Read preprocessing

This step removes low-quality reads and reads that match host background (e.g., human). The host-depleted read set is used for the rest of the pipeline.

| Folder | Contents |
|--------|----------|
| `fix` | Per-sample FASTQ after naming / format fixes (compressed). |
| `rarefaction` | Per-sample FASTQ after optional rarefaction subsampling. |
| `fastp` | Trimmed reads and FASTP reports. |
| `nanoplot` | Read-length and quality summaries (NanoPlot). |
| `multiqc` | **MultiQC** HTML report (`multiqc_report.html`)  |
| `nohuman` | FASTQ of reads **not** mapping to the human reference (human-depleted). |
| `nohost` | Optional: FASTQ after depletion against an extra host genome. |

### Taxonomic classification

Host-depleted reads are then mapped to the metagenomic database and summarized so you can see **which organisms (taxa) are present** in each sample.

| Folder | Contents |
|--------|----------|
| `metamaps` | MetaMaps mapping and classification outputs (`mapDirectly`, `Classify`). |
| `krona` | Krona HTML charts of taxonomic composition. |

### Per-virus reads

Reads that MetaMaps assigns to a given **virus** are extracted into separate FASTQs **per sample and per virus** for downstream work.

| Folder | Contents |
|--------|----------|
| `seqtk` | FASTQ per sample per virus (reads assigned to that virus). |

### Variant calling (BAMs and references)

For each candidate virus, reads are aligned to the **viral reference** (**BAMs**). **Variant calls** - differences from that reference - are derived from those alignments. 

| Folder | Contents |
|--------|----------|
| `reffix` | Reference FASTA with **cleaned headers** for each virus. |
| `medaka` | **BAMs** of reads aligned to each viral reference and **VCF** files with **variant calls** (differences from that reference). |
| `addingDepth` | Per-virus depth tables (coverage + consensus + classification). |

### Consensus

A **consensus genome** is called per virus (iVar draft, **Homopolish** final genome).

| Folder | Contents |
|--------|----------|
| `ivar` | Consensus sequences from iVar (input to Homopolish). |
| `homopolish` | **Polished** consensus FASTA (typical final genome per virus per sample). |

### Run-level summaries and provenance

Outputs that **summarize the whole run** (combined tables, coverage plots, read-count inputs, and pipeline provenance).

| Folder | Contents |
|--------|----------|
| `final` | Combined **final TSV** across the whole run. |
| `rcoverage` | Coverage **PDFs** for identified viruses (if enabled). |
| `read_count` | Read-count **CSV/PDF** and staged inputs for the read-distribution figure (aggregates material from `fix`, `fastp`, `nohuman`, `nohost`, `metamaps`—folder names must stay as the pipeline expects). |
| `pipeline_info` | Nextflow reports, command lines, and software versions. |

**Note:** For a typical outbreak-style use case, priority outputs are often **`homopolish`** (polished consensus), **`final`**, **`rcoverage`**, and **`read_count`**. 

## 8. High performance computing

You can run Metatropics on a **high-performance cluster** instead of a local workstation. For the Flemish Tier‑1 system **CalcUA** (VSC Antwerp), ready-made **Slurm** submission scripts and **Nextflow** profiles live under [`nf-metatropics/assets/calcua/`](nf-metatropics/assets/calcua/). **Setup, editing the batch scripts, and `sbatch` commands** are documented in the [CalcUA README](nf-metatropics/assets/calcua/README.md).

## 9. Citation

If you use Metatropics in your research, please cite:

```
De Souza Novaes, A., Jansen, D., de Block, T., Vercauteren, K., & Rezende, A. M. (2026). Metatropics: Human viral pathogen identification and consensus genome calling from nanopore metagenomic sequencing data (Version 0.0.5). GitHub. https://github.com/DaanJansen94/Metatropics
```

Also cite the **nf-core** framework, and other tools you rely on; see [`nf-metatropics/assets/citing/CITATIONS.md`](nf-metatropics/assets/citing/CITATIONS.md).