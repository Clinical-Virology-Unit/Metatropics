[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metatropics-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/metatropics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# Metatropics: A Viral Metagenomics ONT Pipeline

**Metatropics (v0.1.1)** is a Nextflow-driven bioinformatics pipeline for metagenomic Oxford Nanopore sequencing data. It is built to detect viral pathogens in complex samples (e.g. blood and swabs from a range of body sites) and, where coverage allows, to generate high-quality viral consensus genomes and to perform variant analysis.

**Metatropics** is an abbreviation of **Metagenomics for Tropical Fevers**, and reflects how the project began, with an emphasis on finding human viral pathogens in patients presenting with tropical fevers. The same pipeline has since been validated and applied outside that first setting, including for other febrile syndromes, for genomic surveillance, and for research and diagnostic questions around viral pathogens relevant to human health.

## Why Metatropics?

- **High sensitivity across viral diversity**: maps against large reference databases (beyond small sets like RefSeq), improving read mapping/quantification and sensitivity for highly diverse human viruses (e.g., arenaviruses such as Lassa virus) and validated in field use.
- **High specificity / low false positives**: uses [Virasign](https://github.com/DaanJansen94/virasign) viral classification (two-step design that filters background-like calls and improves specificity; see Virasign for details) plus additional evidence metrics (consensus breadth, [NOGR](./nf-metatropics/assets/output/NOGR.md), [Z-score](./nf-metatropics/assets/output/Z_SCORE.md)).
- **High-accuracy variant calling and consensus**: uses Clair3 for ONT variant calling; benchmarks show ONT+Clair3 can match or exceed Illumina short-read calling accuracy (see eLife [`98300`](https://elifesciences.org/articles/98300)), producing the highest quality viral consensus genomes.
- **Fast interpretation and visual confirmation**: an interactive, easily shareable HTML report makes it straightforward to confirm true viral hits versus false positives (e.g. amplicon contamination).

For details on outputs and how to read the report, see [`nf-metatropics/assets/output/README.md`](./nf-metatropics/assets/output/README.md).

---

## Pipeline summary

![Figure](./nf-metatropics/assets/logo/MetaTropics_workflow.png)

---

## 1. Clone the repository
```bash
git clone https://github.com/Clinical-Virology-Unit/Metatropics.git
cd Metatropics
```

---

## 2. Java and Nextflow
You need **Java 17+** and **[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) 25.04.6**. On Debian/Ubuntu you can do:
```bash
sudo apt update && sudo apt install -y openjdk-17-jdk curl
curl -sSL https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/
nextflow self-update 25.04.6
nextflow -version
# If `self-update` fails: add `export NXF_VER=25.04.6` to `~/.bashrc`
```

---

## 3. Containers
Use **[Docker](https://docs.docker.com/engine/install/)** on a typical Linux workstation (example below), or **[Singularity](https://sylabs.io/docs/) / [Apptainer](https://apptainer.org/docs/)** on many HPC clusters - then run with the matching Nextflow profile (e.g. `-profile docker`, `-profile singularity`).

**Docker** (example):
```bash
curl -fsSL https://get.docker.com/ | sudo sh
sudo usermod -aG docker "$USER"   # then log out and back in (or `newgrp docker`)
docker run --rm hello-world
```

---

## 4. Configure paths (input, output)

At the repository root, choose the params file to match how you start and edit it using absolute paths.

### POD5 start (basecalling inside the pipeline)

Use **[`params_POD5.yaml`](params_POD5.yaml)**.

| Setting | Purpose |
|-----|-------------------|
| `input` | Copy **[`POD5.csv`](nf-metatropics/assets/submission/POD5.csv)** from [`nf-metatropics/assets/submission/`](nf-metatropics/assets/submission/), edit it, then set `input` to that file’s absolute path. |
| `input_dir` | Directory containing POD5 files. |
| `outdir` | Where results are written. |
| `kit_name` | Dorado `--kit-name` (default: `TWIST-96A-UDI`). |
| `Host` | Optional: host depletion (e.g., `human,pan`). |
| `virasign_ultrasensitive` | Optional: enable ultrasensitive viral identification mode. |

To autogenerate a [`POD5.csv`](nf-metatropics/assets/submission/POD5.csv) template, run **`pip install .`** once at the repo root, then **`metatropics-samplesheet pod5 -i .`** from your POD5 directory. For TWIST UDI plates, run **`metatropics-samplesheet pod5 TWIST-96A-UDI`** to create `run.txt`, edit sample names and wells, then **`metatropics-samplesheet pod5 TWIST-96A-UDI run.txt`** to build `POD5.csv`.

### fastq_pass start (on-device basecalled, demultiplex inside the pipeline)

Use **[`params_fastq_pass.yaml`](params_fastq_pass.yaml)**.

| Setting | Purpose |
|-----|-------------------|
| `input` | Copy **[`POD5.csv`](nf-metatropics/assets/submission/POD5.csv)** (same barcode mapping as POD5 mode), edit it, then set `input` to that file’s absolute path. |
| `input_dir` | Directory with combined basecalled reads (e.g. `fastq_pass`). |
| `outdir` | Where results are written. |
| `basecall` | Must be `false` so the pipeline skips Dorado basecalling and only demultiplexes. |
| `kit_name` | Dorado `--kit-name` (default: `TWIST-96A-UDI`). |
| `Host` | Optional: host depletion (e.g., `human,pan`). |
| `virasign_ultrasensitive` | Optional: enable ultrasensitive viral identification mode. |

Use the same **`metatropics-samplesheet pod5 -i .`** command from your `fastq_pass` folder to create a `POD5.csv` template there. For TWIST UDI plates, run **`metatropics-samplesheet pod5 TWIST-96A-UDI`** to create `run.txt`, edit sample names and wells, then **`metatropics-samplesheet pod5 TWIST-96A-UDI run.txt`** to build `POD5.csv`.

### FASTQ start (basecalled reads)

Use **[`params_fastq.yaml`](params_fastq.yaml)**.

| Setting | Purpose |
|-----|-------------------|
| `input` | Copy **[`fastq.csv`](nf-metatropics/assets/submission/fastq.csv)** from [`nf-metatropics/assets/submission/`](nf-metatropics/assets/submission/), edit it, then set `input` to that file’s absolute path. |
| `outdir` | Where results are written. |
| `Host` | Optional: host depletion (e.g., `human,pan`). |
| `virasign_ultrasensitive` | Optional: enable ultrasensitive viral identification mode. |

To autogenerate a [`fastq.csv`](nf-metatropics/assets/submission/fastq.csv) from a reads folder, run **`pip install .`** once at the repo root, then **`metatropics-samplesheet -i .`** from your FASTQ directory.

Additional options: **[`nf-metatropics/assets/submission/all_options.md`](nf-metatropics/assets/submission/all_options.md)**.

---

## 5. Running Metatropics

With your params file in place, run from the repository root (swap `-profile docker` for e.g. `-profile singularity`, if needed):

```bash
nextflow run nf-metatropics/ -profile docker -params-file params_fastq.yaml -resume
```

---

## 6. Output

Results are written under your chosen `--outdir` and summarized below:

| Group | Role |
|--------|------|
| **`Basecalling/`** | Dorado basecalling and demultiplexing (optional). |
| **`Reads/`** | Read QC, trimming, human / optional host depletion. |
| **`Classification/`** | Virasign viral classification outputs and reports. |
| **`Variant_calling/`** | Clair3 variant calls (VCFs) and HTML report. |
| **`Consensus/`** | bcftools consensus genome. |
| **`Summary/`** | Final Metatropics report (`Summary/metatropics/Metatropics_Summary_RVDB.html`) listing all identified viruses, plus read-count summaries and pipeline provenance. |

For a detailed description of each output subfolder, see [`nf-metatropics/assets/output/README.md`](nf-metatropics/assets/output/README.md).

---

## 7. High performance computing

You can run Metatropics on a high-performance cluster instead of a local workstation. For the Flemish Tier‑1 system CalcUA (VSC Antwerp), ready-made Slurm submission scripts and Nextflow profiles live under [`nf-metatropics/assets/calcua/`](nf-metatropics/assets/calcua/). Setup, editing the batch scripts, and `sbatch` commands are documented in the [CalcUA README](nf-metatropics/assets/calcua/README.md).

---

## 8. Citation

If you use Metatropics in your research, please cite:

```
De Souza Novaes, A.†, Jansen, D.†, de Block, T., Coppens, S., Ariën, K.K., Selhorst, P., Rezende, A.M.*, & Vercauteren, K*. (2026). Metatropics: Human viral pathogen identification and consensus genome calling from nanopore metagenomic sequencing data. (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.20430617
```

† These authors contributed equally (shared first authorship).  
\* These authors contributed equally (shared senior authorship).

Also cite the nf-core framework, and other tools you rely on; see [`nf-metatropics/assets/citing/CITATIONS.md`](nf-metatropics/assets/citing/CITATIONS.md).
