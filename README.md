[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metatropics-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/metatropics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# Metatropics

**Metatropics** is a nextflow-driven bioinformatics pipeline for metagenomic Oxford Nanopore sequencing data. It is built to detect viral pathogens in complex samples and, where coverage allows, to generate high-quality viral consensus genomes.

**Metatropics** is an abbreviation of **Metagenomics for Tropical Fevers**, and reflects how to project began, with an emphasis on finding human viral pathogens in patients presenting with tropical fevers. The same pipeline has since been validated and applied outside that first setting, including for other febrile syndromes, for genomic surveillance, and for research and diagnostic questions around viral pathogens relevant to human health.

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
Below one can see the output directories and their description. `basecalling` and `demultiplexing` will exist only in case the user has used POD5 files as input.

1. [`basecalling`] - fastq files after the basecalling without being demultiplexed
2. [`demultiplexing`] - directories and fastq files produced after demultiplexing
3. [`fix`] - gziped fastq files for each sample of the run
4. [`rarefaction`] - gziped fastq files for each rarefied sample
5. [`fastp`] - results after trimming analysis performed by FASTP
6. [`nanoplot`] - quality results for the sequencing data just after demultiplexing
7. [`minimap2`] - BAM files about mapping against host genome
8. [`nohuman`] - gziped fastq files without reads mapping to human genome
9. [`nohost`] - gziped fastq files without reads mapping to host genome (-optional)
10. [`metamaps`] - results from both steps of Metamaps execution for read classification (mapDirectly and Classify)
11. [`r`] - intermediate table report and graphical PDF report for each sample
12. [`ref`] - header of the reads and fasta reference genomes for each virus found for each sample
13. [`krona`] - HTML files for each sample with interactive composition pie chart
14. [`reffix`] - fasta refence genomes with fixed header for each virus found during the run
15. [`seqtk`] - gziped fastq file for each set of read classified to a virus for each sample
16. [`medaka`] - BAM file for each virus with mapping results from the virus genome reference for each sample
17. [`samtools`] - mapping statistics calculated to BAM files present in the `medaka` directory
18. [`ivar`] - consensus sequences produced for each virus found in each sample
19. [`bam`] - detailed statistics for the BAM files from `medaka` directory for each position of virus refence genome
20. [`homopolish`] - consensus sequence for each virus in each sample polished for the indel variations
21. [`addingDepth`] - table report for each virus in each sample
22. [`final`] - final table report for all the run
23. [`pipeline_info`] - reports on the execution of the pipeline produced by NextFlow
24. [`rcoverage`] - PDF files including coverage figures of identified viruses
25. [`read_count`] - PDF and CSV files representing read distribution. These figures visualize the distribution of all reads, including trimmed, human, viral, and other reads.

**Note:** For the INRB mpox analysis, the most important files are the polished consensus sequences (20), the final report (22), the coverage (24) and read distribution figures (25). 

Tip 1: If you have limited space, you can delete the 'work' directory and, after selecting the necessary output files, also remove the 'output' directory.

Tip 2: When you encounter errors, make sure to double-check the memory allocated to your processes. This is often the cause, or alternatively, consider including rarefaction.

## 8. High performance computing

You can also run this pipeline on an HPC cluster; for example, on the Flemish Tier-1 system **CalcUA** (VSC), use the Slurm submission scripts and Nextflow profile under [`nf-metatropics/assets/calcua/`](nf-metatropics/assets/calcua/). For setup and usage, see the [README in that folder](nf-metatropics/assets/calcua/README.md).

## 9. Citation

If you use Metatropics in your research, please cite:

```
De Souza Novaes, A., Jansen, D., de Block, T., Vercauteren, K., & Rezende, A. M. (2026). Metatropics: Human viral pathogen identification and consensus genome calling from nanopore metagenomic sequencing data (Version 0.0.5). GitHub. https://github.com/DaanJansen94/Metatropics
```

Also cite the **nf-core** framework, and other tools you rely on; see [`nf-metatropics/assets/citing/CITATIONS.md`](nf-metatropics/assets/citing/CITATIONS.md).