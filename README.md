[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metatropics-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/metatropics)

# Metatropics
The metatropics pipeline is a [Nextflow](https://www.nextflow.io/)-driven workflow designed for viral identification and the creation of consensus genomes from nanopore (metagenomic) sequencing data. It leverages container systems like [Docker](https://www.docker.com) and [Singularity](https://sylabs.io/docs/), utilizing one container per process to avoid software dependency conflicts and simplifies maintainenance. This container-based approach ensures that installation is straightforward and results are highly reproducible.

### Pipeline summary

![Figure](./nf-metatropics/assets/logo/Metatropics.jpg)

### 1. Clone the repository
```bash
git clone https://github.com/DaanJansen94/Metatropics.git
cd Metatropics
```

### 2. Java and Nextflow
You need **Java 17+** and **[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) ≥ 22.10.1**. On Debian/Ubuntu you can do:
```bash
sudo apt update && sudo apt install -y openjdk-17-jdk curl
curl -sSL https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/
nextflow -version
```

### 3. Containers
Use **[Docker](https://docs.docker.com/engine/install/)** on a typical Linux workstation (example below), or **[Singularity](https://sylabs.io/docs/) / [Apptainer](https://apptainer.org/docs/)** on many HPC clusters - then run with the matching Nextflow profile (e.g. `-profile docker`, `-profile singularity`).

**Docker** (example):
```bash
curl -fsSL https://get.docker.com/ | sudo sh
sudo usermod -aG docker "$USER"   # then log out and back in (or `newgrp docker`)
docker run --rm hello-world
```

### 4. Download databases
Download and unpack the required database ((Viral RefSeq and human host genomes):

```
mkdir -p Databases && cd Databases
wget -c https://zenodo.org/records/13132915/files/combined_databases.tar.gz
tar -xzvf combined_databases.tar.gz
rm combined_databases.tar.gz
```

### 5. Set PATHs

**Note:** Edit **[`params_fastq.yaml`](params_fastq.yaml)** for FASTQ input and **[`params_POD5.yaml`](params_POD5.yaml)** for POD5 / raw squiggle data (both at the **Metatropics** repo root). Example **samplesheets only** are in **[`nf-metatropics/assets/submission/`](nf-metatropics/assets/submission/)** — see [`README.md` there](nf-metatropics/assets/submission/README.md).

```
# from the repository root (after step 1)
nano params_fastq.yaml   # or params_POD5.yaml

input: absolute path to your samplesheet — use or copy [`nf-metatropics/assets/submission/fastq.csv`](nf-metatropics/assets/submission/fastq.csv) or [`POD5.csv`](nf-metatropics/assets/submission/POD5.csv)
outdir: absolute path to your output directory
Human_host_fasta: optional path to the human background genome
Other_host_fasta: optional path to any additional host genome
dbmeta: path to your ViralRefseq (MetaMaps) database
quality: 30 # for high-quality genomes
depth: 20 # for high-quality genomes
```

### 6. Set Input

**Note:** Ensure that each sample’s FASTQ is consolidated (e.g. one file per barcode) before you list paths in the samplesheet. **Example CSVs** are under **[`nf-metatropics/assets/submission/`](nf-metatropics/assets/submission/)** — see that [README](nf-metatropics/assets/submission/README.md) (`fastq.csv` vs `POD5.csv`). **Params** stay at the repo root: **[`params_fastq.yaml`](params_fastq.yaml)** / **[`params_POD5.yaml`](params_POD5.yaml)**. The CSV columns differ by starting data:
- For raw reads (FASTQ): **[`params_fastq.yaml`](params_fastq.yaml)** + samplesheet **[`fastq.csv`](nf-metatropics/assets/submission/fastq.csv)** (copy and edit paths)
```
The params file contains the most important paths
sample,single_end,barcode
sample_name01,True,/home/antonio/metatropics/nf-metatropics/fastq/barcode01.fastq
sample_name02,True,/home/antonio/metatropics/nf-metatropics/fastq/barcode02.fastq
```

- For squiggle data (POD5): **[`params_POD5.yaml`](params_POD5.yaml)** + samplesheet **[`POD5.csv`](nf-metatropics/assets/submission/POD5.csv)** (copy and edit; set `input_dir` to your POD5 directory)
```
sample,single_end,barcode
sample_name01,True,barcode01
sample_name02,True,barcode02
```

### 7. Runing pipeline

```
nextflow run nf-metatropics/ -profile docker -params-file params_fastq.yaml -resume
```

   ```
   nextflow run nf-metatropics/ --help

   Input/output options
    --input                       [string]  Path to comma-separated file containing information about the samples in the experiment.
    --input_dir                   [string]  Input directory with POD5 [default: None]
    --outdir                      [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
Reference genome options
    --Human_host_fasta            [string]  Optional FASTA for the human background removal step.
    --Other_host_fasta            [string]  Optional FASTA for an additional host background (e.g. mosquito, primate).
    --dbmeta                      [string]  Path for the MetaMaps database for read classification [default: None]
   Generic options
    --basecall                    [boolean] In case POD5 is the input, that option shoud be true [default: false]
    --model                       [string]  In case POD5 is the input, the dorado model for basecalling should be provided. Choose from: fast, hac, or sup [default: hac]
    --kit_name                    [string]  In case POD5 is the input, the kit name should be provided. Available options include: EXP-NBD103, EXP-NBD104, EXP-NBD114, EXP-NBD114-24, EXP-NBD196, EXP-PBC001, EXP-PBC096, SQK-16S024, SQK-16S114-24, SQK-LWB001, SQK-MLK111-96-XL, SQK-MLK114-96-XL, SQK-NBD111-24, SQK-NBD111-96, SQK-NBD114-24, SQK-NBD114-96, SQK-PBK004, SQK-PCB109, SQK-PCB110, SQK-PCB111-24, SQK-PCB114-24, SQK-RAB201, SQK-RAB204, SQK-RBK001, SQK-RBK004, SQK-RBK110-96, SQK-RBK111-24, SQK-RBK111-96, SQK-RBK114-24, SQK-RBK114-96, SQK-RLB001, SQK-RPB004, SQK-RPB114-24, TWIST-16-UDI, TWIST-96A-UDI, VSK-PTC001, VSK-VMK001, VSK-VMK004, VSK-VPS001 [default: TWIST-96A-UDI]
    --minLength                   [integer] Minimum length for a read to be analyzed. [default: 200]
    --minVirus                    [number]  Minimum virus data frequency in the raw data to be part of the output. [default: 0.01]
    --usegpu                      [boolean] In case POD5 is the input, the use of GPU Nvidia should be true.
    --pair                        [boolean] If barcodes were added at both sides of a read (true) or only at one side (false).
    --quality                     [integer] Minimum quality for a base to build the consensus [default: 7]
    --agreement                   [number]  Minimum base frequency to be called without ambiguit to build the consensus [default: 0.7]
    --depth                       [integer] Minimum depth of a position to build the consensus [default: 5]
    --front                       [integer] Number of bases to delete at 5 prime of the read [default: 0]
    --tail                        [integer] Number of bases to delete at 3 prime of the read [default: 0]
    --rcoverage                   [string]  Coverage figures [default: false]
    --horizontal_coverage         [integer] Minimum horizontal coverage threshold [default: 1]
   Rarefaction options
    --perform_rarefaction         [boolean] Option to perform rarefaction to a specified number of bases [default: false]
    --target_bases                [number]  Number of bases to which you want to rarefy each sample [default: 1 billion bases, equivalent to 500,000 reads of 2kb each]
   Docker cleanup 
    --enable_docker_cleanup       [boolean] Removes all downloaded Docker images to free up root space [default: false]
   ```
**Note:** If internet access is unavailable, disable the docker cleanup option to retain images after the initial download, allowing the pipeline to run without internet access.

### 8. Output
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

### 9. High performance computing

You can also run this pipeline on an HPC cluster; for example, on the Flemish Tier-1 system **CalcUA** (VSC), use the Slurm submission scripts and Nextflow profile under [`nf-metatropics/assets/calcua/`](nf-metatropics/assets/calcua/). For setup and usage, see the [README in that folder](nf-metatropics/assets/calcua/README.md).

## Citation

If you use Metatropics in your research, please cite:

```
De Souza Novaes, A., Jansen, D., de Block, T., Vercauteren, K., & Rezende, A. M. (2026). Metatropics: Human viral pathogen identification and consensus genome calling from nanopore metagenomic sequencing data (Version 0.0.5). GitHub. https://github.com/DaanJansen94/Metatropics
```

Also cite the **nf-core** framework, and other tools you rely on; see [`nf-metatropics/assets/Citing/CITATIONS.md`](nf-metatropics/assets/Citing/CITATIONS.md).

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome. Please open a pull request or an issue to discuss larger changes first.

## Support

If you run into problems or have questions, open an issue on [GitHub](https://github.com/DaanJansen94/Metatropics/issues).
