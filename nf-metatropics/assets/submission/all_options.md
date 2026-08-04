# Pipeline parameters

## Input/output options

### FASTQ start (basecalled reads)

| Option | Description |
|--------|-------------|
| `--input` | Path to the samplesheet CSV (example: [`fastq.csv`](./fastq.csv)) |
| `--outdir` | Output directory for results. |

To autogenerate a `samplesheet.csv` from a folder of demultiplexed FASTQ files, run `pip install .` once at the repo root, then `metatropics-samplesheet -i .` from your FASTQ directory.

### POD5 start (basecalling inside the pipeline)

| Option | Description |
|--------|-------------|
| `--input` | Path to the samplesheet CSV (example: [`POD5.csv`](./POD5.csv)). |
| `--input_dir` | Input directory with POD5. Default: None. |
| `--outdir` | Output directory for results. Use absolute paths on cloud storage. |
| `--model` | Dorado model (`fast`, `hac`, or `sup`). Default: `hac`. |
| `--kit_name` | Dorado `--kit-name`. Allowed: `EXP-NBD103`,`EXP-NBD104`,`EXP-NBD114`,`EXP-NBD114-24`,`EXP-NBD196`,`EXP-PBC001`,`EXP-PBC096`,`SQK-16S024`,`SQK-16S114-24`,`SQK-LWB001`,`SQK-MLK111-96-XL`,`SQK-MLK114-96-XL`,`SQK-NBD111-24`,`SQK-NBD111-96`,`SQK-NBD114-24`,`SQK-NBD114-96`,`SQK-PBK004`,`SQK-PCB109`,`SQK-PCB110`,`SQK-PCB111-24`,`SQK-PCB114-24`,`SQK-RAB201`,`SQK-RAB204`,`SQK-RBK001`,`SQK-RBK004`,`SQK-RBK110-96`,`SQK-RBK111-24`,`SQK-RBK111-96`,`SQK-RBK114-24`,`SQK-RBK114-96`,`SQK-RLB001`,`SQK-RPB004`,`SQK-RPB114-24`,`TWIST-16-UDI`,`TWIST-96A-UDI`,`TWIST-96B-UDI`,`TWIST-96C-UDI`,`TWIST-96D-UDI`,`VSK-PTC001`,`VSK-VMK001`,`VSK-VMK004`,`VSK-VPS001`. Default: `TWIST-96A-UDI`. |

To autogenerate a `POD5.csv` template, run `metatropics-samplesheet pod5 -i .` from your POD5 directory. For TWIST UDI plates, run `metatropics-samplesheet pod5 TWIST-96A-UDI` to create `run.txt`, edit sample names and wells, then `metatropics-samplesheet pod5 TWIST-96A-UDI run.txt` to build `POD5.csv`.

### fastq_pass start (on-device basecalled, demultiplex inside the pipeline)

| Option | Description |
|--------|-------------|
| `--input` | Path to the samplesheet CSV (example: [`POD5.csv`](./POD5.csv)). |
| `--input_dir` | Path to the `fastq_pass` folder with basecalled reads (not yet demultiplexed). |
| `--basecall` | Must be `false` to skip pipeline basecalling and run demultiplexing only. |
| `--outdir` | Output directory for results. |
| `--kit_name` | Dorado `--kit-name` (e.g. `TWIST-96C-UDI`). Pipeline container: Dorado 2.0. Default: `TWIST-96A-UDI`. |

Use `metatropics-samplesheet pod5 -i .` from your `fastq_pass` folder to create a `POD5.csv` template there. For TWIST UDI plates, run `metatropics-samplesheet pod5 TWIST-96A-UDI` to create `run.txt`, edit sample names and wells, then `metatropics-samplesheet pod5 TWIST-96A-UDI run.txt` to build `POD5.csv`.

Basecalling on the instrument should keep barcodes intact (Dorado `--no-trim`); otherwise demultiplexing may classify most reads as unclassified.

## Read processing options

| Option | Description |
|--------|-------------|
| `--minLength` | Min read length to analyse. Default: 200. |
| `--quality` | Min base quality used by QC and variant filtering. Default: 15. |
| `--front` | Bases to trim at 5′. Default: 25. |
| `--tail` | Bases to trim at 3′. Default: 25. |

## Rarefaction options

| Option | Description |
|--------|-------------|
| `--rarefaction` | Rarefy each sample to a target number of bases. Default: true. |
| `--target_bases` | Target bases per sample when rarefying (e.g. ~1M reads × 500 bp). Default: 5e8 bases. |

## Host depletion option

| Option | Description |
|--------|-------------|
| `--Human_host_fasta` | Optional FASTA for the human background removal step. |
| `--Other_host_fasta` | Optional FASTA for an additional host background (e.g. mosquito, primate). |
| `--Host` | Optional host keyword(s) to auto-download FASTA under `Metatropics/Databases`. Multiple hosts are supported and will be merged (e.g. `human,aedes,culex`). Supported keywords include `human`, `pan`, `gorilla`, `orangutan`, `macaque`, `aedes`, `anopheles`, `culex`, `bat`, `rat`, `dog`, `cat`, `camel`, `goat`, `pig`, `cow`, `mouse`, `chicken`. |

## Viral classifier options (Virasign)

### Important options

| Option | Description |
|--------|-------------|
| `--virasign_database` | Database(s) (e.g. `RVDB`, `RefSeq`, or `RVDB,RefSeq`). Default: `RVDB`. |
| `--virasign_db_dir` | Where Virasign databases are stored (defaults to `Databases/`). |
| `--virasign_ultrasensitive` | Enable Virasign ultrasensitive mode. Default: false. |
| `--virasign_blind` | Blind specific viral species from analysis (not in any output). Use Virasign abbreviations (e.g. `HEP,HIV,HTLV,EBV,CMV,HPV`) or full species names (comma-separated). To list abbreviations: `virasign --blinding`. |

### Advanced / less common options

| Option | Description |
|--------|-------------|
| `--virasign_rvdb_version` | RVDB release to use (optional). |
| `--virasign_accessions` | Extra accessions to include (optional). |
| `--virasign_min_identity` | Min alignment identity threshold (optional). |
| `--virasign_min_mapped_reads` | Min mapped reads to report a hit (optional). |
| `--virasign_coverage_depth` | Min per-position depth for coverage filtering (optional). |
| `--virasign_coverage_breadth` | Min breadth (fraction) for coverage filtering (optional). |
| `--virasign_min_nogr` | Min number of non-overlapping genomic regions (NoGR) required (optional). |
| `--virasign_zscore` | Enable/disable Z-score background correction (optional; **default: enabled / `true`**). Set to `false` to disable. |
| `--virasign_zscore_controls` | Override auto-detected Z-score controls with **sample names and/or FASTQ paths** (≥2). Comma-separated or a text file (one entry per line). |
| `--virasign_threads` | Threads for Virasign (optional; defaults to the task CPUs). |
| `--virasign_ram_gb` | Minimap2 RAM/GB hint (`-I`) for Virasign (optional). |
| `--virasign_enable_clustering` | Enable clustering of references in database prep. Default: false. |
| `--virasign_cluster_identity` | Clustering identity threshold (optional). |
| `--virasign_max_ambiguous_fraction` | Max allowed ambiguous fraction when preparing DB (optional). |

## Variant calling options and consensus thresholds (Clair3)

| Option | Description |
|--------|-------------|
| `--clair3_model` | Optional override for the Nanopore model. |
| `--clair3_min_mq` | Min mapping quality (MAPQ). Default: 15. |
| `--clair3_min_bq` | Min base quality (BQ). Default: 15. |
| `--clair3_min_alt_reads` | Min ALT-supporting reads. Default: 10. |
| `--depth` | Min per-position depth. Default: 25. |
| `--agreement` | Min VAF for applying variants into the consensus. Default: 0.7. |
