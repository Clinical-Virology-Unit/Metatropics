# Pipeline parameters

## Input/output options

### FASTQ start (basecalled reads)

| Option | Description |
|--------|-------------|
| `--input` | Path to the samplesheet CSV (example: [`fastq.csv`](./fastq.csv)) |
| `--outdir` | Output directory for results. |

To autogenerate a `samplesheet.csv` from a reads folder, run `pip install .` once at the repo root, then `metatropics-samplesheet -i .` from your FASTQ directory.

### POD5 start (basecalling inside the pipeline)

| Option | Description |
|--------|-------------|
| `--input` | Path to the samplesheet CSV (example: [`POD5.csv`](./POD5.csv)). |
| `--input_dir` | Input directory with POD5. Default: None. |
| `--outdir` | Output directory for results. Use absolute paths on cloud storage. |
| `--model` | Dorado model (`fast`, `hac`, or `sup`). Default: `hac`. |
| `--kit_name` | Dorado `--kit-name`. Allowed: `EXP-NBD103`, `EXP-NBD104`, `EXP-NBD114`, `EXP-NBD114-24`, `EXP-NBD196`, `EXP-PBC001`, `EXP-PBC096`, `SQK-16S024`, `SQK-16S114-24`, `SQK-LWB001`, `SQK-MLK111-96-XL`, `SQK-MLK114-96-XL`, `SQK-NBD111-24`, `SQK-NBD111-96`, `SQK-NBD114-24`, `SQK-NBD114-96`, `SQK-PBK004`, `SQK-PCB109`, `SQK-PCB110`, `SQK-PCB111-24`, `SQK-PCB114-24`, `SQK-RAB201`, `SQK-RAB204`, `SQK-RBK001`, `SQK-RBK004`, `SQK-RBK110-96`, `SQK-RBK111-24`, `SQK-RBK111-96`, `SQK-RBK114-24`, `SQK-RBK114-96`, `SQK-RLB001`, `SQK-RPB004`, `SQK-RPB114-24`, `TWIST-16-UDI`, `TWIST-96A-UDI`, `VSK-PTC001`, `VSK-VMK001`, `VSK-VMK004`, `VSK-VPS001`. Default: `TWIST-96A-UDI`. |

## Read processing options

| Option | Description |
|--------|-------------|
| `--minLength` | Minimum read length to analyse. Default: 200. |
| `--quality` | Minimum base quality for consensus. Default: 7. |
| `--agreement` | Minimum base frequency for unambiguous consensus calls. Default: 0.7. |
| `--depth` | Minimum per-position depth for consensus. Default: 5. |
| `--front` | Bases to trim at 5′. Default: 0. |
| `--tail` | Bases to trim at 3′. Default: 0. |

## Rarefaction options

| Option | Description |
|--------|-------------|
| `--perform_rarefaction` | Rarefy each sample to a target number of bases. Default: false. |
| `--target_bases` | Target bases per sample when rarefying (e.g. ~500k reads × 2 kb). Default: 1e9 bases. |

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
| `--virasign_min_identity` | Minimum alignment identity threshold (optional). |
| `--virasign_min_mapped_reads` | Minimum mapped reads to report a hit (optional). |
| `--virasign_coverage_depth` | Minimum per-position depth for coverage filtering (optional). |
| `--virasign_coverage_breadth` | Minimum breadth (fraction) for coverage filtering (optional). |
| `--virasign_min_nogr` | Min number of non-overlapping genomic regions (NoGR) required (optional). |
| `--virasign_zscore` | Enable/disable z-score filtering (optional). |
| `--virasign_zscore_controls` | Path to negative-control FASTQs for z-score (optional). |
| `--virasign_threads` | Threads for Virasign (optional; defaults to the task CPUs). |
| `--virasign_ram_gb` | Minimap2 RAM/GB hint (`-I`) for Virasign (optional). |
| `--virasign_enable_clustering` | Enable clustering of references in database prep. Default: false. |
| `--virasign_cluster_identity` | Clustering identity threshold (optional). |
| `--virasign_max_ambiguous_fraction` | Max allowed ambiguous fraction when preparing DB (optional). |
