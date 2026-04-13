# Pipeline parameters

## Input/output options

| Option | Type | Description |
|--------|------|-------------|
| `--input` | string | Path to the samplesheet CSV. *(FASTQ: run `pip install .` at the repository root, then **`metatropics-samplesheet -i .`** in the folder that contains your reads)* |
| `--input_dir` | string | Input directory with POD5. Default: None. |
| `--outdir` | string | Output directory for results. Use absolute paths on cloud storage. |

## Reference genome options

| Option | Type | Description |
|--------|------|-------------|
| `--Human_host_fasta` | string | Optional FASTA for the human background removal step. |
| `--Other_host_fasta` | string | Optional FASTA for an additional host background (e.g. mosquito, primate). |
| `--dbmeta` | string | Path for the MetaMaps database for read classification. Default: None. |

## Generic options

| Option | Type | Description |
|--------|------|-------------|
| `--basecall` | boolean | If POD5 is the input, set true. Default: false. |
| `--model` | string | POD5: Dorado model (`fast`, `hac`, or `sup`). Default: hac. |
| `--kit_name` | string | POD5: Dorado `--kit-name`. Allowed: `EXP-NBD103`, `EXP-NBD104`, `EXP-NBD114`, `EXP-NBD114-24`, `EXP-NBD196`, `EXP-PBC001`, `EXP-PBC096`, `SQK-16S024`, `SQK-16S114-24`, `SQK-LWB001`, `SQK-MLK111-96-XL`, `SQK-MLK114-96-XL`, `SQK-NBD111-24`, `SQK-NBD111-96`, `SQK-NBD114-24`, `SQK-NBD114-96`, `SQK-PBK004`, `SQK-PCB109`, `SQK-PCB110`, `SQK-PCB111-24`, `SQK-PCB114-24`, `SQK-RAB201`, `SQK-RAB204`, `SQK-RBK001`, `SQK-RBK004`, `SQK-RBK110-96`, `SQK-RBK111-24`, `SQK-RBK111-96`, `SQK-RBK114-24`, `SQK-RBK114-96`, `SQK-RLB001`, `SQK-RPB004`, `SQK-RPB114-24`, `TWIST-16-UDI`, `TWIST-96A-UDI`, `VSK-PTC001`, `VSK-VMK001`, `VSK-VMK004`, `VSK-VPS001`. Default: `TWIST-96A-UDI`. |
| `--minLength` | integer | Minimum read length to analyse. Default: 200. |
| `--minVirus` | number | Minimum virus frequency in raw data to report. Default: 0.01. |
| `--usegpu` | boolean | POD5: use NVIDIA GPU for basecalling. |
| `--pair` | boolean | Barcodes on both read ends (true) or one end (false). |
| `--quality` | integer | Minimum base quality for consensus. Default: 7. |
| `--agreement` | number | Minimum base frequency for unambiguous consensus calls. Default: 0.7. |
| `--depth` | integer | Minimum per-position depth for consensus. Default: 5. |
| `--front` | integer | Bases to trim at 5′. Default: 0. |
| `--tail` | integer | Bases to trim at 3′. Default: 0. |
| `--rcoverage` | string | Coverage figures. Default: false. |
| `--horizontal_coverage` | integer | Minimum horizontal coverage threshold. Default: 1. |

## Rarefaction options

| Option | Type | Description |
|--------|------|-------------|
| `--perform_rarefaction` | boolean | Rarefy each sample to a target number of bases. Default: false. |
| `--target_bases` | number | Target bases per sample when rarefying (e.g. ~500k reads × 2 kb). Default: 1e9 bases. |

## Docker cleanup

| Option | Type | Description |
|--------|------|-------------|
| `--enable_docker_cleanup` | boolean | Remove downloaded Docker images after the run to free disk space. Default: false. |
