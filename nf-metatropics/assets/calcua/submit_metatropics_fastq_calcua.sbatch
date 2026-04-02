#!/bin/bash -l
#SBATCH --job-name=mt_fastq
#SBATCH --partition=zen3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=72:00:00
#SBATCH --account=ap_itg_cvu
#SBATCH --output=metatropics_fastq_%j.out
#SBATCH --error=metatropics_fastq_%j.err

set -euo pipefail

# Pipeline checkout and params: FASTQ workflow (params must have basecall: false, no input_dir).
PIPELINE_DIR="${VSC_SCRATCH}/Metatropics_runs/Metatropics/nf-metatropics"
PARAMS_FILE="${VSC_SCRATCH}/Metatropics_runs/Metatropics/params_fastq.yaml"

unset JAVA_HOME || true
unset JAVA_CMD || true

module load Java
module load Nextflow

export APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/cache"
export APPTAINER_TMPDIR="${TMPDIR:-${VSC_SCRATCH}/apptainer/tmp}"
export NXF_APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/nextflow_cache"

mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$NXF_APPTAINER_CACHEDIR"

export NXF_HOME="${VSC_SCRATCH}/.nextflow"
mkdir -p "$NXF_HOME"

nextflow run "$PIPELINE_DIR" \
  -profile vsc_calcua \
  -params-file "$PARAMS_FILE" \
  -resume

# Optional: single node (all tasks on this allocation; no per-task Slurm jobs). Raise #SBATCH
# resources, then comment the block above and uncomment:
#
# nextflow run "$PIPELINE_DIR" \
#   -profile vsc_calcua,single_node \
#   -params-file "$PARAMS_FILE" \
#   -resume
