#!/bin/bash -l
#SBATCH --job-name=metatropics
#SBATCH --partition=zen3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=72:00:00
#SBATCH --account=ap_itg_cvu
#SBATCH --output=metatropics_%j.out
#SBATCH --error=metatropics_%j.err

set -euo pipefail

# Edit these if your Metatropics checkout is not under this path on scratch.
PIPELINE_DIR="${VSC_SCRATCH}/Metatropics_runs/Metatropics/nf-metatropics"
PARAMS_FILE="${VSC_SCRATCH}/Metatropics_runs/Metatropics/params.yaml"

unset JAVA_HOME || true
unset JAVA_CMD || true

module load Java
module load Nextflow

export APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/cache"
export APPTAINER_TMPDIR="${VSC_SCRATCH}/apptainer/tmp"
export NXF_APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/nextflow_cache"

mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$NXF_APPTAINER_CACHEDIR"

export NXF_HOME="${VSC_SCRATCH}/.nextflow"
mkdir -p "$NXF_HOME"

# Default: Slurm schedules each pipeline task as its own job; keep this head job small.
nextflow run "$PIPELINE_DIR" \
  -profile vsc_calcua \
  -params-file "$PARAMS_FILE" \
  -resume

# Optional: run on a single node (all processes on that allocation; no per-task Slurm jobs).
# Raise #SBATCH cpus/mem/time above to cover the whole run (e.g. a full node). Then comment
# out the nextflow block above and uncomment:
#
# nextflow run "$PIPELINE_DIR" \
#   -profile vsc_calcua,single_node \
#   -params-file "$PARAMS_FILE" \
#   -resume
