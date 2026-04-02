#!/bin/bash -l
#SBATCH --job-name=mt_pod5
#SBATCH --partition=zen3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=72:00:00
#SBATCH --account=ap_itg_cvu
#SBATCH --output=metatropics_pod5_%j.out
#SBATCH --error=metatropics_pod5_%j.err

set -euo pipefail

# Pipeline checkout and params: POD5 + Dorado basecalling (GPU). Params: basecall: true, input_dir, usegpu: true.
PIPELINE_DIR="${VSC_SCRATCH}/Metatropics_runs/Metatropics/nf-metatropics"
PARAMS_FILE="${VSC_SCRATCH}/Metatropics_runs/Metatropics/params_POD5.yaml"

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

# `-profile vsc_calcua_gpu` includes the full CPU profile (conf/vsc_calcua_cpu.config) and only
# sends DORADO_* processes to the GPU partition (ampere_gpu). All later steps (same as FASTQ input)
# submit to the CPU partition inferred from this job, e.g. zen3 via SLURM_JOB_PARTITION.

nextflow run "$PIPELINE_DIR" \
  -profile vsc_calcua_gpu \
  -params-file "$PARAMS_FILE" \
  -resume
