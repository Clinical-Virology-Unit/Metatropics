#!/bin/bash -l
#SBATCH --job-name=mt_img_dl
#SBATCH --partition=zen3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --account=ap_itg_cvu
#SBATCH --output=initial_images_download_%A_%a.out
#SBATCH --error=initial_images_download_%A_%a.err
# One Slurm array task per image (18 tasks → indices 0–17). Tasks may run on different nodes.
# To limit concurrent pulls (e.g. Docker Hub): use --array=0-17%8
#SBATCH --array=0-17

set -euo pipefail

# docker:// URIs (singularity_pull_docker_container mode). Keep array length in sync with #SBATCH --array above.
METATROPICS_DOCKER_IMAGES=(
  "docker://daanjansen94/python3:v3"
  "docker://daanjansen94/bbmap:38.86"
  "docker://daanjansen94/fastp:v0.23.1"
  "docker://daanjansen94/nanoplot:v1.41.0"
  "docker://daanjansen94/minimap:v2.28"
  "docker://daanjansen94/samtools:v1.17"
  "docker://daanjansen94/metamaps:v0.1"
  "docker://daanjansen94/multiqc:v1.12"
  "docker://daanjansen94/ivar:v1.4.1"
  "docker://daanjansen94/krona:v2.7.1"
  "docker://daanjansen94/seqtk:v1.3"
  "docker://daanjansen94/medaka:v1.4.3"
  "docker://daanjansen94/nf_r_plots:v4.2.2"
  "docker://daanjansen94/rcoverage:v4.0.2"
  "docker://daanjansen94/bam-readcount:v0.8"
  "docker://daanjansen94/homopolish:v0.4.1"
  "docker://daanjansen94/rocker-tidyverse:latest"
  "docker://nanoporetech/dorado:sha4644018526d3644d92a9e680ab8f2d1eeff2e272"
)

unset JAVA_HOME || true
unset JAVA_CMD || true

module load Java
module load Nextflow

export APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/cache"
export APPTAINER_TMPDIR="${VSC_SCRATCH}/apptainer/tmp"
export NXF_APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/nextflow_cache"

mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$NXF_APPTAINER_CACHEDIR"

IDX="${SLURM_ARRAY_TASK_ID:?}"
URI="${METATROPICS_DOCKER_IMAGES[$IDX]}"

REF="${URI#docker://}"
SAFE=$(echo "$REF" | tr '/:' '-')
OUT="${NXF_APPTAINER_CACHEDIR}/${SAFE}.img"

if [[ -s "$OUT" ]]; then
  echo "[${IDX}] Already present, skipping: $OUT"
  exit 0
fi

echo "[${IDX}] Pulling -> $OUT"
singularity pull "$OUT" "$URI"
echo "[${IDX}] Done."
