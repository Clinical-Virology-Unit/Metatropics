# CalcUA (VSC Antwerp)

This guide is for **running Metatropics on CalcUA**, the Tier‑1 HPC at **VSC Antwerp**. On a cluster you normally **submit a job** with a script instead of starting Nextflow by hand on a login node; the `.sbatch` files here do that and point Metatropics at storage and containers in a way that fits CalcUA. Under the hood they use Slurm (the scheduler), Apptainer (containers on HPC), and the Nextflow profiles `vsc_calcua` or `vsc_calcua_gpu`.

**Before you submit:** clone the Metatropics repository and set up `params_fastq.yaml` or `params_POD5.yaml` as in the **repository root README**.

## How to run

1. **Edit the `.sbatch` you need** — set `#SBATCH` (account, partition, walltime), `PIPELINE_DIR`, and `PARAMS_FILE`
2. **Submit from the Metatropics repository root**

| Your input | Script | Params file (typical) | Nextflow profile |
|------------|--------|------------------------|------------------|
| FASTQ | [submit_metatropics_fastq_calcua.sbatch](submit_metatropics_fastq_calcua.sbatch) | `params_fastq.yaml` | `vsc_calcua` (CPU) |
| POD5 | [submit_metatropics_pod5_calcua.sbatch](submit_metatropics_pod5_calcua.sbatch) | `params_POD5.yaml` | `vsc_calcua_gpu` (GPU + CPU) |

```bash
# FASTQ
sbatch nf-metatropics/assets/calcua/submit_metatropics_fastq_calcua.sbatch

# POD5 + Dorado (GPU)
sbatch nf-metatropics/assets/calcua/submit_metatropics_pod5_calcua.sbatch
```

Configs used by those profiles: [`../../conf/vsc_calcua_cpu.config`](../../conf/vsc_calcua_cpu.config), [`../../conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config).

---

## Additional details

Optional reading: partitions, caches, and where to tune resources.

**Apptainer caches** — Images are pulled when needed and cached under scratch (`NXF_APPTAINER_CACHEDIR` in the `.sbatch`). The script prefers the job’s `$TMPDIR` for `APPTAINER_TMPDIR` when present (else scratch) to reduce tmp-quota issues during Docker→SIF conversion.

**GPU partitions (Dorado / CUDA)** - Usual GPU queues on CalcUA (Vaughan / Leibniz). Confirm with `sinfo -o '%P %G %l %m'` and the [VSC Antwerp hardware docs](https://docs.vscentrum.be/antwerp/tier2_hardware.html); access depends on your project.

| Partition | GPU | GPU memory | GPUs/node | Max walltime |
|-----------|-----|------------|-----------|--------------|
| `ampere_gpu` | NVIDIA A100 | 40 GB | 4 | 1 day |
| `pascal_gpu` | NVIDIA P100 | 16 GB | 2 | 1 day |
| `arcturus_gpu` | AMD MI100 | 32 GB | 2 | 1 day |

Dorado is **NVIDIA/CUDA**: use **`ampere_gpu`** or **`pascal_gpu`**, not `arcturus_gpu`. To change GPU queue or Slurm GPU flags, edit `calcua_gpu_slurm_partition` and `calcua_gpu_cluster_options` in [`../../conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config), then submit with [`submit_metatropics_pod5_calcua.sbatch`](submit_metatropics_pod5_calcua.sbatch).

**CPU/RAM for Dorado** — `DORADO_*` processes use the `process_gpu` label. Defaults in [`../../conf/base.config`](../../conf/base.config); CalcUA POD5 tuning in [`../../conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config) (with GPU partition / `--gres`). GPU allocation follows `calcua_gpu_cluster_options`. The Dorado image is pulled from Docker Hub and cached as a SIF under `NXF_APPTAINER_CACHEDIR` like other images.

**CPU partitions** — With profile `vsc_calcua`, these CPU partitions are supported (max per-task resources):

| Partition | Max CPU | Max RAM | Max walltime |
|-----------|---------|---------|--------------|
| `zen2` | 64 | 240 GB | 3 days |
| `zen3` | 64 | 240 GB | 3 days |
| `zen3_512` | 64 | 496 GB | 3 days |
| `broadwell` | 28 | 112 GB | 3 days |
| `broadwell_256` | 28 | 240 GB | 3 days |
| `skylake` | 28 | 176 GB | 7 days |

In Slurm-scheduled mode these are limits for individual pipeline tasks; for `single_node` runs they are effectively up to what you requested in `sbatch` (see comments in the FASTQ `.sbatch`). Process labels are tuned conservatively in [`../../conf/vsc_calcua_cpu.config`](../../conf/vsc_calcua_cpu.config) so the same run can work across these partitions.
