# CalcUA (VSC Antwerp)

To run Metatropics on CalcUA, use Slurm and keep Apptainer caches on `$VSC_SCRATCH` (the scripts below set that). Edit `#SBATCH` (account, partition, walltime) and set `PIPELINE_DIR` / `PARAMS_FILE` in each `.sbatch` if your paths differ from the defaults.

There are **two** submission scripts:

| Script | Input | Params file | Nextflow profile |
|--------|--------|-------------|------------------|
| [submit_metatropics_fastq_calcua.sbatch](submit_metatropics_fastq_calcua.sbatch) | FASTQ | `params_fastq.yaml` | `vsc_calcua` (CPU) |
| [submit_metatropics_pod5_calcua.sbatch](submit_metatropics_pod5_calcua.sbatch) | POD5 + Dorado | `params_POD5.yaml` | `vsc_calcua_gpu` (CPU tasks + GPU queue for Dorado) |

CPU-focused config: [`conf/vsc_calcua_cpu.config`](../../conf/vsc_calcua_cpu.config). GPU extras: [`conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config).

Container images are pulled automatically when needed and cached under your scratch Apptainer paths (notably `NXF_APPTAINER_CACHEDIR`, set in the `.sbatch` file), so later runs reuse existing images. The script also prefers the job’s `$TMPDIR` for `APPTAINER_TMPDIR` when available (fallback to scratch) to avoid scratch tmp quota issues during Docker→SIF conversion.

### GPU partitions (Dorado / CUDA)

These are the usual **GPU** queues on CalcUA (Vaughan / Leibniz). Confirm with `sinfo -o '%P %G %l %m'` and [VSC Antwerp hardware docs](https://docs.vscentrum.be/antwerp/tier2_hardware.html); exact quotas and access depend on your project.

| Partition | GPU | memory | node | Max walltime |
|-----------|-----|--------|------|--------------|
| `ampere_gpu` | NVIDIA A100 | 40 GB | 4 | 1 day |
| `pascal_gpu` | NVIDIA P100 | 16 GB | 2 | 1 day |
| `arcturus_gpu` | AMD MI100 | 32 GB | 2 | 1 day |

Dorado is **NVIDIA/CUDA**: use **`ampere_gpu`** or **`pascal_gpu`**, not `arcturus_gpu`.

To use another GPU queue or Slurm GPU flags, edit `calcua_gpu_slurm_partition` and `calcua_gpu_cluster_options` in [`conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config), then run via [`submit_metatropics_pod5_calcua.sbatch`](submit_metatropics_pod5_calcua.sbatch) (same `nextflow` line is already in that script).

**CPU/RAM for Dorado:** `DORADO_*` use the `process_gpu` label. Default limits are in [`conf/base.config`](../../conf/base.config); CalcUA tuning for POD5 runs is in [`conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config) (together with the GPU partition / `--gres`). GPU allocation is what you set via `calcua_gpu_cluster_options`. The Dorado image in the modules is pulled from Docker Hub and cached as a SIF under `NXF_APPTAINER_CACHEDIR` like other images.

On CalcUA, the `vsc_calcua` profile supports these **CPU** partitions (max per-task resources):

| Partition | Max CPU | Max RAM | Max walltime |
|----------|---------|---------|--------------|
| `zen2` | 64 | 240 GB | 3 days |
| `zen3` | 64 | 240 GB | 3 days |
| `zen3_512` | 64 | 496 GB | 3 days |
| `broadwell` | 28 | 112 GB | 3 days |
| `broadwell_256` | 28 | 240 GB | 3 days |
| `skylake` | 28 | 176 GB | 7 days |

In Slurm-scheduled mode, these are limits for individual pipeline tasks; for `single_node` runs they’re “up to what you requested with `sbatch`” (see commented block in the FASTQ `.sbatch`).

The CalcUA profile also tunes process label resources conservatively so the same pipeline run can work on any of the partitions above; see [`conf/vsc_calcua_cpu.config`](../../conf/vsc_calcua_cpu.config).

**Submit** (from the **Metatropics** repo root, the folder that contains `nf-metatropics/`):

```bash
# FASTQ
sbatch nf-metatropics/assets/calcua/submit_metatropics_fastq_calcua.sbatch

# POD5 + Dorado (GPU)
sbatch nf-metatropics/assets/calcua/submit_metatropics_pod5_calcua.sbatch
```

If the repo lives elsewhere on scratch, use the full path to the `.sbatch` file, e.g. `sbatch $VSC_SCRATCH/Metatropics_runs/submit_metatropics_fastq_calcua.sbatch`.
