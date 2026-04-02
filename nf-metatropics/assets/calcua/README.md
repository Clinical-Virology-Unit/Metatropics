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

| Partition | GPU | Approx. GPU memory (per GPU) | GPUs / node (typical) | Node RAM (typical) | Max walltime |
|-----------|-----|------------------------------|------------------------|-------------------|--------------|
| `ampere_gpu` | NVIDIA A100 | 40 GB | 4 | ~240 GB | 1 day |
| `pascal_gpu` | NVIDIA P100 | 16 GB | 2 | ~112 GB | 1 day |
| `arcturus_gpu` | AMD MI100 | 32 GB | 2 | ~240 GB | 1 day |

Dorado is **NVIDIA/CUDA**: use **`ampere_gpu`** or **`pascal_gpu`**, not `arcturus_gpu`.

**Choosing a GPU queue:** edit [`conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config) and set `calcua_gpu_slurm_partition` (e.g. `pascal_gpu` instead of `ampere_gpu`), or pass at launch:

```bash
nextflow run ... -profile vsc_calcua_gpu --calcua_gpu_slurm_partition pascal_gpu --calcua_gpu_cluster_options '--gres=gpu:1'
```

**How many GPUs:** by default `calcua_gpu_cluster_options = '--gres=gpu:1'`, so each **Dorado** Slurm job requests **one** GPU. The container uses Singularity `--nv`; with a single GPU allocated, that job typically sees one device. To request more GPUs per job (uncommon for Dorado), change the string (e.g. `--gres=gpu:2`) only if your workflow and cluster policy allow it.

**CPU/RAM for basecalling:** `DORADO_ONT` is labelled `process_high` in the pipeline, so on CalcUA it still gets **CPU and memory** from [`conf/vsc_calcua_cpu.config`](../../conf/vsc_calcua_cpu.config) (same labels as other heavy steps). GPU allocation is **only** the Slurm `clusterOptions` above; Dorado’s internal settings (e.g. model in `params_POD5.yaml`, batch size in the module) control how hard the **GPU** works, not Slurm’s CPU count.

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
