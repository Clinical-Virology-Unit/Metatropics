# CalcUA (VSC Antwerp)

To run Metatropics on CalcUA, use Slurm and keep Apptainer caches on `$VSC_SCRATCH` (the scripts below set that). Edit `#SBATCH` (account, partition, walltime) and set `PIPELINE_DIR` / `PARAMS_FILE` in each `.sbatch` if your paths differ from the defaults.

There are **two** submission scripts:

| Script | Input | Params file | Nextflow profile |
|--------|--------|-------------|------------------|
| [submit_metatropics_fastq_calcua.sbatch](submit_metatropics_fastq_calcua.sbatch) | FASTQ | `params_fastq.yaml` | `vsc_calcua` (CPU) |
| [submit_metatropics_pod5_calcua.sbatch](submit_metatropics_pod5_calcua.sbatch) | POD5 + Dorado | `params_POD5.yaml` | `vsc_calcua_gpu` (CPU tasks + GPU queue for Dorado) |

CPU-focused config: [`conf/vsc_calcua_cpu.config`](../../conf/vsc_calcua_cpu.config). GPU extras: [`conf/vsc_calcua_gpu.config`](../../conf/vsc_calcua_gpu.config).

For GPU steps (Dorado), `calcua_gpu_slurm_partition` / `calcua_gpu_cluster_options` are set in the GPU config (edit or override via CLI if your partition differs).

Container images are pulled automatically when needed and cached under your scratch Apptainer paths (notably `NXF_APPTAINER_CACHEDIR`, set in the `.sbatch` file), so later runs reuse existing images. The script also prefers the job’s `$TMPDIR` for `APPTAINER_TMPDIR` when available (fallback to scratch) to avoid scratch tmp quota issues during Docker→SIF conversion.

On CalcUA, the `vsc_calcua` profile supports these Slurm partitions (max per-task resources):

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
