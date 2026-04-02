# CalcUA (VSC Antwerp)

To run Metatropics on CalcUA, use Slurm and keep Apptainer caches on `$VSC_SCRATCH` (the script below sets that). Edit `#SBATCH` (account, partition, walltime) and set `PIPELINE_DIR` / `PARAMS_FILE` in [submit_metatropics_calcua.sbatch](submit_metatropics_calcua.sbatch) if your paths differ from the defaults. Then submit the pipeline job; it runs Nextflow with `-profile vsc_calcua` (see [`conf/vsc_calcua.config`](../../conf/vsc_calcua.config)).

Container images are pulled automatically when needed and cached under your scratch Apptainer paths (notably `NXF_APPTAINER_CACHEDIR`, set in the `.sbatch` file), so later runs reuse existing images. The script also prefers the job’s `$TMPDIR` for `APPTAINER_TMPDIR` when available (fallback to scratch) to avoid scratch tmp quota issues during Docker→SIF conversion.

**Run the pipeline** (from the **Metatropics** repo root, the folder that contains `nf-metatropics/`):

```bash
sbatch nf-metatropics/assets/calcua/submit_metatropics_calcua.sbatch
```

If the repo lives elsewhere on scratch, call `sbatch` with the full path to the `.sbatch` file. Copies next to a run directory work too, e.g. `sbatch $VSC_SCRATCH/Metatropics_runs/submit_metatropics_calcua.sbatch`.
