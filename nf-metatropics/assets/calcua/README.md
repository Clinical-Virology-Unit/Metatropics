# CalcUA (VSC Antwerp)

To run Metatropics on CalcUA, use Slurm and keep Apptainer caches on `$VSC_SCRATCH` (the scripts below set that). **Step 1:** edit `#SBATCH` (account, partition, walltime) in both [initial_images_download_calcua.sbatch](initial_images_download_calcua.sbatch) and [submit_metatropics_calcua.sbatch](submit_metatropics_calcua.sbatch), and set `PIPELINE_DIR` / `PARAMS_FILE` in [submit_metatropics_calcua.sbatch](submit_metatropics_calcua.sbatch) if your paths differ from the defaults. **Step 2:** submit [initial_images_download_calcua.sbatch](initial_images_download_calcua.sbatch) once and wait until every array task has finished (check `squeue` / log files)—this pulls the pipeline’s Docker images into your cache so later runs are faster. **Step 3:** submit [submit_metatropics_calcua.sbatch](submit_metatropics_calcua.sbatch); it runs Nextflow with `-profile vsc_calcua` (see [`conf/vsc_calcua.config`](../../conf/vsc_calcua.config)) and reuses those images.

Images are stored under your scratch Apptainer paths (notably `NXF_APPTAINER_CACHEDIR`, set in the `.sbatch` files). **Later runs** use the same variables, so existing images are picked up automatically and are not downloaded again. You can skip Step 2 and rely on Nextflow to pull any missing image during the pipeline, but the first run will be slower and may hit pull timeouts; Step 2 avoids that by filling the cache up front.

**Step 2 — initial container images** (from the **Metatropics** repo root, the folder that contains `nf-metatropics/`):

```bash
sbatch nf-metatropics/assets/calcua/initial_images_download_calcua.sbatch
```

**Step 3 — run the pipeline:**

```bash
sbatch nf-metatropics/assets/calcua/submit_metatropics_calcua.sbatch
```

If the repo lives elsewhere on scratch, call `sbatch` with the full path to each `.sbatch` file. Copies next to a run directory work too, e.g. `sbatch $VSC_SCRATCH/Metatropics_runs/initial_images_download_calcua.sbatch`.
