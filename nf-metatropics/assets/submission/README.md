# Samplesheets (`fastq.csv` / `POD5.csv`)

This folder holds **example CSV samplesheets** only. Edit paths or barcode names, then point Nextflow at them with `--input` (usually via a params file).

| File | Use when |
|------|----------|
| [`fastq.csv`](fastq.csv) | **FASTQ** already available: `barcode` column = path to each sample’s FASTQ (e.g. `.fastq.gz`). |
| [`POD5.csv`](POD5.csv) | **POD5** runs: `barcode` column = barcode label (e.g. `barcode01`); set `input_dir` in your params to the directory that contains the POD5 data. |

**Parameter files** (paths, `basecall`, resources, etc.) live at the **Metatropics repo root**: [`params_fastq.yaml`](../../../params_fastq.yaml) and [`params_POD5.yaml`](../../../params_POD5.yaml). Copy or edit those; they reference a samplesheet path in `input:` — aim that at your edited `fastq.csv` / `POD5.csv` (or a copy elsewhere).

```bash
nextflow run nf-metatropics/ -profile docker -params-file params_fastq.yaml -resume
```

### HPC (CalcUA / Slurm)

See **[`../calcua/`](../calcua/)** for Slurm scripts and `-profile vsc_calcua` / `vsc_calcua_gpu`.
