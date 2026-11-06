# Samplesheets

These files are **example CSV samplesheets**. Each has a header row and one row per sample:

| Column | Meaning |
|--------|---------|
| `sample` | Sample name used in the pipeline. |
| `barcode` | Depends on the file (see below). |

**[`fastq.csv`](fastq.csv)** — for runs that start from **FASTQ** reads: `barcode` is the **full path** to that sample’s read file. To **auto-generate** a samplesheet CSV from a folder of FASTQ files, follow the **`--input`** note in **[`all_options.md`](all_options.md)**.

**[`POD5.csv`](POD5.csv)** — use this samplesheet when starting from **raw signal** (`pod5`) data, or from **basecalled but not demultiplexed** reads (`fastq_pass`). In both cases, `barcode` is the **barcode label** (e.g. `barcode01`) for that sample. The data directory is set separately in your params file (`input_dir`), not in the CSV.

**[`wells.csv`](wells.csv)** — TWIST lab plate (`sample,well`). Generate with **`metatropics-samplesheet pod5 TWIST-96A-UDI`** (or `fastq_pass` / `demultiplexed`). Usable **alone** as `input` when `kit_name` matches the plate (well → `barcodeNN` is mapped automatically), or convert to `POD5.csv` / `fastq_pass.csv` with **`metatropics-samplesheet … TWIST-96A-UDI wells.csv`**.

**[`all_options.md`](all_options.md)** — pipeline parameters as tables.
