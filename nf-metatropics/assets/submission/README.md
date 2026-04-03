# Samplesheets

These files are **example CSV samplesheets**. Each has a header row and one row per sample:

| Column | Meaning |
|--------|---------|
| `sample` | Sample name used in the pipeline. |
| `single_end` | `True` or `False`. |
| `barcode` | Depends on the file (see below). |

**[`fastq.csv`](fastq.csv)** - for runs that start from **FASTQ** reads: `barcode` is the **full path** to that sample’s `.fastq.gz` file.

**[`POD5.csv`](POD5.csv)** — for runs that start from **POD5** data: `barcode` is the **barcode label** (e.g. `barcode01`) for that sample. Which directory holds the POD5 files is **not** listed in the CSV; set that separately in your params next to the samplesheet path.

**[`all_options.md`](all_options.md)** — pipeline parameters as tables.
